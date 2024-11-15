#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

#include "header.h"

void smooth_to_mesh(long long NumPart, float *smoothing_length, 
		    float *x, float *y, float *z,
		    float xc, float yc, float zc, float lbox,
		    float theta,
		    int width, int height, float *data)
{
  int i,j,k,l,m;
  int ix,iy,iz;
  int jx,jy,jz;
  int nx,ny,nz;
  
  float dx,dy,dz;
  float hx0,hx1,hy0,hy1;
  float dxmin,dxmax,dymin,dymax,drmin,drmax;

  double dmax=0.0,dmin=1.e10;  
  double *dens;

  int dn;
  float dh, dh_inv;

  float h_inv;

  float delx,dely,delz;
  float fmin, fmax;

  dens = (double*) malloc(sizeof(double)*width*height);

  for(i=0;i<width*height;i++)
    dens[i]=0;
    
    printf("%d %d\n",width,height);

#ifdef KERNEL_SMOOTHING    
  tabulate_integral();
#endif
  
#ifdef ENABLE_OPENMP  
  int NThreads=4;
  int ThisThread;
  
    printf("Number of threads: %d\n",omp_get_num_threads());
    printf("Number of processors: %d\n",omp_get_num_procs());

#pragma omp parallel num_threads(NThreads) private(i,ix,iy,iz,j,jx,jy,jz,nx,ny,nz,dn,dx,dy,dz,dxmin,dxmax,dymin,dymax,drmin,drmax,dh,h_inv,dh_inv,k,l,m,hx0,hx1,hy0,hy1,delx,dely,delz,fmin,fmax,ThisThread)
  {
    ThisThread=omp_get_thread_num();    
    
    for(i=ThisThread*NumPart/NThreads;i<(ThisThread+1)*NumPart/NThreads;i++)
#else  
  for(i=0; i<NumPart; i++)
#endif
    {
      dx=x[i]-xc;
      dy=y[i]-yc;
      dz=z[i]-zc;

//      delx=cos(theta)*dx-sin(theta)*dz;
//      dely=dy;
//      delz=sin(theta)*dx+cos(theta)*dz;
      
      nx=(int)width*(dx/lbox+0.5);
      ny=(int)height*(dy/lbox+0.5);

#ifdef NONPERIODIC      
      if(nx<0 || nx>=width || ny<0 || ny>=height) continue;
#else
      if(nx<0) nx+=width;
      if(nx>=width) nx-=width;
      if(ny<0) ny+=height;
      if(ny>=height) ny-=height;            
#endif
      
#ifdef KERNEL_SMOOTHING  
      dh=width*smoothing_length[i]/lbox;
      
      h_inv=1./smoothing_length[i];
      dh_inv=1./dh;
      
      if(dh>=30) dh=30;
      
      dn=1+(int)dh;
      
      for(iy=ny-dn;iy<ny+dn;iy++)
      	for(ix=nx-dn;ix<nx+dn;ix++)
	  {
	    dx=width*(delx/lbox+0.5)-(ix+0.5);
	    dy=height*(dely/lbox+0.5)-(iy+0.5);
	    
  	    if(dx>0)
  	      {
  		dxmin=dx-0.5;
  		dxmax=dx+0.5;
  	      }
  	    else
  	      {
  		dxmin=dx+0.5;
  		dxmax=dx-0.5;
  	      }
	    
  	    if(dy>0)
  	      {
  		dymin=dy-0.5;
  		dymax=dy+0.5;
  	      }
  	    else
  	      {
  		dymin=dy+0.5;
  		dymax=dy-0.5;
  	      }
	    
  	    drmin=sqrt(dxmin*dxmin+dymin*dymin);
  	    drmax=sqrt(dxmax*dxmax+dymax*dymax);
	    
  	    jx=ix;
  	    jy=iy;

#ifdef NONPERIODIC      
	    if(jx<0 || jx>=width || jy<0 || jy>=height) continue;
#else
	    if(jx<0) jx+=width;
	    if(jx>=width) jx-=width;
	    if(jy<0) jy+=height;
	    if(jy>=height) jy-=height;      
#endif
	    j=jx+width*jy;

	    fmin=get_integral_value(drmin*dh_inv);
	    fmax=get_integral_value(drmax*dh_inv);
	    
  	    dens[j]+=(fmax-fmin)*h_inv*h_inv*h_inv;	    
	  }
#else
      dx=width*(dx/lbox+0.5)-(nx+0.5);
      dy=height*(dy/lbox+0.5)-(ny+0.5);

      if(dx>0)
	jx=nx+1;
      else
	jx=nx-1;

      if(dy>0)
	jy=ny+1;
      else
	jy=ny-1;

#ifdef NONPERIODIC      
      if(jx<0 || jx>=width || jy<0 || jy>=height) continue;
#else
      if(jx<0) jx+=width;
      if(jx>=width) jx-=width;
      if(jy<0) jy+=height;
      if(jy>=height) jy-=height;      
#endif
      
      hx0=1.-fabs(dx);
      hx1=fabs(dx);
      
      hy0=1.-fabs(dy);
      hy1=fabs(dy);
        
      j=nx+width*ny;
      dens[j]+=hx0*hy0*width*width/pow(lbox,2);
      j=jx+width*ny;
      dens[j]+=hx1*hy0*width*width/pow(lbox,2);
      j=nx+width*jy;
      dens[j]+=hx0*hy1*width*width/pow(lbox,2);
      j=jx+width*jy;
      dens[j]+=hx1*hy1*width*width/pow(lbox,2);

#endif
    }
#ifdef ENABLE_OPENMP  
#pragma omp barrier
  }
#endif    

  for(i=0; i<height*width; i++)
    {
        printf("%lf\n",dens[i]);

      if(dens[i]>dmax) dmax=dens[i];
      if(dens[i]<dmin) dmin=dens[i];
    }
    
    if(dmin==0.0) dmin=1.e0;
  
  printf("Log10 dmin:%g Log10 dmax:%g\n",log10(dmin),log10(dmax));
  
  for(j=0;j<height;j++)
    for(i=0;i<width;i++)
      {
	*((data+j*width)+i)=dens[i+j*width];
      }  
}

