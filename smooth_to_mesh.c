#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "omp.h"

#include "header.h"

#define LOG10DENSMAX 13.5
#define LOG10DENSMIN 0

void smooth_to_mesh(long long NumPart, float *smoothing_length, int num_ngb,
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

#ifdef KERNEL_SMOOTHING    
  tabulate_integral();
#endif
  
#ifdef ENABLE_OPENMP  
  int NThreads=4;
  int ThisThread;
  
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

      delx=cos(theta)*dx-sin(theta)*dz;
      dely=dy;
      delz=sin(theta)*dx+cos(theta)*dz;

      nx=(int)width*(delx/lbox+0.5);
      ny=(int)height*(dely/lbox+0.5);

      if(nx<0 || nx>=width || ny<0 || ny>=height) continue;

#ifdef KERNEL_SMOOTHING  
      dh=height*smoothing_length[i]/lbox;
      
      h_inv=1./smoothing_length[i];
      dh_inv=1./dh;
      
      if(dh>=20) dh=20;
      
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
  	    //if(jz<0) jz+=height;
  	    //if(jz>=width) jz-=height;

	    if(jx<0 || jx>=width || jy<0 || jy>=width) continue;
	    
	    j=jx+width*jy;

	    fmin=get_integral_value(drmin*dh_inv);
	    fmax=get_integral_value(drmax*dh_inv);
	    
  	    dens[j]+=(fmax-fmin)*h_inv*h_inv*h_inv;	    
	  }
#else
      dx=width*(delx/lbox+0.5)-(nx+0.5);
      dy=height*(dely/lbox+0.5)-(ny+0.5);

      if(dx>0)
	jx=nx+1;
      else
	jx=nx-1;

      if(dy>0)
	jy=ny+1;
      else
	jy=ny-1;

      if(jx<0 || jx>=width || jy<0 || jy>=width) continue;
      
      hx0=1.-fabs(dx);
      hx1=fabs(dx);
      
      hy0=1.-fabs(dy);
      hy1=fabs(dy);

      j=nx+width*ny;
      dens[j]+=hx0*hy0*width*height/pow(lbox,2);
      j=jx+width*ny;
      dens[j]+=hx1*hy0*width*height/pow(lbox,2);
      j=nx+width*jy;
      dens[j]+=hx0*hy1*width*height/pow(lbox,2);
      j=jx+width*jy;
      dens[j]+=hx1*hy1*width*height/pow(lbox,2);
#endif      
    }
#ifdef ENABLE_OPENMP  
#pragma omp barrier
  }
#endif    

  for(i=0; i<height*width; i++)
    {
      if(dens[i]>dmax) dmax=dens[i];
      if(dens[i]<dmin) dmin=dens[i];
    }

  printf("Log10 dmin:%g Log10 dmax:%g\n",log10(dmin),log10(dmax));
  for(i=0; i<height*width; i++)
    {
      if(log10(dens[i])<LOG10DENSMIN) dens[i]=pow(10,LOG10DENSMIN);
      dens[i]=(log10(dens[i])-LOG10DENSMIN)/(LOG10DENSMAX-LOG10DENSMIN);
    }

    for(j=0;j<height;j++)
    for(i=0;i<width;i++)
      {
	*((data+j*width)+i)=dens[i+j*width];
      }  
}

