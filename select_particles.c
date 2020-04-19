#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "header.h"

void select_particles(float *x, float *y, float *z, int *ptype, double BoxSize, long long NumPart, float xc, float yc, float zc, float lbox, int ptype_keep, long long *NumPartKeep)
{
  int i, j;
  int ix,iy,iz;
  
  float dx[3];
  float dxmin[3]={1.e10,1.e10,1.e10},dxmax[3]={-1e10,-1e10,-1e10};

  for(i=0,*NumPartKeep=0;i<NumPart;i++)
    {
      dx[0]=x[i]-xc;
      dx[1]=y[i]-yc;
      dx[2]=z[i]-zc;

      for(j=0;j<3;j++)
	{
	  if(dx[j]>0.5*BoxSize) dx[j]-=BoxSize;
	  if(dx[j]<-0.5*BoxSize) dx[j]+=BoxSize;      
	}

      ix=(int)(2*dx[0]/lbox);
      iy=(int)(2*dx[1]/lbox);
      iz=(int)(2*dx[2]/lbox);      

      if(ix==0 && iy==0 && iz==0 && ptype[i]==ptype_keep)
	{
	  x[*NumPartKeep]=x[i];
	  y[*NumPartKeep]=y[i];
	  z[*NumPartKeep]=z[i];
	  if(x[i]<dxmin[0]) dxmin[0]=x[i];
	  if(y[i]<dxmin[1]) dxmin[1]=y[i];
	  if(z[i]<dxmin[2]) dxmin[2]=z[i];	  
	  if(x[i]>dxmax[0]) dxmax[0]=x[i];
	  if(y[i]>dxmax[1]) dxmax[1]=y[i];
	  if(z[i]>dxmax[2]) dxmax[2]=z[i];	  	  
	  (*NumPartKeep)++;
	}
    }

  /*
  for(j=0;j<3;j++)
    xc[j]/=(float)NumPart;

  xc[0]+=x[0];
  xc[1]+=y[0];
  xc[2]+=z[0];  
  
  fprintf(stdout,"Weighted centre: %g %g %g (Task %d)\n",xc[0],xc[1],xc[2],ThisTask);
  fflush(stdout);

  // Shift particles so that weighted centre lies at box centre
  
  for(i=0;i<3;i++)
    {
      dxmin[i]=1.e10;
      dxmax[i]=-1.e10;
    }
  
  for(i=0;i<NumPart;i++)
    {
      dx[0]=x[i]-x[0];
      dx[1]=y[i]-y[0];
      dx[2]=z[i]-z[0];

      for(j=0;j<3;j++)
	{
	  if(dx[j]>0.5*BoxSize) dx[j]-=BoxSize;
	  if(dx[j]<-0.5*BoxSize) dx[j]+=BoxSize;      
	  if(dx[j]<dxmin[j]) dxmin[j]=dx[j];
	  if(dx[j]>dxmax[j]) dxmax[j]=dx[j];	  	  
	}      
    }

  for(j=0;j<3;j++)
    {
      xc[j]+=0.5*(dxmin[j]+dxmax[j]);
      if(xc[j]>BoxSize) xc[j]-=BoxSize;
      if(xc[j]<0) xc[j]+=BoxSize;
    }

  fprintf(stdout,"Bounding box centre: %g %g %g (Task %d)\n",xc[0],xc[1],xc[2],ThisTask);
  fflush(stdout);  
  
  float LocalBoundingBoxLength=-1.e10;
  
  for(i=0;i<3;i++)
    if((dxmax[i]-dxmin[i])>LocalBoundingBoxLength) LocalBoundingBoxLength=dxmax[i]-dxmin[i];

  LocalBoundingBoxLength*=1.1;

  MPI_Allreduce(&LocalBoundingBoxLength,&(*BoundingBoxLength),1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
  fprintf(stdout,"Bounding Box Length: %g (Task %d)\n",*BoundingBoxLength,ThisTask);
  fflush(stdout);
  
  for(i=0;i<NumPart;i++)
    {
      if(ptype[i]>2) continue;
      
  */
}
