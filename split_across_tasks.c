#include <stdio.h>
#include <stdlib.h>
#include "header.h"
#include "mpi.h"

void split_across_tasks_as_slabs(float *x, float *y, float *z, int **pid, long long NumPart, float BoxSize)
{
  int i,j;
  int nx;
  long long num_tot, num_proc, *num_x, *num_gather;
  /*  
  num_x=(long long *)malloc(sizeof(long long)*NTask);
  num_gather=(long long *)malloc(sizeof(long long)*NTask);
  
  for(i=0;i<NTask;i++)
    num_x[i]=0;
  
  for(i=0;i<NumPart;i++)
    {
      nx=(int)(NTask*(x[i]/BoxSize));   //Note that (int) always rounds down
#ifdef CIC
      if(nx<0) nx+=NTask;
      if(nx>=NTask) nx-=NTask;
#else
      if(nx<0) nx=0;
      if(nx>=NTask) nx=NTask-1;
#endif      
      pid[nx][num_x[nx]]=i;
      num_x[nx]++;
    }

  for(i=0;i<NTask;i++)
    num_gather[i]=0;
  
  for(i=0;i<NTask;i++)
    {
      if(i!=ThisTask)
        {
	  MPI_Send(&num_x[i],1,MPI_LONG_LONG_INT,i,ThisTask,MPI_COMM_WORLD);
	  MPI_Recv(&num_gather[i],1,MPI_LONG_LONG_INT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
  
  num_proc=0;

  for(i=0;i<NTask;i++)
    if(i!=ThisTask) num_proc+=num_gather[i];
  
  num_proc+=num_x[ThisTask];
  num_gather[ThisTask]=num_x[ThisTask];
  
  MPI_Reduce(&num_proc,&num_tot,1,MPI_LONG_LONG_INT,MPI_SUM,0,MPI_COMM_WORLD);
  
  MPI_Bcast(&num_tot,1,MPI_LONG_LONG_INT,0,MPI_COMM_WORLD);
  
  fprintf(stdout,"Allocating memory for %lld particles on PE %d\n\n",num_proc,ThisTask);
  fflush(stdout);  

  float *posx, *posy, *posz;
  int *partid;
  
  posx = (float*)malloc(sizeof(float)*num_proc);
  posy = (float*)malloc(sizeof(float)*num_proc);
  posz = (float*)malloc(sizeof(float)*num_proc);
  partid = (int*)malloc(sizeof(int)*num_proc);  

  long long noffset=0;

  for(j=0;j<ThisTask;j++)
    noffset+=num_gather[j];
  
  for(j=0;j<num_x[ThisTask];j++)
    {
      posx[j+noffset]=x[pid[ThisTask][j]];
      posy[j+noffset]=y[pid[ThisTask][j]];
      posz[j+noffset]=z[pid[ThisTask][j]];
      partid[j+noffset]=pid[ThisTask][j];      
    }

  int nsend,nrec, *isend, *irec;
  float *send, *rec;

  for(i=0;i<NTask;i++)
    {
      if(i!=ThisTask)
        {
	  nsend=num_x[i];
	  nrec=num_gather[i];
          
	  send = (float*)malloc(sizeof(float)*nsend);
	  rec = (float*)malloc(sizeof(float)*nrec);
	  	  
	  noffset=0;
	  for(j=0;j<i;j++)
	    noffset+=num_gather[j];
	  
	  // First do x
	  for(j=0;j<nsend;j++)
            send[j]=x[pid[i][j]];
	  
	  MPI_Sendrecv(send,nsend,MPI_FLOAT,i,ThisTask,
		       rec,nrec,MPI_FLOAT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  
	  for(j=0;j<nrec;j++)
	    posx[j+noffset]=rec[j];
	  
	  // ... then do y
	  for(j=0;j<nsend;j++)
	    send[j]=y[pid[i][j]];
	  
	  MPI_Sendrecv(send,nsend,MPI_FLOAT,i,ThisTask,
		       rec,nrec,MPI_FLOAT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  
	  for(j=0;j<nrec;j++)
	    posy[j+noffset]=rec[j];
	  
	  // ... then do z
	  for(j=0;j<nsend;j++)
	    send[j]=z[pid[i][j]];
	  
	  MPI_Sendrecv(send,nsend,MPI_FLOAT,i,ThisTask,
		       rec,nrec,MPI_FLOAT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  
	  for(j=0;j<nrec;j++)
	    posz[j+noffset]=rec[j];

	  free(send);
	  free(rec);

	  isend = (int*)malloc(sizeof(int)*nsend);
	  irec = (int*)malloc(sizeof(int)*nrec);
	  
	  // ... then do IDs
	  for(j=0;j<nsend;j++)
	    isend[j]=pid[i][j];
	  
	  MPI_Sendrecv(isend,nsend,MPI_INT,i,ThisTask,
		       irec,nrec,MPI_INT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  
	  for(j=0;j<nrec;j++)
	    partid[j+noffset]=irec[j];

	  free(isend);
	  free(irec);
        }
      
    }

  float sumx,sumy,sumz;
  float xmin,ymin,zmin;
  float xmax,ymax,zmax;  

  sumx=sumy=sumz=0.0;
  xmin=ymin=zmin=1.e10;
  xmax=ymax=zmax=-1.e10;
	
  for(i=0;i<num_proc;i++)
    {
      sumx+=posx[i];
      sumy+=posy[i];
      sumz+=posz[i];
        
      if(posx[i]>xmax) xmax=posx[i];
      if(posy[i]>ymax) ymax=posy[i];
      if(posz[i]>zmax) zmax=posz[i];
      if(posx[i]<xmin) xmin=posx[i];
      if(posy[i]<ymin) ymin=posy[i];
      if(posz[i]<zmin) zmin=posz[i];

    }
    
  sumx=sumx/(double)num_proc;
  sumy=sumy/(double)num_proc;
  sumz=sumz/(double)num_proc;

  fprintf(stdout,"Number of particles on Task %d: %lld \n",ThisTask,num_proc);  
  fprintf(stdout,"Average position on Task %d: (%f,%f,%f) \n",ThisTask,sumx,sumy,sumz);
  fprintf(stdout,"Boundaries on Task %d: (%g,%g) \n",ThisTask,xmin,xmax);
  fflush(stdout);

  free(x);
  free(y);
  free(z);
  
  free(num_x);
  */

}
