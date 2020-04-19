#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <hdf5.h>

#ifdef DOUBLE_FFTW
#include <drfftw_mpi.h>
#else
#include <rfftw_mpi.h>
#endif

#include "header.h"

#define BUFSIZE 500

#define rhocrit 27.755

#include <png.h>

void find_neighbours(int num_part,float *smoothing_length, int num_ngb,
		     float *posx, float *posy, float *posz,
		     float xmin, float ymin, float zmin);


float finterpolant(float xmin, float xmax);

/*
void setRGB(png_byte *ptr, float val);

void setRGB(png_byte *ptr, float val)
{
  int v = (int)(val * 768);
  if (v < 0) v = 0;
  if (v > 768) v = 768;
  int offset = v % 256;
  
  if (v<256) {
    ptr[0] = 0; ptr[1] = 0; ptr[2] = 0; ptr[3]=offset;
  }
  else if (v>256&&v<512) {
    ptr[0] = 0; ptr[1]=0; ptr[2] = offset; ptr[3] = 255-offset;
  }
  else {
    ptr[0] = offset; ptr[1] = 0; ptr[2]= 255-offset; ptr[3] = offset;
  }
}
*/
#define IMAGE_DIMENSION 1024
#define MAX_COLOUR_COMPONENT_VALUE 255
#define FMIN 2.5e-10

static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;
static int fftsize, maxfftsize;
static fftw_complex *fft_of_rhogrid;

struct linking_list {
    int num;
    struct linking_list *ptr;
};

typedef struct linking_list LIST;

int main(int argc, char *argv[])
{
  int i,j,k,l,m,n,buf[BUFSIZE];
  char filename[128],image_file[128],image_file_root[128],debug_file[128];
  FILE *infile,*outfile;
  char buffer[BUFSIZE];
  double SliceSize,ImageSliceSize,f_Image=0.1;
  double mp;

  int nsend,nrec;

  float *send, *rec;
  int *isend, *irec;


  char junk[256-24];
  long long NumPart=0, NumPartInFile=0;
  long long fileoffset;

  float pos[3];
  float *x,*y,*z;
  int *ptype;
  float *posx,*posy,*posz;
  int *partid;
  float *smoothing_length;
  float xmin,ymin,zmin;
  float xmax,ymax,zmax;
  double sum,sumx,sumy,sumz;

  int PMGRID;
	
  float fslice=0.1;
  
  int NSlab;
  int Nmesh;
  int nx,ny,nz;
  int ix,iy,iz;
  int jx,jy,jz;
  int kx,ky,kz;  
  long long num_x[BUFSIZE]={0},num_gather[BUFSIZE]={0};
  long long num_proc=0, num_tot=0;

  double dx,dy,dz;
    
  int isHDF5=0;
  int rank;
  float *dummy;

  float hx0,hy0,hz0,hx1,hy1,hz1;
  int lcount;
  hsize_t dims[2],start[2],count[2];
  hid_t hdf5_file, hdf5_headergrp, hdf5_grp, hdf5_attribute, hdf5_datatype, hdf5_dataset, hdf5_dataspace_in_file, hdf5_dataspace_in_memory;

  struct sim_info header;

  long long NumPartRead=0;

  int ptype_keep;

  float BoundingBox;

  int itmax;
  
  ThisTask=0;
  NTask=0;
  
  // Now check that user input is correct
  if(argc<2)
    {
      fprintf(stdout,"Usage: %s <filename> (<isHDF5>)\n",argv[0]); fflush(stdout);
      exit(0);
    }
  
  i=1;

  double xc=0.0,yc=0.0,zc=0.0,lbox=0.0;

  while(i<argc)
    {
      if(strcmp(argv[i],"-input")==0)
	sprintf(filename,"%s",argv[++i]);
      
      if(strcmp(argv[i],"-output")==0)
	sprintf(image_file_root,"%s",argv[++i]);      
      
      if(strcmp(argv[i],"-xc")==0)
	xc=atof(argv[++i]);

      if(strcmp(argv[i],"-yc")==0)
	yc=atof(argv[++i]);

      if(strcmp(argv[i],"-zc")==0)
	zc=atof(argv[++i]);

      if(strcmp(argv[i],"-lbox")==0)
	lbox=atof(argv[++i]);

      if(strcmp(argv[i],"-itmax")==0)
	itmax=atoi(argv[++i]);      

      if(strcmp(argv[i],"-isHDF5")==0) isHDF5=1;
      if(strcmp(argv[i],"-stars")==0) ptype_keep=5;
      if(strcmp(argv[i],"-gas")==0) ptype_keep=0;
      if(strcmp(argv[i],"-dark_matter")==0) ptype_keep=1;      

      i++;
    }

  if(isHDF5==1){
    fprintf(stdout,"Assuming HDF5 input...\n");
    fflush(stdout);
    sprintf(filename,"%s.hdf5",filename);
  }

  if((infile=fopen(filename,"r"))==0)
    {
      sprintf(filename,"%s.0",argv[1]);
      if((infile=fopen(filename,"r"))==0)
        {
	  fprintf(stdout,"Error: cannot find input file as either %s or %s\n ",argv[1],filename);
	  fflush(stdout);
	  exit(0);
        }
    }
  
  fprintf(stdout,"Filename: %s\n",filename);
  
  if(xc>0 || lbox>0)
    {
      fprintf(stdout,"Centre: (%g|%g|%g) \n",xc,yc,zc);
      fprintf(stdout,"Box Length: %g\n",lbox);
    }
	  
  
  // File exists - now let's get some basic information. For simplicity, do this on Task 0.

  char *var = (char*)filename;
  char temp[1024]={};
  
  i=0;
  
  while (*var)
    {
      temp[i++]=*var;
      if(strncmp(var,".",1)==0)
        {
	  break;
        } else {
	var++;
      }
    }

  if(isHDF5==1) 
    read_hdf5_header(filename, &header, &NumPart);
  else
    read_gadget_binary_header(filename, &header, &NumPart);

  if(header.BoxSize==0) header.BoxSize=1.e6;
  
  // Echo the cosmological parameters to stdout
  fprintf(stdout,"Omega0 %g OmegaLambda %g HubbleParam %g\n",header.Omega0,header.OmegaLambda,header.HubbleParam);
  fprintf(stdout,"NumPart %lld NumFiles %d BoxSize %g time %g mp %g\n",NumPart,header.NumFiles,header.BoxSize,header.time,mp);
  fflush(stdout);

  // Allocate memory....

  x = (float*)malloc(sizeof(float)*NumPart);
  y = (float*)malloc(sizeof(float)*NumPart);
  z = (float*)malloc(sizeof(float)*NumPart);
  ptype = (int*)malloc(sizeof(int)*NumPart);  

  if(isHDF5==1)
    read_particles_from_hdf5(temp,x,y,z,ptype,header.NumFiles,&NumPartRead);    
  else
    read_particles_from_gadget_binary(temp,x,y,z,ptype,header.NumFiles,&NumPartRead);    

  fprintf(stdout,"Read particles...\n");
  fflush(stdout);
  
  select_particles(x,y,z,ptype,header.BoxSize,NumPartRead,xc,yc,zc,lbox,ptype_keep,&NumPart);
  
  //  split_across_tasks_as_slabs(x,y,z,pid,NumPartRead,BoxSize);

  smoothing_length = (float*)malloc(2*sizeof(float)*(NumPart));  

  float dxmin[3]={0,0,0};
  int num_ngb=40;

  clock_t tstart=clock();  
  find_neighbours(NumPart,smoothing_length,num_ngb,x,y,z,dxmin[0],dxmin[1],dxmin[2]);
  clock_t tfinish=clock()-tstart;

  fprintf(stdout,"find_neighbours - time: %lu s...\n",tfinish/CLOCKS_PER_SEC); 
  fflush(stdout);
  
  float data[IMAGE_DIMENSION][IMAGE_DIMENSION];

  BoundingBox=lbox;
  float theta=0.;
  int iter=0;
  
  while(BoundingBox>0.01*lbox)
    {
      tstart=tfinish;
      smooth_to_mesh(NumPart,smoothing_length,num_ngb,x,y,z,xc,yc,zc,BoundingBox,theta,IMAGE_DIMENSION,IMAGE_DIMENSION,*data);
      tfinish=clock()-tstart;    
      fprintf(stdout,"smooth_to_mesh - time: %lu s...\n",tfinish/CLOCKS_PER_SEC); 
      fflush(stdout);        
      
      sprintf(image_file,"%s.%d.%s",image_file_root,iter,"png");  
      
      //write_to_ppm(image_file,IMAGE_DIMENSION,IMAGE_DIMENSION,MAX_COLOUR_COMPONENT_VALUE,*data);
      write_to_png(image_file,IMAGE_DIMENSION,IMAGE_DIMENSION,*data);      
      BoundingBox/=1.025;
      theta+=2*M_PI/72.0;
      iter++;
      if(iter>itmax-1) break;
    }


  /*

*/  
  fprintf(stdout,"Finished...\n");
  fflush(stdout);
}
