#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>
#include "header.h"

int rank;
hsize_t dims[2], count[2], start[2];

hid_t hdf5_file, hdf5_headergrp, hdf5_grp;
hid_t hdf5_attribute; 
hid_t hdf5_datatype, hdf5_dataset, hdf5_dataspace_in_file, hdf5_dataspace_in_memory;
hid_t hdf5_status;

void read_hdf5_header(char *filename, struct sim_info *header, long long *NumPart)
{
  int i,j;
  long long NumPartInFile=0;
  
  //  sprintf(filename,"%s.hdf5",temp);

  SnapFormat=3;
  
  hdf5_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5_headergrp = H5Gopen(hdf5_file, "/Header");

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
  hdf5_status=H5Aread(hdf5_attribute, H5T_NATIVE_INT, header->npart);
  hdf5_status=H5Aclose(hdf5_attribute);
  
  for(i=0;i<6;i++) {
    NumPartInFile+=header->npart[i];
    fprintf(stdout,"Species %u : n_in_file %d\n",i,header->npart[i]);
    fflush(stdout);
  }
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
  hdf5_status=H5Aread(hdf5_attribute, H5T_NATIVE_INT, &(header->NumFiles));
  hdf5_status=H5Aclose(hdf5_attribute);


  if(header->NumFiles>1) {
    hdf5_status=H5Gclose(hdf5_headergrp);
    hdf5_status=H5Fclose(hdf5_file);
    
    for(i=0;i<header->NumFiles;i++)
      {
	//	sprintf(filename,"%s%d.hdf5",temp,i);
        
	hdf5_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	hdf5_headergrp = H5Gopen(hdf5_file, "/Header");
        
	hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
	hdf5_status=H5Aread(hdf5_attribute, H5T_NATIVE_INT, header->npart);
	hdf5_status=H5Aclose(hdf5_attribute);

	
	for(j=0;j<6;j++) {
	  *NumPart+=header->npart[j];
	}
      }
    
    hdf5_status=H5Gclose(hdf5_headergrp);
    hdf5_status=H5Fclose(hdf5_file);
    
    //    sprintf(filename,"%s",argv[1]);
    
    hdf5_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    hdf5_headergrp = H5Gopen(hdf5_file, "/Header");
  } else 
    *NumPart = NumPartInFile;
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total");
  hdf5_status=H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header->nall);
  hdf5_status=H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "MassTable");
  hdf5_status=H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header->massarr);
  hdf5_status=H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Time");
  hdf5_status=H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &(header->time));
  hdf5_status=H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "BoxSize");
  hdf5_status=H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &(header->BoxSize));
  hdf5_status=H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Omega0");
  hdf5_status=H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &(header->Omega0));
  hdf5_status=H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "OmegaLambda");
  hdf5_status=H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &(header->OmegaLambda));
  hdf5_status=H5Aclose(hdf5_attribute);
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "HubbleParam");
  hdf5_status=H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &(header->HubbleParam));
  hdf5_status=H5Aclose(hdf5_attribute);

  hdf5_status=H5Gclose(hdf5_headergrp);
  hdf5_status=H5Fclose(hdf5_file);
}

void read_gadget_binary_header(char *filename, struct sim_info *header, long long *NumPart)
{
  int i,j;
  long long NumPartInFile;
  
  FILE *infile;
  char label[4];
  long long fileoffset;
  
  // First check what SnapFormat the file is; SnapFormat=1 implies that there is a small
  // buffer of size int that gives the number of particles in the next block, while
  // SnapFormat=2 implies that there are label blocks that describe the contents of the next
  // block, e.g. "POS" for positions, "VEL" for velocities, etc...

  infile=fopen(filename,"r");
  fseek(infile,sizeof(int),0);  // Skip the first block - an single int
  fread(&label, sizeof(char), 4, infile);
  SnapFormat=1;
  if(strcmp(label,"HEAD")==0) SnapFormat=2;
  rewind(infile);
  
  if(SnapFormat==2)
    {
      fileoffset=sizeof(int)+4*sizeof(char)+2*sizeof(int)+sizeof(int);
    }
  else
    {
      fileoffset=sizeof(int);
    }
  
  fseek(infile,fileoffset,0);
  
  fread(header->npart,sizeof(int),6,infile);

  for(i=0;i<6;i++) {
    NumPartInFile+=header->npart[i];
  }
  
  fileoffset+=6*sizeof(int);
  
  fseek(infile,fileoffset,0);
  
  fread(&(header->massarr),sizeof(double),6,infile);
  
  fread(&(header->time),sizeof(double),1,infile);
  
  fileoffset+=6*sizeof(double)+2*sizeof(double)+2*sizeof(int);
  
  fseek(infile,fileoffset,0);
  fread(header->nall,sizeof(int),6,infile);
  
  for(i=0;i<6;i++) {
    fprintf(stdout,"Species %d : n_in_file %d, n_in_snapshot %d\n",i,header->npart[i],header->nall[i]);
    fflush(stdout);
    *NumPart+=header->nall[i];
  }
  
  fileoffset+=7*sizeof(int);
  
  fseek(infile,fileoffset,0);
  fread(&(header->NumFiles),sizeof(int),1,infile);
  fread(&(header->BoxSize),sizeof(double),1,infile);
  fread(&(header->Omega0),sizeof(double),1,infile);
  fread(&(header->OmegaLambda),sizeof(double),1,infile);
  fread(&(header->HubbleParam),sizeof(double),1,infile);
  
  fileoffset+=sizeof(int)+4*sizeof(double)+9*sizeof(int)+(256-196)+sizeof(int);
  
  fclose(infile);
}

void read_particles_from_hdf5(char *temp, float *x, float *y, float *z, int *ptype, int NumFiles, long long *NThisTask)
{
  int i,j,k;
  int nfile;
  float pos[3];
  int nx,ny,nz;
  int npart[6];
  char filename[132],buffer[32];
  float *dummy;

  for(nfile=0;nfile<NumFiles;nfile++)
    {
      if(NumFiles>1)
	{
	  sprintf(filename,"%s%d.hdf5",temp,nfile);
	} else {
	sprintf(filename,"%shdf5",temp);
      }

      hdf5_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

      hdf5_headergrp = H5Gopen(hdf5_file, "/Header");
      hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
      hdf5_status=H5Aread(hdf5_attribute, H5T_NATIVE_INT, npart);
      hdf5_status=H5Aclose(hdf5_attribute);
      hdf5_status=H5Gclose(hdf5_headergrp);

      for(i = 0; i < 6; i++)
	{
	  if(npart[i] > 0)
	    {
	      sprintf(buffer, "/PartType%d", i);
	      hdf5_grp = H5Gopen(hdf5_file, buffer);
	      hdf5_dataset = H5Dopen(hdf5_grp, "Coordinates");
	      hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
              
	      rank = 2;

	      dims[1] = 3;              
	      dims[0] = npart[i];

	      hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
              
	      dims[0] = npart[i];
	      hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);
	      
	      start[0]=0;
	      start[1]=0;
              
	      count[0]= npart[i];
	      count[1]= 3;
	      
	      hdf5_status=H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,
				  start, NULL, count, NULL);
	      
	      dummy = (float*)malloc(3*sizeof(float)*count[0]);
              
	      hdf5_status=H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory,
		      hdf5_dataspace_in_file, H5P_DEFAULT, dummy);
	      
	      k=0;
	      for(j=0;j<count[0];j++)
		{
		  x[*NThisTask]=dummy[k];
		  y[*NThisTask]=dummy[k+1];
		  z[*NThisTask]=dummy[k+2];
		  ptype[*NThisTask]=i;
		  k+=3;
                  (*NThisTask)++;
		}
	      free(dummy);
              
	      hdf5_status=H5Sclose(hdf5_dataspace_in_memory);
	      hdf5_status=H5Sclose(hdf5_dataspace_in_file);
	      hdf5_status=H5Tclose(hdf5_datatype);
	      hdf5_status=H5Dclose(hdf5_dataset);
	      hdf5_status=H5Gclose(hdf5_grp);	      
	    }
	}
      
      hdf5_status=H5Fclose(hdf5_file);

    }

}

void read_particles_from_gadget_binary(char *temp, float *x, float *y, float *z, int *ptype, int NumFiles, long long *NThisTask)
{
  int i,j,k;
  int nfile;
  float pos[3];
  int nx,ny,nz;
  int npart[6];
  char filename[32],buffer[32];
  float *dummy;

  long long NumPartInFile;
  long long fileoffset;
  FILE *infile;

  *NThisTask=0;
  
  for(nfile=0;nfile<NumFiles;nfile++)
    {  
      if(NumFiles>1)
	sprintf(filename,"%s.%d",temp,nfile);
      else
	sprintf(filename,"%s",temp);
      
      infile=fopen(filename,"r");
      
      // Now we need to read in the data...
      if(SnapFormat==2)
	{
	  fileoffset=sizeof(int)+4*sizeof(char)+2*sizeof(int)+sizeof(int);
	}
      else
	{
	  fileoffset=sizeof(int);
	}
      
      fseek(infile,fileoffset,0);
      
      NumPartInFile=0;
      
      fread(&npart,sizeof(int),6,infile);
      for(i=0;i<6;i++) {
	NumPartInFile+=npart[i];
      }
      
      fileoffset+=(256+sizeof(int));
      
      if(SnapFormat==2)
	{
	  fileoffset=+sizeof(int)+4*sizeof(char)+2*sizeof(int)+sizeof(int);
	}
      else
	{
	  fileoffset+=sizeof(int);
	}
      
      fseek(infile,fileoffset,0);
      
      fprintf(stdout,"Reading %lld particles from %s on Task %d...\n",NumPartInFile,filename,ThisTask);
      fflush(stdout);
      
      fileoffset+=(3*ThisTask*(NumPartInFile/NTask)*sizeof(float));
      
      fseek(infile,fileoffset,0);
      
      for(i=ThisTask*(NumPartInFile/NTask);i<(ThisTask+1)*(NumPartInFile/NTask);i++)
	{
	  fread(&pos,3*sizeof(float),1,infile);
	  x[*NThisTask]=pos[0];
	  y[*NThisTask]=pos[1];
	  z[*NThisTask]=pos[2];
	  ptype[*NThisTask]=0;
	  (*NThisTask)++;
	}
      fclose(infile);
    }
}
