#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"

int main(int argc, char *argv[])
{
  int i,j,k,l,m,n;
  char filename[128];
  FILE *infile,*outfile;

  int isHDF5=0;

  // Now check that user input is correct
  if(argc<2)
    {
      fprintf(stdout,"Usage: %s <filename> (<isHDF5>)\n",argv[0]);
      fflush(stdout);
      exit(0);
    }
  
  sprintf(filename,"%s",argv[1]);
  
  if(argc>2)
    if(strcmp(argv[2],"-isHDF5")==0)
      isHDF5=1;
  
  if(isHDF5==1){
    fprintf(stdout,"Assumung HDF5 input...\n");
    fflush(stdout);
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

  // Read in data from file...
  
  fprintf(stdout,"Reading %s\n",filename);
  fflush(stdout);

  if(isHDF5==1)
    read_header_hdf5(filename);
  else
    read_header_bin(filename);

  NumPart=0;
  
  for(i=0;i<6;i++)
    {
      NumPart+=nall[i];
      fprintf(stdout,"Species %d : n_in_file %ld\n",i,npart[i]);
      fflush(stdout);
    }

  // Allocate memory
  fprintf(stdout,"Allocating memory for %ld particles...\n",NumPart);
  fflush(stdout);  

  x = (float*)malloc(sizeof(float)*(NumPart));
  y = (float*)malloc(sizeof(float)*(NumPart));
  z = (float*)malloc(sizeof(float)*(NumPart));
  h = (float*)malloc(sizeof(float)*(NumPart));  
  id = (int*)malloc(sizeof(int)*(NumPart));  
  
  if(isHDF5==1)
    {
      read_block_hdf5(filename,"Coordinates");  // Positions
      read_block_hdf5(filename,"ParticleIDs");  // IDs
    }
  else
    {
      read_block_bin(filename,1);   // Positions
      read_block_bin(filename,3);   // IDs
    }

  fprintf(stdout,"Carrying out neighbour search on %d particles...\n",NumPart);
  fflush(stdout);    
  find_neighbours(512);
  scatter_to_mesh(64);
  // Finish up...
  
  return(0);
  
}


