#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

#include "header.h"

void find_neighbours(int num_part,float *smoothing_length,int num_ngb,
		     float *posx, float *posy, float *posz,
		     float xmin, float ymin, float zmin)
{
  int i;
  float *r2;
  struct tree_node_nd *root_nd=NULL, *node_nd;
  struct point *P;

  clock_t start=clock();
  
  fprintf(stdout,"Initialising search mesh for %d particles... %d\n",num_part,ThisTask);
  fflush(stdout);

  P=(struct point *)malloc(num_part*sizeof(struct point));  

  for(i=0;i<num_part;i++)
    {
      P[i].pos[0]=posx[i]-xmin;
      P[i].pos[1]=posy[i]-ymin;
      P[i].pos[2]=posz[i]-zmin;
      P[i].mass=0.1;
    }

  make_tree(P,num_part,&root_nd);
  printf("Finished making tree...\n");
 
  clock_t finish=clock()-start;

  fprintf(stdout,"Time (%lu s)...\n",finish/CLOCKS_PER_SEC); 
  fflush(stdout);

  start=finish;
  
  get_distance_to_nth_nearest_neighbour(P,num_part,num_ngb,root_nd);
  printf("Finished walking tree...\n");

  for(i=0;i<num_part;i++)
    smoothing_length[i]=P[i].dist_ngb;

  finish=clock()-start;

  fprintf(stdout,"Time (%lu s)...\n",finish/CLOCKS_PER_SEC); 
  fflush(stdout);

  start=finish;

  return;
  
}
 
