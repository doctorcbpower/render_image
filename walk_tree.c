#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "omp.h"

#include "header.h"

int cmpfunc(const void* elem1, const void* elem2)
{
  if(*(const float*)elem1 < *(const float*)elem2)
    return -1;
  return *(const float*)elem1 > *(const float*)elem2;
}

void build_interaction_list(struct tree_node_nd *node, struct interaction_list **list, int *index)
{
  struct interaction_list *tmp;
  
  if(node!=NULL)
    {
      build_interaction_list(node->p0,&(*list),index);

      tmp=(struct interaction_list *)malloc(sizeof(struct interaction_list));
      tmp->index=(*index)++;
      tmp->node=node;
      tmp->left=NULL;
      tmp->right=NULL;

      if(*list==NULL)
	{
	  (*list)=tmp;
	}
      else
	{
	  (*list)->left=tmp;
	  tmp->right=(*list);
	  (*list)=tmp;
	}
      
      build_interaction_list(node->p1,&(*list),index);
#if NDIM==2 || NDIM==3
      build_interaction_list(node->p2,&(*list),index);
      build_interaction_list(node->p3,&(*list),index);
#endif
#if NDIM==3
      build_interaction_list(node->p4,&(*list),index);
      build_interaction_list(node->p5,&(*list),index);
      build_interaction_list(node->p6,&(*list),index);
      build_interaction_list(node->p7,&(*list),index);
#endif
    }      
}

void get_distance_to_nth_nearest_neighbour(struct point *P, int num_pts, int num_ngb, struct tree_node_nd *root)
{
  int i,j;
  int count, ncount;
  float x[NDIM];
  float search, fsearch;
  float dx2;
  float r2[100000];
  
  struct link_list *ngb, *tmp;
  //struct btree_node *broot;
  /*
  struct interaction_list *int_list, *tmp_int;
  int id_int=0;
  */
  // Estimate mean interparticle separation using bounding box
  search=1.0;
  for(i=0;i<NDIM;i++)
    search*=(root->xmax[i]-root->xmin[i]);
  
  search=pow(search,1./(float)NDIM);
  
  search*=pow((float)num_pts,-1./(float)NDIM);      
#if NDIM==2
  search*=M_PI/4.;
#endif
#if NDIM==3
  search*=M_PI/6.;
#endif

  printf("Using a search radius of %g units...\n",search);
  /*
  build_interaction_list(root,&int_list,&id_int);
  printf("Built interaction list...\n");
  */

#ifdef ENABLE_OPENMP  
  int NThreads=4;
  int ThisThread;

#pragma omp parallel num_threads(NThreads) private(x,r2,count,ncount,fsearch,ngb,tmp,i,j,ThisThread)
  {
    ThisThread=omp_get_thread_num();    

    for(i=ThisThread*num_pts/NThreads;i<(ThisThread+1)*num_pts/NThreads;i++)
#else
    for(i=0;i<num_pts;i++)      
#endif
      {
#if NDIM==3
	memcpy(x,&(P[i].pos[0]),3*sizeof(float));
#elif NDIM==2
	memcpy(x,&(P[i].pos[0]),2*sizeof(float));
#else
	memcpy(x,&(P[i].pos[0]),sizeof(float));      
#endif
	fsearch=0.1;
	count=0;

	while(count<=num_ngb)
	  {
	    count=0;
	    ncount=0;
	    ngb=NULL;
	    get_multiple_nodes_nd(root,x,fsearch*search,&count,&ncount,&ngb);
	    if(count>num_ngb) break;
	    fsearch*=1.5;
	  }

	tmp=ngb;
	count=0;

	while(tmp!=NULL)
	  {
	    r2[count]=0.0;
	    for(j=0;j<NDIM;j++)
	      r2[count]+=(x[j]-tmp->x[j])*(x[j]-tmp->x[j]);
	    count++;
	    tmp=tmp->ptr;
	  }

	qsort(r2,count,sizeof(float),cmpfunc);

	P[i].dist_ngb=sqrt(r2[num_ngb]);
      }
#ifdef ENABLE_OPENMP  
#pragma omp barrier
  }
#endif

  return;
}

void get_kernel_density_estimate(struct point *P, int num_pts, struct tree_node_nd *root)
{
  int i,j;
  int count,ncount;
  float x[NDIM];
  float dx2, xvar;
  float hfac;
  
  struct link_list *ngb, *tmp;

  double density;
  float h;
  
  for(i=0;i<num_pts;i++)
    {
#if NDIM==3
      memcpy(x,&(P[i].pos[0]),3*sizeof(float));      
#elif NDIM==2
      memcpy(x,&(P[i].pos[0]),2*sizeof(float));
#else
      memcpy(x,&(P[i].pos[0]),sizeof(float));      
#endif
      
      h=P[i].dist_ngb;
      hfac=1./h;
#if NDIM==2
      hfac/=h;
#elif NDIM==3
      hfac/=pow(h,2.);
#endif      
      count=0;
      ncount=0;

      ngb=NULL;
      get_multiple_nodes_nd(root,x,h,&count,&ncount,&ngb);
      
      tmp=ngb;
      density=0.0;

      while(tmp!=NULL)
	{
	  dx2=0.0;
	  for(j=0;j<NDIM;j++)
	    dx2+=(x[j]-tmp->x[j])*(x[j]-tmp->x[j]);
	  xvar=sqrt(dx2)/h;

	  density+=P[i].mass*get_kernel_value(xvar)*hfac;

	  tmp=tmp->ptr;
	}
      P[i].density=density;

      }
  
  return;
}

void get_potential_estimate(struct point *P, int num_pts, struct tree_node_nd *root)
{
  int i,j;
  int count,ncount;
  float x[NDIM];
  float theta=0.3;
  
  for(i=0;i<1;i++)
    {
#if NDIM==3
      memcpy(x,&(P[i].pos[0]),3*sizeof(float));      
#elif NDIM==2
      memcpy(x,&(P[i].pos[0]),2*sizeof(float));
#else
      memcpy(x,&(P[i].pos[0]),sizeof(float));      
#endif
      count=0;
      ncount=0;
      sweep_over_nodes_nd(root,x,theta,&count,&ncount);
      printf("%d %d\n",count,ncount);
    }
  
  return;
}

void scan_nodes(struct tree_node *node, int *count)
{

  if(node!=NULL)
    {
      scan_nodes(node->left, count);
      (*count)++;
      scan_nodes(node->right, count);
    }

}

void get_node(struct tree_node *node, float x, int *count, struct link_list **ngb)
{
  float xmin,xmax,xmid;
  
  int i;
  struct link_list *tmp;
  
  if(node!=NULL)
    {
      xmin=node->xmin;
      xmax=node->xmax;
      
      xmid=0.5*(xmin+xmax);
      
      if(x<=xmid)
	if(node->left!=NULL)
	  get_node(node->left,x,count,&(*ngb));
	else
	  {
	    (*count)+=node->num_members;
	    for(i=0;i<node->num_members;i++)
	      {
		tmp=(struct link_list *)malloc(sizeof(struct link_list));
		tmp->x[0]=node->x[i];
		if(*ngb==NULL)
		  tmp->ptr=NULL;
		else
		  tmp->ptr=*ngb;
		*ngb=tmp;
	      }
	  }
      else
	if(node->right!=NULL)	
	  get_node(node->right,x,count,&(*ngb));
	else
	  {
	    *count+=node->num_members;	  
	    for(i=0;i<node->num_members;i++)
	      {
		tmp=(struct link_list *)malloc(sizeof(struct link_list));
		tmp->x[0]=node->x[i];
		if(*ngb==NULL)
		  tmp->ptr=NULL;
		else
		  tmp->ptr=*ngb;
		*ngb=tmp;
	      }
	  }
    }
  
}

void get_multiple_nodes(struct tree_node *node, float x, float search, int *count, struct link_list **ngb)
{
  float xmin,xmax;
  int i;
  struct link_list *tmp;
  
  if(node!=NULL)
    {
      xmin=node->xmin;
      xmax=node->xmax;
      
      if((x-search)<xmax && (x+search)>xmin)
	{
	  get_multiple_nodes(node->left,x,search,count,&(*ngb));
	  
	  if(node->split==0)
	    {
	      *count+=node->num_members;
	      for(i=0;i<node->num_members;i++)
		{
		  tmp=(struct link_list *)malloc(sizeof(struct link_list));
		  tmp->x[0]=node->x[i];
		  if(*ngb==NULL)
		    tmp->ptr=NULL;
		  else
		    tmp->ptr=*ngb;
		  *ngb=tmp;
		}
	    }
	  get_multiple_nodes(node->right,x,search,count,&(*ngb));
	}
    }      
}

void add_node_index_nd(struct tree_node_nd *node, int *count)
{
  if(node!=NULL)
    {
      node->index=++(*count);
      
      add_node_index_nd(node->p0,count);
      add_node_index_nd(node->p1,count);

#if NDIM==2 || NDIM==3
      add_node_index_nd(node->p2,count);
      add_node_index_nd(node->p3,count);
#endif      

#if NDIM==3
      add_node_index_nd(node->p4,count);
      add_node_index_nd(node->p5,count);
      add_node_index_nd(node->p6,count);
      add_node_index_nd(node->p7,count);
#endif      
    }
}

void get_node_nd(struct tree_node_nd *node, float *x, int *count, struct link_list **ngb)
{
  int i,j,k;
  
  float xmin[NDIM],xmax[NDIM];
  float xmid[NDIM];
  
  struct link_list *tmp;
  
  if(node!=NULL)
    {
      for(i=0;i<NDIM;i++)
	xmid[i]=0.5*(node->xmin[i]+node->xmax[i]);
      
#ifdef DEBUG
      for(i=0;i<NDIM;i++)
	printf("(xmin, xmax)[%d]: (%g,%g) x[%d]: %g, N: %d\n",i,node->xmin[i],node->xmax[i],i,x[i],node->num_members);
      printf("\n");
#endif
      
      j=1,k=0;
      if(x[0]<=xmid[0]) j=0;
      k+=j;
#if NDIM==2 || NDIM==3      
      j=1;
      if(x[1]<=xmid[1]) j=0;
      k+=2*j;
#endif
#if NDIM==3      
      j=1;
      if(x[2]<=xmid[2]) j=0;
      k+=4*j;
#endif
      int open_node=1;
      
      if(k==0)
	{
	  if(node->p0!=NULL)
	    {
	      get_node_nd(node->p0,x,count,&(*ngb));
	      open_node=0;
	    }
	}
      else if(k==1)
	{
	  if(node->p1!=NULL)
	    {
	      get_node_nd(node->p1,x,count,&(*ngb));
	      open_node=0;
	    }
	}
#if NDIM==2 || NDIM==3
      else if(k==2)
	{
	  if(node->p2!=NULL)
	    {
	      get_node_nd(node->p2,x,count,&(*ngb));
	      open_node=0;
	    }
	}
      else if(k==3)
	{
	  if(node->p3!=NULL)
	    {
	      get_node_nd(node->p3,x,count,&(*ngb));
	      open_node=0;
	    }
	}
#endif
#if NDIM==3
      else if(k==4)
	{
	  if(node->p4!=NULL)
	    {
	      get_node_nd(node->p4,x,count,&(*ngb));
	      open_node=0;
	    }
	}
      else if(k==5)
	{
	  if(node->p5!=NULL)
	    {
	      get_node_nd(node->p5,x,count,&(*ngb));
	      open_node=0;
	    }
	}
      else if(k==6)
	{
	  if(node->p6!=NULL)
	    {
	      get_node_nd(node->p6,x,count,&(*ngb));
	      open_node=0;
	    }
	}
      else if(k==7)
	{
	  if(node->p7!=NULL)
	    {
	      get_node_nd(node->p7,x,count,&(*ngb));
	      open_node=0;
	    }
	}
#endif      
      if(open_node==1)
	{
	  (*count)+=node->num_members;
	  for(i=0;i<node->num_members;i++)
	    {
	      tmp=(struct link_list *)malloc(sizeof(struct link_list));
	      for(j=0;j<NDIM;j++)
		tmp->x[j]=node->x[j][i];
	      if(*ngb==NULL)
		tmp->ptr=NULL;
	      else
		tmp->ptr=*ngb;
	      *ngb=tmp;
	    }
	}
    }
}

void get_multiple_nodes_nd(struct tree_node_nd *node, float *x, float search, int *count, int *ncount, struct link_list **ngb)
{
  int i,j;
  float xmin[NDIM],xmax[NDIM];
  struct link_list *tmp;
  
  if(node!=NULL)
    {
      for(i=0;i<NDIM;i++)
	{
	  xmin[i]=node->xmin[i];
	  xmax[i]=node->xmax[i];
	}

      if((x[0]-search)<xmax[0] && (x[0]+search)>xmin[0]
#if NDIM==2 || NDIM==3
	 && (x[1]-search)<xmax[1] && (x[1]+search)>xmin[1]
#endif
#if NDIM==3
	 && (x[2]-search)<xmax[2] && (x[2]+search)>xmin[2]
#endif
	 )
	{
	  get_multiple_nodes_nd(node->p0,x,search,count,ncount,&(*ngb));      
	  
	  if(node->split==0)
	    {
	      (*ncount)++;
	      *count+=node->num_members;
	      for(i=0;i<node->num_members;i++)
		{
		  tmp=(struct link_list *)malloc(sizeof(struct link_list));
		  for(j=0;j<NDIM;j++)
		    tmp->x[j]=node->x[j][i];
		  if(*ngb==NULL)
		    tmp->ptr=NULL;
		  else
		    tmp->ptr=*ngb;
		  *ngb=tmp;
		}
	    }
	  get_multiple_nodes_nd(node->p1,x,search,count,ncount,&(*ngb));      
#if NDIM==2 || NDIM==3
	  get_multiple_nodes_nd(node->p2,x,search,count,ncount,&(*ngb));
	  get_multiple_nodes_nd(node->p3,x,search,count,ncount,&(*ngb));      	  
#endif      
#if NDIM==3
	  get_multiple_nodes_nd(node->p4,x,search,count,ncount,&(*ngb));
	  get_multiple_nodes_nd(node->p5,x,search,count,ncount,&(*ngb));
	  get_multiple_nodes_nd(node->p6,x,search,count,ncount,&(*ngb));
	  get_multiple_nodes_nd(node->p7,x,search,count,ncount,&(*ngb));
#endif      
	}
    }
}

void sweep_over_nodes_nd(struct tree_node_nd *node, float *x, float theta, int *count, int *ncount)
{
  int i,j;
  float xmin[NDIM],xmid[NDIM],xmax[NDIM];
  float xmean[NDIM];
  float dist_to_cell,cell_size,cell_offset;
  struct link_list *tmp;

  if(node!=NULL)
    {

      cell_size=1.0;
      cell_offset=0.0;
      
      for(i=0;i<NDIM;i++)
	{
	  xmin[i]=node->xmin[i];
	  xmax[i]=node->xmax[i];
	  xmid[i]=0.5*(node->xmin[i]+node->xmax[i]);
	  xmean[i]=node->xmean[i];
	  cell_size*=(xmax[i]-xmin[i]);
	  cell_offset+=(xmean[i]-xmid[i])*(xmean[i]-xmid[i]);
	}

      cell_size=pow(cell_size,1./(float)NDIM);
      cell_offset=sqrt(cell_offset);
      
      dist_to_cell=0.0;
      
      for(j=0;j<NDIM;j++)
	dist_to_cell+=(x[j]-xmean[j])*(x[j]-xmean[j]);

      if(dist_to_cell<(cell_size/theta+cell_offset))
	{
	  sweep_over_nodes_nd(node->p0,x,theta,count,ncount);      
	  printf("%g %g %g %g %d\n",x[0],x[1],xmid[0],xmid[1],node->num_members);
	  
	  *count+=node->num_members;

	  /*      
	      if(node->split==0)
	      {
	      (*ncount)++;

	      for(i=0;i<node->num_members;i++)
	      {
	      tmp=(struct link_list *)malloc(sizeof(struct link_list));
	      for(j=0;j<NDIM;j++)
	      tmp->x[j]=node->x[j][i];
	      if(*ngb==NULL)
	      tmp->ptr=NULL;
	      else
	      tmp->ptr=*ngb;
	      *ngb=tmp;
	      }
	      }
      */
      
	  sweep_over_nodes_nd(node->p1,x,theta,count,ncount);      
#if NDIM==2 || NDIM==3
	  sweep_over_nodes_nd(node->p2,x,theta,count,ncount);
	  sweep_over_nodes_nd(node->p3,x,theta,count,ncount);      	  
#endif      
#if NDIM==3
	  sweep_over_nodes_nd(node->p4,x,theta,count,ncount);
	  sweep_over_nodes_nd(node->p5,x,theta,count,ncount);
	  sweep_over_nodes_nd(node->p6,x,theta,count,ncount);
	  sweep_over_nodes_nd(node->p7,x,theta,count,ncount);
#endif
	}
    }

  
}


