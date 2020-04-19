#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "header.h"

void make_tree(struct point *P, int num_pts, struct tree_node_nd **root)
{
  int i, j;
  int count=0;
  float x[NDIM];
  float xmin[NDIM],xmax[NDIM];

  struct tree_node_nd *root_nd;

  for(j=0;j<NDIM;j++)
    {
      xmin[j]=1.e10;
      xmax[j]=-1.e10;
    }
  
  for(i=0;i<num_pts;i++)
    for(j=0;j<NDIM;j++)
      {
	if(xmin[j]>P[i].pos[j]) xmin[j]=P[i].pos[j];
	if(xmax[j]<P[i].pos[j]) xmax[j]=P[i].pos[j];
      }
  
  for(j=0;j<NDIM;j++)
      printf("Extent: xmin[%d]: %g, xmax[%d]: %g\n",j,xmin[j],j,xmax[j]);
  
  root_nd=NULL;
  for(i=0;i<num_pts;i++)
    {
#if NDIM==3
      memcpy(x,&(P[i].pos[0]),3*sizeof(float));
#elif NDIM==2
      memcpy(x,&(P[i].pos[0]),2*sizeof(float));
#else
      memcpy(x,&(P[i].pos[0]),sizeof(float));
#endif      
      root_nd=add_to_node_nd(root_nd,x,xmin,xmax);      
    }      
  
  *root=root_nd;
  
  return;
}

struct tree_node *add_to_node(struct tree_node *node, float x, float xmin, float xmax) 
{
  int i;
  float xmid;

  xmid=0.5*(xmin+xmax);
  
  if(node==NULL)
    {
      node=(struct tree_node *)malloc(sizeof(struct tree_node));
      node->num_members=1;
      node->x[node->num_members-1]=x;
      node->xmin=xmin;
      node->xmax=xmax;
      node->split=0;
      node->left=node->right=NULL;
    }
  else
    {
      if(node->split==0)
	{
	  (node->num_members)++;

	  if(node->num_members>MAXNODE)
	    {
	      node->split=1;
	      for(i=0;i<MAXNODE;i++)
		{
		  if(node->x[i]<=xmid)
		    node->left=add_to_node(node->left,node->x[i],xmin,xmid);
		  else
		    node->right=add_to_node(node->right,node->x[i],xmid,xmax);	
		}
	      if(x<=xmid)
		node->left=add_to_node(node->left,x,xmin,xmid);
	      else
		node->right=add_to_node(node->right,x,xmid,xmax);	
	    }
	  else
	    node->x[node->num_members-1]=x;
	}
      else
	{
	  node->num_members++;
	  
	  if(x<=xmid)
	    node->left=add_to_node(node->left,x,xmin,xmid);
	  else
	    node->right=add_to_node(node->right,x,xmid,xmax);	
	  
	}
    }
  
  return node;
}

struct tree_node_nd *add_to_node_nd(struct tree_node_nd *node, float *x, float *xmin, float *xmax) 
{
  int i,j,k;
  float xmid[NDIM];
  float y[NDIM];
  float node_xmin[NDIM],node_xmax[NDIM];
  
  for(i=0;i<NDIM;i++)
    xmid[i]=0.5*(xmin[i]+xmax[i]);
#ifdef VERBOSE  
  for(i=0;i<NDIM;i++)
    printf("(xmin, xmax)[%d]: (%g,%g) x[%d]: %g\n",i,xmin[i],xmax[i],i,x[i]);
  printf("\n");
#endif
  if(node==NULL)
    {
      node=(struct tree_node_nd *)malloc(sizeof(struct tree_node_nd));
      node->num_members=1;
      for(i=0;i<NDIM;i++)
	{
	  node->x[i][node->num_members-1]=x[i];
	  node->xmin[i]=xmin[i];
	  node->xmax[i]=xmax[i];
	  node->xmean[i]=x[i];
	}
      node->split=0;
      node->p0=node->p1=NULL;
#if NDIM==2 || NDIM==3      
      node->p2=node->p3=NULL;
#endif
#if NDIM==3            
      node->p4=node->p5=node->p6=node->p7=NULL;
#endif
    }
  else  // node==NULL
    {
      if(node->split==0)
	{
	  (node->num_members)++;

	  if(node->num_members>MAXNODE)
	    {
	      // First process existing particles in cell
	      node->split=1;	  
	      for(i=0;i<MAXNODE;i++)
		{
		  for(j=0;j<NDIM;j++)
		    {
		      node_xmin[j]=xmin[j];
		      node_xmax[j]=xmax[j];
		    }
		  
		  j=1,k=0;

		  if(node->x[0][i]<=xmid[0])
		    {
		      j=0;
		      node_xmax[0]=xmid[0];
		    }
		  else
		    node_xmin[0]=xmid[0];
		  k+=j;

		  y[0]=node->x[0][i]; 
#if NDIM==2 || NDIM==3
		  j=1;
		  
		  if(node->x[1][i]<=xmid[1])
		    {
		      j=0;
		      node_xmax[1]=xmid[1];
		    }
		  else
		    node_xmin[1]=xmid[1];
		  
		  k+=2*j;

		  y[1]=node->x[1][i]; 		  
#endif
#if NDIM==3            
		  j=1;
		  
		  if(node->x[2][i]<=xmid[2])
		    {
		      j=0;
		      node_xmax[2]=xmid[2];
		    }
		  else
		    node_xmin[2]=xmid[2];
		  
		  k+=4*j;

		  y[2]=node->x[2][i];		  
#endif
		  if(k==0)
		    node->p0=add_to_node_nd(node->p0,y,node_xmin,node_xmax);
		  else if(k==1)
		    node->p1=add_to_node_nd(node->p1,y,node_xmin,node_xmax);
#if NDIM==2 || NDIM==3
		  else if(k==2)
		    node->p2=add_to_node_nd(node->p2,y,node_xmin,node_xmax);
		  else if(k==3)
		    node->p3=add_to_node_nd(node->p3,y,node_xmin,node_xmax);      
#endif
#if NDIM==3
		  else if(k==4)
		    node->p4=add_to_node_nd(node->p4,y,node_xmin,node_xmax);
		  else if(k==5)
		    node->p5=add_to_node_nd(node->p5,y,node_xmin,node_xmax);      
		  else if(k==6)
		    node->p6=add_to_node_nd(node->p6,y,node_xmin,node_xmax);
		  else if(k==7)
		    node->p7=add_to_node_nd(node->p7,y,node_xmin,node_xmax);      
#endif      
		}  // for(i=0;i<MAXNODE;i++)
	      
	      for(i=0;i<NDIM;i++)
		{
		  node_xmin[i]=xmin[i];
		  node_xmax[i]=xmax[i];
		}
	      
	      // Then add new particle
	      j=1,k=0;
	      
	      if(x[0]<=xmid[0])
		{
		  j=0;
		  node_xmax[0]=xmid[0];
		}
	      else
		node_xmin[0]=xmid[0];
	      k+=j;
	      
#if NDIM==2 || NDIM==3
	      j=1;
	      
	      if(x[1]<=xmid[1])
		{
		  j=0;
		  node_xmax[1]=xmid[1];
		}
	      else
		node_xmin[1]=xmid[1];
	      
	      k+=2*j;      
#endif
#if NDIM==3            
	      j=1;
	      
	      if(x[2]<=xmid[2])
		{
		  j=0;
		  node_xmax[2]=xmid[2];
		}
	      else
		node_xmin[2]=xmid[2];
	      
	      k+=4*j;      
#endif
	      if(k==0)
		node->p0=add_to_node_nd(node->p0,x,node_xmin,node_xmax);
	      else if(k==1)
		node->p1=add_to_node_nd(node->p1,x,node_xmin,node_xmax);
#if NDIM==2 || NDIM==3
	      else if(k==2)
		node->p2=add_to_node_nd(node->p2,x,node_xmin,node_xmax);
	      else if(k==3)
		node->p3=add_to_node_nd(node->p3,x,node_xmin,node_xmax);      
#endif
#if NDIM==3
	      else if(k==4)
		node->p4=add_to_node_nd(node->p4,x,node_xmin,node_xmax);
	      else if(k==5)
		node->p5=add_to_node_nd(node->p5,x,node_xmin,node_xmax);      
	      else if(k==6)
		node->p6=add_to_node_nd(node->p6,x,node_xmin,node_xmax);
	      else if(k==7)
		node->p7=add_to_node_nd(node->p7,x,node_xmin,node_xmax);      
#endif      
	    } // if(node->num_members>MAXNODE)
	  else
	    {
	      for(i=0;i<NDIM;i++)
		{
		  node->x[i][node->num_members-1]=x[i];
		  node->xmean[i]=0.5*(node->xmean[i]+x[i]);
		}
	    }
	} // if(node->split==0)
      else
	{
	  (node->num_members)++;
	  
	  for(i=0;i<NDIM;i++)
	    {
	      node_xmin[i]=xmin[i];
	      node_xmax[i]=xmax[i];
	    }
	  
	  // Then add new particle
	  j=1,k=0;
	  
	  if(x[0]<=xmid[0])
	    {
	      j=0;
	      node_xmax[0]=xmid[0];
	    }
	  else
	    node_xmin[0]=xmid[0];
	  k+=j;
	  
#if NDIM==2 || NDIM==3
	  j=1;
	  
	  if(x[1]<=xmid[1])
	    {
	      j=0;
	      node_xmax[1]=xmid[1];
	    }
	  else
	    node_xmin[1]=xmid[1];
	  
	  k+=2*j;      
#endif
#if NDIM==3            
	  j=1;
	  
	  if(x[2]<=xmid[2])
	    {
	      j=0;
	      node_xmax[2]=xmid[2];
	    }
	  else
	    node_xmin[2]=xmid[2];
	  
	  k+=4*j;      
#endif
	  if(k==0)
	    node->p0=add_to_node_nd(node->p0,x,node_xmin,node_xmax);
	  else if(k==1)
	    node->p1=add_to_node_nd(node->p1,x,node_xmin,node_xmax);
#if NDIM==2 || NDIM==3
	  else if(k==2)
	    node->p2=add_to_node_nd(node->p2,x,node_xmin,node_xmax);
	  else if(k==3)
	    node->p3=add_to_node_nd(node->p3,x,node_xmin,node_xmax);      
#endif
#if NDIM==3
	  else if(k==4)
	    node->p4=add_to_node_nd(node->p4,x,node_xmin,node_xmax);
	  else if(k==5)
	    node->p5=add_to_node_nd(node->p5,x,node_xmin,node_xmax);      
	  else if(k==6)
	    node->p6=add_to_node_nd(node->p6,x,node_xmin,node_xmax);
	  else if(k==7)
	    node->p7=add_to_node_nd(node->p7,x,node_xmin,node_xmax);      
#endif      
	}
    }

  return node;
}

struct btree_node *make_btree(struct btree_node *node, float x)
{
  if(node==NULL)
    {
      node=(struct btree_node *)malloc(sizeof(struct btree_node));
      node->x=x;
      node->left=node->right=NULL;
    }
  else
    {
      if(x<=node->x)
	node->left=make_btree(node->left,x);
      else
	node->right=make_btree(node->right,x);
    }
  return node;
}

void scan_btree(struct btree_node *node)
{
  if(node!=NULL)
    {
      scan_btree(node->left);
      printf("%g\n",node->x);
      scan_btree(node->right);
    }
}

void get_values(struct btree_node *node, float *x, int *idx, int ncount)
{
  if(node!=NULL)
    {
      if(*idx<ncount) 
	{
	  get_values(node->left,x,idx,ncount);
	  *x=node->x;
	  (*idx)++;
	  get_values(node->right,x,idx,ncount);
	}
    }
}

