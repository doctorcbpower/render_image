#include <stdint.h>
#include <png.h>

#define MAXNODE 8
#define NDIM 3

struct link_list {
  float x[NDIM];
  struct link_list *ptr;
};

struct interaction_list {
  int index;
  struct tree_node_nd *node;
  struct interaction_list *left, *right;
};

struct tree_node {
  float x[MAXNODE];
  float xmin,xmax;
  int num_members;
  int i;
  int split;
  struct tree_node *left, *right;
};

struct tree_node_nd {
  float x[NDIM][MAXNODE];
  float xmin[NDIM], xmax[NDIM];
  float xmean[NDIM];
  int num_members;
  int index;
  int split;
  struct tree_node_nd *p0,*p1;
#if NDIM==2 || NDIM==3
  struct tree_node_nd *p2,*p3;
#endif
#if NDIM==3
  struct tree_node_nd *p4,*p5,*p6,*p7;
#endif  
};

struct btree_node {
  float x;
  struct btree_node *left, *right;
};

struct point {
  float pos[NDIM];
  float mass;
  float density;
  float dist_ngb;
};

#define NTAB 100
float rh_table[NTAB];
float kernel_table[NTAB];
float integral_table[NTAB];

struct tree_node *add_to_node(struct tree_node *, float, float, float);
void get_node(struct tree_node *, float, int *, struct link_list **);
void get_multiple_nodes(struct tree_node *, float, float, int *, struct link_list **);
void scan_nodes(struct tree_node *, int *);
float get_rand(int);

struct btree_node *make_btree(struct btree_node *, float);
void scan_btree(struct btree_node *);

struct tree_node_nd *add_to_node_nd(struct tree_node_nd *,float *,float *,float *);
void get_node_nd(struct tree_node_nd *, float *, int *, struct link_list **);
void get_multiple_nodes_nd(struct tree_node_nd *, float *, float, int *, int *, struct link_list **);
void sweep_over_nodes_nd(struct tree_node_nd *, float *, float, int *, int *);
void add_node_index_nd(struct tree_node_nd *, int *);
void build_interaction_list(struct tree_node_nd *, struct interaction_list **, int *);
void get_values(struct btree_node *, float *, int *, int);

void get_points(struct point *, int *, int);
void make_tree(struct point *, int, struct tree_node_nd **);
void get_distance_to_nth_nearest_neighbour(struct point *, int, int, struct tree_node_nd *);
void get_kernel_density_estimate(struct point *, int, struct tree_node_nd *);
void get_potential_estimate(struct point *, int, struct tree_node_nd *);

int cmpfunc(const void*, const void*);

void tabulate_kernel(void);
float cubic_spline_kernel(float);
float get_kernel_value(float);

float cumulative_cubic_spline_interpolant(float);
void tabulate_integral(void);
float get_integral_value(float);

struct sim_info {
  int npart[6],nall[6];
  double massarr[6];
  double time;
  int NumFiles;
  double BoxSize;
  double Omega0,OmegaLambda,HubbleParam;  
  int SnapFormat;
};

void read_hdf5_header(char *, struct sim_info *, long long *);
void read_gadget_binary_header(char *, struct sim_info *, long long *);
void read_particles_from_hdf5(char *, float *, float *, float *, int *, int, long long *);
void read_particles_from_gadget_binary(char *, float *, float *, float *, int *, int, long long *);
void select_particles(float *, float *, float *, int *, double, long long, float, float, float, float, int, long long *);
void split_across_tasks_as_slabs(float *, float *, float *, int **, long long, float);
void smooth_to_mesh(long long, float *, int, float *, float *, float *, float, float, float, float, float, int, int, float *);

void write_to_ppm(char *, int, int, int, float *);

typedef struct
{
  uint8_t red;
  uint8_t green;
  uint8_t blue;
}
  pixel_t;

/* A picture. */

typedef struct
{
  pixel_t *pixels;
  size_t width;
  size_t height;
}
  bitmap_t;

void write_to_png(char *, int, int, float *);
pixel_t * pixel_at (bitmap_t *, int, int);
int save_png_to_file (bitmap_t *, const char *);

extern int ThisTask, NTask;
extern int SnapFormat;
extern int NThisTask;
