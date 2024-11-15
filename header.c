#include "header.h"
#include <png.h>

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

void check_input_filenames(char *, char *, int, int *);
void read_hdf5_header(char *, struct sim_info *, long long *);
void read_gadget_binary_header(char *, struct sim_info *, long long *);
void read_particles_from_hdf5(char *, float *, float *, float *, int *, int, long long *);
void read_particles_from_gadget_binary(char *, float *, float *, float *, int *, int, long long *);
void select_particles(float *, float *, float *, int *, double, long long, float, float, float, float, int, long long *);
void smooth_to_mesh(long long, float *, float *, float *, float *, float, float, float, float, float, int, int, float *);
void split_across_tasks_as_slabs(float *, float *, float *, long long *, float, float, float, float);

void write_to_ppm(char *, int, int, int, float *);
void write_to_png(char *, int, int,float *);
pixel_t * pixel_at (bitmap_t *, int, int);
int save_png_to_file (bitmap_t *, const char *);

int ThisTask, NTask;
int SnapFormat;
