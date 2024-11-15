#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <hdf5.h>

#include "header.h"

#define BUFSIZE 500

#define rhocrit 27.755

#include <png.h>

#include <mpi.h>

void find_neighbours(int num_part,float *smoothing_length, int num_ngb,
		     float *posx, float *posy, float *posz,
		     float xmin, float ymin, float zmin);


float finterpolant(float xmin, float xmax);

#define IMAGE_DIMENSION 1024
#define IMAGE_DIMENSIONX 256
#define IMAGE_DIMENSIONY 256
#define MAX_COLOUR_COMPONENT_VALUE 255
#define FMIN 2.5e-10

int main(int argc, char *argv[])
{
    int i,j,k,l,m,n,buf[BUFSIZE];
    char file_root[128],filename[128],image_file[128],image_file_root[128],debug_file[128];
    FILE *infile,*outfile;
    char buffer[BUFSIZE];
    double SliceSize,ImageSliceSize,f_Image=0.1;
    double mp;
    
    char junk[256-24];
    long long NumPart=0, NumPartInFile=0;
    long long fileoffset;
    
    float *x,*y,*z;
    int *ptype;
    
    float *smoothing_length;
    float xmin,ymin,zmin;
    float xmax,ymax,zmax;
    
    
    double dx,dy,dz;
    
    int isHDF5=0;
    int isDistributed=0;
    
    struct sim_info header;
    
    long long NumPartRead=0;
    
    int ptype_keep=1;  // Default - assume we keep only dark matter
    
    float BoundingBox;

    double xcen=0.5, ycen=0.5, zcen=0.5;
    
    double xc=0.0,yc=0.0,zc=0.0,lbox=0.0;

    int itmax;
    
    clock_t tstart, tfinish;
    
#ifdef ENABLE_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD,&NTask);
#else  
    ThisTask=0;
    NTask=1;
#endif
    
        // Now check that user input is correct
    if(ThisTask==0)
        if(argc<2)
        {
            fprintf(stdout,"Usage: %s -input <filename> (Input file)\n",argv[0]);
            fprintf(stdout,"          -isHDF5 (HDF5 format; optional)\n");
            fprintf(stdout,"          -output <filename> (Output file)\n");
            fprintf(stdout,"          -xc (x centre, box units; optional)\n");
            fprintf(stdout,"          -yc (y centre, box units; optional)\n");
            fprintf(stdout,"          -zc (z centre, box units; optional)\n");
            fprintf(stdout,"          -lbox (assumed box size, box units; optional\n");
            fprintf(stdout,"          -itmax (maximum number of iterations; optional\n");
            fprintf(stdout,"          -stars (plot stars only; optional\n");
            fprintf(stdout,"          -gas (plot gas only; optional\n");
            fprintf(stdout,"          -dark_matter (plot dark matter only; optional\n");
            fflush(stdout);
            exit(0);
        }
    
    i=1;
    
    while(i<argc)
    {
        if(strcmp(argv[i],"-input")==0)
            sprintf(file_root,"%s",argv[++i]);
        
        if(strcmp(argv[i],"-output")==0)
            sprintf(image_file_root,"%s",argv[++i]);
        
        if(strcmp(argv[i],"-xc")==0)
            xcen=atof(argv[++i]);
        
        if(strcmp(argv[i],"-yc")==0)
            ycen=atof(argv[++i]);
        
        if(strcmp(argv[i],"-zc")==0)
            zcen=atof(argv[++i]);
        
        if(strcmp(argv[i],"-lbox")==0)
            lbox=atof(argv[++i]);
        
        if(strcmp(argv[i],"-itmax")==0)
            itmax=atoi(argv[++i]);
        
        if(strcmp(argv[i],"-isHDF5")==0) isHDF5=1;
        if(strcmp(argv[i],"-stars")==0) ptype_keep=4;
        if(strcmp(argv[i],"-gas")==0) ptype_keep=0;
        if(strcmp(argv[i],"-dark_matter")==0) ptype_keep=1;
        
        i++;
    }
    
        // Get filename and figure out if it's distributed across multiple files
    check_input_filenames(filename,file_root,isHDF5,&isDistributed);
    
        // Echo user defined inputs to screen
    if(ThisTask==0)
    {
        if(xc>0 || lbox>0)
        {
            fprintf(stdout,"Centre: (%g|%g|%g) \n",xcen,ycen,zcen);
            fprintf(stdout,"Box Length: %g\n",lbox);
        }
    }
    
    if(ThisTask==0)
    {
        fprintf(stdout,"Reading header...\n");
        fflush(stdout);
    }
    
    if(isHDF5==1)
        read_hdf5_header(filename, &header, &NumPart);
    else
        read_gadget_binary_header(filename, &header, &NumPart);
    
    if(header.BoxSize==0) 
        header.BoxSize=1.e6;  // This can happen in non-cosmological runs
    else
    {
        xcen*=header.BoxSize;
        ycen*=header.BoxSize;
        zcen*=header.BoxSize;
        lbox*=header.BoxSize;
    }
    
    xc=xcen; yc=ycen; zc=zcen;
    
    printf("%g %g %g %g\n",xc,yc,zc,lbox);
    
    if(ThisTask==0)
    {
        fprintf(stdout,"Number of files: %d\n",header.NumFiles);
        fprintf(stdout,"Number of particles: %lld\n",NumPart);
        fprintf(stdout,"Time/Expansion Factor: %g\n",header.time);
        fflush(stdout);
    }
    
        // Allocate memory....
    
    x = (float*)malloc(sizeof(float)*NumPart);
    y = (float*)malloc(sizeof(float)*NumPart);
    z = (float*)malloc(sizeof(float)*NumPart);
    ptype = (int*)malloc(sizeof(int)*NumPart);
    
    if(ThisTask==0)
    {
        fprintf(stdout,"Allocated memory for positions...\n");
        fflush(stdout);
    }
    
    if(ThisTask==0)
    {
        fprintf(stdout,"Reading particles...\n");
        fflush(stdout);
    }
    
    if(isHDF5==1)
        read_particles_from_hdf5(file_root,x,y,z,ptype,header.NumFiles,&NumPartRead);
    else
        read_particles_from_gadget_binary(file_root,x,y,z,ptype,header.NumFiles,&NumPartRead);
    
    select_particles(x,y,z,ptype,header.BoxSize,NumPartRead,xc,yc,zc,lbox,ptype_keep,&NumPart);
#ifdef ENABLE_MPI
    split_across_tasks_as_slabs(x,y,z,&NumPart,xc,yc,zc,lbox);
    
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    if(NumPart==0)
    {
        if(ThisTask==0)
        {
            fprintf(stdout,"Error: no particles selected...\n");
            fflush(stdout);
        }
        exit(0);
    }
    else
    {
        if(ThisTask==0)
        {
            fprintf(stdout,"Plotting %d particles...\n",NumPart);
            fflush(stdout);
        }
    }
    
    smoothing_length = (float*)malloc(2*sizeof(float)*(NumPart));
    
#ifdef KERNEL_SMOOTHING
    float dxmin[3]={0,0,0};
    int num_ngb=40;
    
    tstart=clock();
    find_neighbours(NumPart,smoothing_length,num_ngb,x,y,z,dxmin[0],dxmin[1],dxmin[2]);
    tfinish=clock()-tstart;
    
    fprintf(stdout,"find_neighbours - time: %lu s...\n",tfinish/CLOCKS_PER_SEC);
    fflush(stdout);
#endif
    
    float data[IMAGE_DIMENSIONX][IMAGE_DIMENSIONY];
    float global_data[IMAGE_DIMENSIONX][IMAGE_DIMENSIONY];
    
    BoundingBox=lbox;
    float theta=0.;
    int iter=0;
    
#define LOG10DENSMAX 8
#define LOG10DENSMIN 3
    
    
    while(BoundingBox>0.1*lbox)
    {
        for(j=0; j<IMAGE_DIMENSIONY; j++)
            for(i=0; i<IMAGE_DIMENSIONX; i++)
                data[i][j]=global_data[i][j]=0.0;
        
        tstart=clock();
        smooth_to_mesh(NumPart,smoothing_length,x,y,z,xc,yc,zc,BoundingBox,theta,IMAGE_DIMENSIONX,IMAGE_DIMENSIONY,*data);
        tfinish=clock()-tstart;
        fprintf(stdout,"smooth_to_mesh - time: %lu s...\n",tfinish/CLOCKS_PER_SEC);
        fflush(stdout);
        
#ifdef ENABLE_MPI
        MPI_Reduce(data, global_data, IMAGE_DIMENSIONX*IMAGE_DIMENSIONY, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        for(j=0; j<IMAGE_DIMENSIONY; j++)
            for(i=0; i<IMAGE_DIMENSIONX; i++)
            {
                global_data[i][j]=data[i][j];
                printf("%lf\n",global_data[i][j]);
            }
#endif
        
        for(j=0; j<IMAGE_DIMENSIONY; j++)
            for(i=0; i<IMAGE_DIMENSIONX; i++)
            {
                if(log10f(global_data[i][j])<LOG10DENSMIN) global_data[i][j]=pow(10,LOG10DENSMIN);
                global_data[i][j]=(log10f(global_data[i][j])-LOG10DENSMIN)/(LOG10DENSMAX-LOG10DENSMIN);
            }
        
        if(ThisTask==0)
        {
            sprintf(image_file,"%s.%04d.%s",image_file_root,iter,"png");
            
//            write_to_ppm(image_file,IMAGE_DIMENSION,IMAGE_DIMENSION,MAX_COLOUR_COMPONENT_VALUE,*data);
            write_to_png(image_file,IMAGE_DIMENSIONX,IMAGE_DIMENSIONY,*global_data);
        }
        
            //      BoundingBox/=1.01;
        if(BoundingBox>0.55*lbox)
            theta+=2*M_PI/180.0;
        if(BoundingBox<0.5*lbox)
            theta+=2*M_PI/180.0;
        iter++;
        if(iter>itmax-1) break;
    }
    
    
    fprintf(stdout,"Finished...\n");
    fflush(stdout);
    
#ifdef ENABLE_MPI
    MPI_Finalize();
#endif
}
