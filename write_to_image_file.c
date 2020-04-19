#include <stdio.h>
#include <stdlib.h>

#include "header.h"

void write_to_ppm(char *image_file, int width, int height, int MaxColorComponentValue, float *data)
{
  FILE *fp = fopen(image_file,"wb");
  fprintf(fp,"P5\t %d \t %d \t %d\n",width,height,MaxColorComponentValue);
  /*
  for (int j = 0; j < height; ++j)
    for (int i = 0; i < width; ++i)
      {
	static unsigned char color[3];
	color[0] = 256-*((data+j*width)+i); 
	color[1] = 0;  
	color[2] = *((data+j*width)+i);  
	printf("%u\n",color[2]);
	
	  (void) fwrite(color, 1, 3, fp);
      }
  */
  //fwrite((unsigned char)data,width*height*sizeof(unsigned char),1,fp);
  fclose(fp);
}

void write_to_png(char *image_file, int width, int height, float *data)
{
  int i,j;
  bitmap_t density_map;

  density_map.width=width;
  density_map.height=height;
  density_map.pixels=calloc(width*height,sizeof(pixel_t));

  for(j=0;j<height;j++)
    for(i=0;i<width;i++)
      {
	pixel_t *pixel=pixel_at(&density_map,i,j);
	pixel->blue=(int)(*((data+j*width)+i)*255);
	pixel->green=(int)(*((data+j*width)+i)*255);
	pixel->red=(int)(*((data+j*width)+i)*255);
      }

  save_png_to_file(&density_map,image_file);
}

/* Given "bitmap", this returns the pixel of bitmap at the point
   ("x", "y"). */

pixel_t * pixel_at (bitmap_t * bitmap, int x, int y)
{
  return bitmap->pixels + bitmap->width * y + x;
}

/* Write "bitmap" to a PNG file specified by "path"; returns 0 on
   success, non-zero on error. */

int save_png_to_file (bitmap_t *bitmap, const char *path)
{
  FILE * fp;
  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;
  size_t x, y;
  png_byte ** row_pointers = NULL;
  /* "status" contains the return value of this function. At first
     it is set to a value which means 'failure'. When the routine
     has finished its work, it is set to a value which means
     'success'. */
  int status = -1;
  /* The following number is set by trial and error only. I cannot
     see where it it is documented in the libpng manual.
  */
  int pixel_size = 3;
  int depth = 8;
  
  fp = fopen (path, "wb");
  if (! fp) {
    goto fopen_failed;
  }
  
  png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    goto png_create_write_struct_failed;
  }
  
  info_ptr = png_create_info_struct (png_ptr);
  if (info_ptr == NULL) {
    goto png_create_info_struct_failed;
  }
  
  /* Set up error handling. */
  
  if (setjmp (png_jmpbuf (png_ptr))) {
    goto png_failure;
  }
  
  /* Set image attributes. */
  
  png_set_IHDR (png_ptr,
		info_ptr,
		bitmap->width,
		bitmap->height,
		depth,
		PNG_COLOR_TYPE_RGB,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
  
  /* Initialize rows of PNG. */
  
  row_pointers = png_malloc (png_ptr, bitmap->height * sizeof (png_byte *));
  for (y = 0; y < bitmap->height; y++) {
    png_byte *row =
      png_malloc (png_ptr, sizeof (uint8_t) * bitmap->width * pixel_size);
    row_pointers[y] = row;
    for (x = 0; x < bitmap->width; x++) {
      pixel_t * pixel = pixel_at (bitmap, x, y);
      *row++ = pixel->red;
      *row++ = pixel->green;
      *row++ = pixel->blue;
    }
  }
  
  /* Write the image data to "fp". */
  
  png_init_io (png_ptr, fp);
  png_set_rows (png_ptr, info_ptr, row_pointers);
  png_write_png (png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
  
  /* The routine has successfully written the file, so we set
     "status" to a value which indicates success. */
  
  status = 0;
  
  for (y = 0; y < bitmap->height; y++) {
    png_free (png_ptr, row_pointers[y]);
  }
  png_free (png_ptr, row_pointers);
  
 png_failure:
 png_create_info_struct_failed:
  png_destroy_write_struct (&png_ptr, &info_ptr);
 png_create_write_struct_failed:
  fclose (fp);
 fopen_failed:
  return status;
}

