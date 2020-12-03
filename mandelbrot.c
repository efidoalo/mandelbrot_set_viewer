/*=============================;
 *
 * File: mandelbrot.c
 * Content: The mandelbrot set
 * GUI allows selection of an area
 * and then a zoom in on that area
 * Date: 29/11/2020
 *
 ******************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gtk/gtk.h>
#include <gmp.h>
#include <mpfr.h>


char *filename = "mandelbrot.ppm";
mpfr_t init_x;
mpfr_t init_y;
mpfr_t width;
mpfr_t height;
mpfr_t dx;
mpfr_t dy;
mpfr_t temp;

GtkWidget *grid = NULL;
GtkWidget *image = NULL;
GtkWidget *window = NULL;

void compute_mandelbrot_set(mpfr_t x,
                            mpfr_t y,
                            mpfr_t width,
                            mpfr_t height,
                            char *filename);

static void
reset_mandelbrot(GtkWidget *widget,
             gpointer   image)
{
  
  /*init_x = -2.0;
  init_y = 1.0;
  width = 3.0;
  height = 2.0;
  dx = 3.0/599.0;
  dy = 2.0/399.0; */
  mpfr_set_d(init_x, -2.0, MPFR_RNDU);
  mpfr_set_d(init_y, 1.0, MPFR_RNDD);
  mpfr_set_d(width, 3.0, MPFR_RNDD);
  mpfr_set_d(height, 2.0, MPFR_RNDD);
  mpfr_set_d(dx,  3.0/599.0, MPFR_RNDD);
  mpfr_set_d(dy, 2.0/399.0, MPFR_RNDD);
  compute_mandelbrot_set(init_x, init_y, width, height, filename);
  gtk_grid_attach( GTK_GRID(grid), image, 0 , 0, 600, 400);
  gtk_image_set_from_file(image, filename);
}


void drag_begin_callback(GtkGestureDrag *gesture,
		         double start_x,
			 double start_y,
			 gpointer data)
{
  // init_x = init_x + (start_x*dx);
  mpfr_mul_d(temp, dx, start_x, MPFR_RNDD);
  mpfr_add(init_x, init_x, temp, MPFR_RNDD);
  
  //init_y = init_y - (start_y*dy);
  mpfr_mul_d(temp, dy, start_y, MPFR_RNDD);
  mpfr_sub(init_y, init_y, temp, MPFR_RNDD); 
}

void drag_end_callback(GtkGestureDrag *gesture,
		       double offset_x,
		       double offset_y,
		       gpointer image)
{
/*  width = offset_x*dx;
  height = offset_y*dy;
  dx = width/599.0;
  dy = height/399.0;
  */
  mpfr_mul_d(width, dx, offset_x, MPFR_RNDD);
  mpfr_mul_d(height, dy, offset_y, MPFR_RNDD);  
  mpfr_div_d(dx, width, 599.0, MPFR_RNDD);
  mpfr_div_d(dy, height, 399.0, MPFR_RNDD);


  //mpfr_printf("init_x:%.10Rf init_y:%.10Rf width:%.10Rf  height:%.10Rf\n", init_x, init_y, width, height);
  compute_mandelbrot_set(init_x, init_y, width, height, filename);
  gtk_image_set_from_file(image, filename);

}

static void
activate (GtkApplication *app,
          gpointer        user_data)
{
  GtkWidget *button;

  /* create a new window, and set its title */
  window = gtk_application_window_new (app);
  gtk_window_set_title (GTK_WINDOW (window), "Mandelbrot Set");

  gtk_container_set_border_width (GTK_CONTAINER (window), 10);

  /* Here we construct the container that is going pack our buttons */
  grid = gtk_grid_new ();
  
  /* Pack the container in the window */
  gtk_container_add (GTK_CONTAINER (window), grid);

  image = gtk_image_new_from_file("mandelbrot.ppm");
  gtk_grid_attach( GTK_GRID(grid), image, 0 , 0, 600, 400);

  button = gtk_button_new_with_label ("Reset");
  g_signal_connect (button, "clicked", G_CALLBACK (reset_mandelbrot), image);

  /* Place the first button in the grid cell (0, 0), and make it fill
   * just 1 cell horizontally and vertically (ie no spanning)
   */
  gtk_grid_attach (GTK_GRID (grid), button, 600, 0, 1, 20);

  GtkGesture *drag_gesture = gtk_gesture_drag_new(window);
  g_signal_connect(drag_gesture, 
		   "drag-begin",
		   G_CALLBACK(drag_begin_callback),
		   NULL);

  g_signal_connect(drag_gesture,
		   "drag-end",
		   G_CALLBACK(drag_end_callback),
		   image);

  gtk_window_resize( GTK_WINDOW(window), 630, 430);
  /* Now that we are done packing our widgets, we show them all
   * in one go, by calling gtk_widget_show_all() on the window.
   * This call recursively calls gtk_widget_show() on all widgets
   * that are contained in the window, directly or indirectly.
   */
  gtk_widget_show_all (window);

}

// Computes the mandelbrot set within the rectangle
// having top left corner (x,y), width and height.
// The mandelbrot set is written to the file given
// by filename (which must have extension .ppm)
// and is a PPM image. If there already is a file
// named filename in the current directory then the
// file is overwritten with the ppm mandelbrot set
// image.
// The image has 600 pixels in width and
// 400 pixels in height.
// (x + y*i)*(x + y*i) = (x^2 - y^2) + 2*x*y*i
void compute_mandelbrot_set(mpfr_t x, 
		            mpfr_t y, 
			    mpfr_t width, 
			    mpfr_t height,
			    char *filename)
{
  int NoOfIterations = 200;

  mpfr_div_d(dx, width, 599.0, MPFR_RNDD);
  mpfr_div_d(dy, height, 399.0, MPFR_RNDD);
  
  // ensure dx == dy so that images are scaled correctly
  if ( mpfr_cmp(dx, dy) >0 ) {
    mpfr_set(dy, dx, MPFR_RNDD);
  } 
  else {
    mpfr_set(dx, dy, MPFR_RNDD);
  }

  mpfr_t max, curr_x, curr_y, constant_x, constant_y,
	 temp_x, temp_x_squared, temp_y_squared,
	 abs_val, uniform_pixel_val;
  
  mpfr_init2(max, 200);
  mpfr_set_d(max, -1.0, MPFR_RNDU);
  mpfr_init2(curr_x, 200);
  mpfr_set_d(curr_x, 0.0, MPFR_RNDU);
  mpfr_init2(curr_y, 200);
  mpfr_set_d(curr_y, 0.0, MPFR_RNDU);
  mpfr_init2(temp_x_squared, 200);
  mpfr_set_d(temp_x_squared, 0.0, MPFR_RNDU);
  mpfr_init2(temp_y_squared, 200);
  mpfr_set_d(temp_y_squared, 0.0, MPFR_RNDU);
 
  mpfr_t array[400][600];
  for (int i=0; i<400; ++i) {
    for (int j=0; j<600; ++j) {
      mpfr_init2(array[i][j], 200);
      mpfr_set_d(array[i][j], 0.0, MPFR_RNDU);
    }
  }
  
  mpfr_init2(constant_x, 200);
  mpfr_set_d(constant_x, 0.0, MPFR_RNDU);
  mpfr_init2(constant_y, 200);
  mpfr_set_d(constant_y, 0.0, MPFR_RNDU);
  mpfr_init2(temp_x, 200);
  mpfr_set_d(temp_x, 0.0, MPFR_RNDU);
  mpfr_init2(abs_val, 200);
  mpfr_set_d(abs_val, 0.0, MPFR_RNDU);
  
  mpfr_init2(uniform_pixel_val, 200);
  mpfr_set_d(uniform_pixel_val, 0, MPFR_RNDU);
  // Generate absolute values after iterations and
  // store in array
  double cut_off = 5.0;
  for (int i=0; i<400; ++i) {
    for (int j=0; j<600; ++j) {
      
      // curr_x = x + ((double)j)*dx;	    
      mpfr_mul_d(curr_x, dx, ((double)j), MPFR_RNDU);
      mpfr_add(curr_x, curr_x, x, MPFR_RNDD);
      
      //curr_y = y - ((double)i)*dy;
      mpfr_mul_d(curr_y, dy, ((double)i), MPFR_RNDU);
      mpfr_sub(curr_y, y, curr_y, MPFR_RNDD);
      
      mpfr_set(constant_x, curr_x, MPFR_RNDD); // constant_x = curr_x;
      mpfr_set(constant_y, curr_y, MPFR_RNDD); // constant_y = curr_y;
      
      mpfr_set_d(temp_x, 0.0, MPFR_RNDU); // temp_x = 0.0;
      
      
      for (int k=1; k<NoOfIterations; ++k) {
	mpfr_set(temp_x, curr_x, MPFR_RNDD); // temp_x = curr_x;
	
	// curr_x = (curr_x*curr_x) - (curr_y*curr_y) + constant_x;
	mpfr_mul(temp_x_squared, curr_x, curr_x, MPFR_RNDU);
	mpfr_mul(temp_y_squared, curr_y, curr_y, MPFR_RNDU);
	mpfr_sub(curr_x, temp_x_squared, temp_y_squared, MPFR_RNDD);
	mpfr_add(curr_x, curr_x, constant_x, MPFR_RNDD);
 
        // curr_y = (2*temp_x*curr_y) + constant_y;
        mpfr_mul_d(curr_y, curr_y, 2.0, MPFR_RNDU);
        mpfr_mul(curr_y, temp_x, curr_y, MPFR_RNDU);
        mpfr_add(curr_y, curr_y, constant_y, MPFR_RNDD);	

        /* if ( sqrt( (curr_x*curr_x) + (curr_y*curr_y) ) > cut_off ) {
             break;
           } */
	mpfr_mul(temp_x_squared, curr_x, curr_x, MPFR_RNDU);
	mpfr_mul(temp_y_squared, curr_y, curr_y, MPFR_RNDU);
	mpfr_add(temp_x_squared, temp_x_squared, temp_y_squared,
		 MPFR_RNDD);

	mpfr_sqrt(temp_x_squared, temp_x_squared, MPFR_RNDU);
	
        if (mpfr_cmp_d(temp_x_squared, cut_off) > 0) {
          break;
	}	
      }
      
      //double abs_val = sqrt( (curr_x*curr_x) + (curr_y*curr_y) );
      mpfr_mul(temp_x_squared, curr_x, curr_x, MPFR_RNDU);
      mpfr_mul(temp_y_squared, curr_y, curr_y, MPFR_RNDU);
      mpfr_add(abs_val, temp_x_squared, temp_y_squared, MPFR_RNDD);
      mpfr_sqrt(abs_val, abs_val, MPFR_RNDU);
      
      
      // if (abs_val > cut_off) {
      //   array[i][j] = cut_off;
      //   max = cut_off;
      // }
      // else {
      //   array[i][j] = abs_val
      //   if (abs_val > max)
      //     max = abs_val;
      // }

      if (mpfr_cmp_d(abs_val, cut_off) > 0) {
        mpfr_set_d( array[i][j], cut_off, MPFR_RNDD);
	mpfr_set_d(max, cut_off, MPFR_RNDD);
      }
      else {
        mpfr_set(array[i][j], abs_val, MPFR_RNDD);
	if (mpfr_cmp(abs_val, max) > 0)
          mpfr_set(max , abs_val, MPFR_RNDD);
      }
    }
  }

  unsigned long int get_red_byte_bitmask = 0x0000000000ff0000,
                  get_green_byte_bitmask = 0x000000000000ff00,
                   get_blue_byte_bitmask = 0x00000000000000ff;

  FILE *fp = fopen(filename, "w");  
  char *header = "P3\n600 400\n255";
  fwrite(header, 1, strlen(header), fp);
  for (int i=0; i<400; ++i) { 
    for (int j=0; j<600; ++j) {
      // int uniform_pixel_val = (int)((array[i][j])*( (0xffffff)/max));
      mpfr_d_div(uniform_pixel_val, 16777215.0, max, MPFR_RNDD);  
      mpfr_mul(uniform_pixel_val, uniform_pixel_val, array[i][j], MPFR_RNDD);
      mpfr_rint(uniform_pixel_val, uniform_pixel_val, MPFR_RNDD);
      unsigned long int rgb = mpfr_get_ui(uniform_pixel_val, MPFR_RNDD);
      unsigned char red_byte = (rgb & get_red_byte_bitmask) >> 16;
      unsigned char green_byte = (rgb & get_green_byte_bitmask) >> 8;
      unsigned char blue_byte = rgb & get_blue_byte_bitmask;

      //printf("colour:%u %u %u ", red_byte, green_byte, blue_byte);
      //mpfr_printf("%.0Rf\n", uniform_pixel_val);
      if (j==0)
        fprintf(fp, "\n%u %u %u", red_byte, 
			                        blue_byte,
			                        green_byte);
      else
	fprintf(fp, " %u %u %u", red_byte, 
	  		                  blue_byte,
			                  green_byte);
    }
  }
  
  fclose(fp);
  mpfr_clear(max);
  mpfr_clear(curr_x);
  mpfr_clear(curr_y);
  mpfr_clear(constant_x);
  mpfr_clear(constant_y);
  mpfr_clear(temp_x);
  mpfr_clear(temp_x_squared);
  mpfr_clear(temp_y_squared);
  mpfr_clear(abs_val);
  mpfr_clear(uniform_pixel_val);

  for (int i=0; i<400; ++i) {
    for (int j=0; j<600; ++j) {
      mpfr_clear( array[i][j] );
    }
  }
}

int
main (int    argc,
      char **argv)
{
  GtkApplication *app;
  int status;
  
  mpfr_init2(init_x, 200);
  mpfr_set_d(init_x, -2.0, MPFR_RNDU);
  mpfr_init2(init_y, 200);
  mpfr_set_d(init_y, 1.0, MPFR_RNDD);
  mpfr_init2(width, 200);
  mpfr_set_d(width, 3.0, MPFR_RNDD);  
  mpfr_init2(height, 200);
  mpfr_set_d(height, 2.0, MPFR_RNDD);
  mpfr_init2( dx, 200);
  mpfr_set_d(dx,  3.0/599.0, MPFR_RNDD);
  mpfr_init2(dy, 200);
  mpfr_set_d(dy, 2.0/399.0, MPFR_RNDD);
  mpfr_init2(temp, 200);
  mpfr_set_d(temp, 0.0, MPFR_RNDU);

  compute_mandelbrot_set(init_x, init_y, width, height, filename); 
  app = gtk_application_new ("mandelbrot.set", G_APPLICATION_FLAGS_NONE);
  g_signal_connect (app, "activate", G_CALLBACK (activate), NULL);
  status = g_application_run (G_APPLICATION (app), argc, argv);
  g_object_unref (app);
  
  mpfr_clear(init_x);
  mpfr_clear(init_y);
  mpfr_clear(width);
  mpfr_clear(height);
  mpfr_clear(dx);
  mpfr_clear(dy);
  mpfr_clear(temp);
  mpfr_free_cache();

  return status;
}
