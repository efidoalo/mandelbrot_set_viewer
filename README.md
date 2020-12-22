# mandelbrot_set_viewer
A program that enables selection and zoom-in on areas of the mandelbrot set.
The program uses GTK+ (https://www.gtk.org/) for the GUI and the multiple 
precision floating numbers with accurate roumdimg (mpfr) library for additional 
precision when calculating points in (and out of) the set. The MPFR library can 
be downloaded here https://www.mpfr.org/.

Compile with: 
gcc -c -I/home/andy/Documents/projects/mandelbrot/gmp-6.2.1 `'pkg-config --cflags gtk+-3.0'` mandelbrot.c .
Link with:
gcc mandelbrot.o -o mandelbrot `'pkg-config --libs gtk+-3.0'` -lm -lmpfr . Where the -I compiler option specifies where to look for the gmp heder (<gmp.h>). The mpfr lirbary is based on the gnu multiple precision (gmp) library.
   
