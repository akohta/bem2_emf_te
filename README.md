# bem2_emf_te  

This is the two-dimensional electromagnetic field analysis program for arbitrary objects irradiated by a TE plane wave (transverse electric wave).  
This is based on the boundary element method, the own developed numerical solution is used.
Intel Math Kernel Library and libpng are required. 
Gmsh is used for create a mesh data of object.  

![analysis model](model_te.png "analysis model")  

## Usage of example code  

1. type 'make' command to compile  
   The executable d2te_bv_solver, example1.out, example2.out, example3.out are created. 
   The executable d2te_bv_solver is the main solver of the boundary integral equations. 
   The example1.out is the executable of source code example1.c, it shows a simplest example using "bem2_emf_te". 
   The example2.out is the execubable of source code example2.c, it shows a example of electromagnetic field intensity analysis. 
   The example3.out is the executable of source code example3.c, it shows a example of outputting the instantaneous value of electromagnetic field as an image.

2. type './d2te_bv_solver' with arguments of incident field datafile name, medium datafile name, mesh datafile name and output dafafile name.  
   For example, './d2te_bv_solver ifd.txt medium_data.txt circle_1.msh ex.dat'. 
   The ifd.txt is the sample of incident field datafile, a TE plane wave is defined in it.
   The medium_data.txt is the sample of medium datafile, one medium is defined in it. The domain numbers are assigned to the medium from 1 in order. 
   The circle_1.msh is the example of mesh datafile, it is an object with a circular cross sectios. 
   It was created by using Gmsh geometry file circle_1.geo in the mesh_sample folder. 
   The d2te_bv_solver solves boundary integral equations with the specified datafiles, outputs the results to a binary file with the output datafile name. 
   It has optional arguments for rotation and translation of the object.
   For the rotation angle around the z-axis is theta and the translation vector is (tx, ty, tz), the arguments are './d2te_bv_solver ifd.txt medium_data.txt circle_1.msh ex.dat theta tx ty tz'.
   As a simple representation of the analysis model, the nodes used for the surface integral are output as a point cloud data. 
   In this example, the file ex.particles is output and the visualization result is particles.png (using Gnuplot script gscript_particles.plt).  
   
3. type './example1.out' with an argument of datafile name output by d2te_bv_solver.  
   For example, './example1.out ex.dat'. 
   This executable calculates electromagnetic field, radiation force and torque.  
  
4. type './example2.out' with an argument of datafile name output by d2te_bv_solver.  
   For example, './example2.out ex.dat'. 
   This executable calculates electromagnetic field intensity distributions, outputs them to text files. 
   The I_example2.png is the visualization result of intensity distributions, created by Gnuplot script gscript_example2.plt.  

5. type './example3.out' with an argument of datafile name output by d2te_bv_solver.  
   For example, './example3.out ex.dat'. 
   This executable calculates instantaneous value of the electromagnetic fields, outputs them to png image files. 
   The image files are output to the folder which has a name adding "images" to the datafile name specified in the argument (file-extension is excluded). 
   Each image file has a name that indicates the cross section, field component and number of time steps (ex. xy_Ez_014.png). 
   The color bar is output as color_bar.png in the same folder. 
   The range of color bar in each cross section is output to the xy_info.txt file. 
   The xy_Ez.gif, xy_Hx.gif and xy_Hy.gif are animated gifs that concatenate the png files created by using the shell script gif_animation.sh.  
   
Please see d2te_src/bem2_emf_te.h for detail of functions. 
The main parts of the code are parallelized by using OpenMP. 
The number of threads is controlled by the environment variable OMP_NUM_THREADS. 
The additional analysis examples are in the folder analysis_sample1 ~ analysis_sample4.

![intensity distributions](I_example2.png "intensity distributions (I_example2.png)")  
![model particles](particles.png "image of the object (particles.png)")![Ez gif](xy_Ez.gif "instantaneous value of the E_z (xy_Ez.gif)")  
![Hx_gif](xy_Hx.gif "instantaneous value of the H_x (xy_Hx.gif)")![Hy gif](xy_Hy.gif "instantaneous value of the H_y (xy_Hy.gif)")  


## Analysis sample 3 (in the folder analysis_sample3)  

![intensity distributions 3](analysis_sample3/I_example2.png "intensity distributions (analysis_sample3/I_example2.png)")  
![model particles 3](analysis_sample3/particles.png "image of the object (analysis_sample3/particles.png)")![Ez gif 3](analysis_sample3/xy_Ez.gif "instantaneous value of the E_z (analysis_sample3/xy_Ez.gif)")  
![Hx_gif 3](analysis_sample3/xy_Hx.gif "instantaneous value of the H_x (analysis_sample3/xy_Hx.gif)")![Hy gif 3](analysis_sample3/xy_Hy.gif "instantaneous value of the H_y (analysis_sample3/xy_Hy.gif)")  


## References  

1. Intel Math Kernel Library [MKL](https://software.intel.com/mkl)  
2. The official PNG reference library [libpng](http://www.libpng.org/pub/png/libpng.html)  
3. Three-dimensional mesh generator [Gmsh](https://gmsh.info/)  
4. The command-line driven graphing utility [gnuplot](http://www.gnuplot.info/)  
5. The utilities for manipulating images [ImageMagick](https://imagemagick.org/)  
