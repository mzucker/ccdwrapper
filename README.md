# ccdwrapper

Minimal example demonstrating how to fake separation queries with libccd.

## Building

Required packages:

  - Eigen
  - libccd
  - GLUT
  
 To compile
 
     cd /path/to/ccdwrapper
     mkdir build
     cd build
     cmake -DCMAKE_BUILD_TYPE=Release ..
     make
     
Running the demo

     ./ccddemo
     
Keys do stuff:

   - `ENTER` toggle animation
   - `+` or `-` change dilation buffer size (visualize using points, below)
   - `p` toggle point drawing
   - `w` toggle wireframe
   - `s` toggle bounding spheres
   - `a` toggle penetration algorithm (MPR works, GJK gives odd results)
   - `q` toggle query type (separation does not seem very useful right now)
   - `ESC` quits
