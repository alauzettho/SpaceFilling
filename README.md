# SpaceFilling
Optimization Algorithm for Space Filling Designs

A quick start : source start.sh

Your main concern is the input.txt file:


ndim                  An int indicating the problem's dimension
npoint                An int indicating the number of points
optim_method 		      The optimization method. Choose between : NM	BFGS	RS. Best method is RS
estimation_method 	  The estimation method. Choose between   : MC	MST		NN	ROS



src/ folder contains all c++ files and main algorithms. Use:

  - "make" to compile
  - "make run" to start program
  - "make clean" tto clean all
  
To launch the program you can use ./spacefilling in src/


python/ folder contains simple matplotlib functions to plot our design
Functions are called at the end of the .cpp main
