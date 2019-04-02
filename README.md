# SpaceFilling
Optimization Algorithm for generating Space Filling Designs

A quick start : source start.sh



The src/ folder contains all c++ files and the main algorithm. Use:

  - "make" to compile
  - "make run" to start the program with generic values
  - "make clean" to clean all files



To launch the program with more options you can use the spacefilling binary in src/
To do so, you have to specify values:

  - ndim		to specify the problem's dimension
  - npoint		to specify the number of points in the design
  - optimMethod		to specify the method. Choose between : NM BFGS RS. Default is RS.
  - estimationMethod	to specify the criteria. Choose between : MC MST NN MM AE KL. Default is MST.


NM is the nelder Mead method
BFGS is a gradient based method (not suited here)
RS is for simulated annealing

MC is the Monte Carlo method
MST stands for Minimal Spanning Tree
NN is Nearest Neighbor
MM is for Maximin
AE is a method based on particules
KL is your typical Shannon Entropy



The python/ folder contains simple matplotlib functions to plot generated designs
