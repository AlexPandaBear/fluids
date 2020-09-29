# HEAT EQUATION SOLVER FOR UNSTEADY INCOMPRESSIBLE FLOWS

This project is a solver for the Navier-Stokes model in the case of 2D incompressible flows. Since the thermal and dynamic aspects of the flow are separable (see incompressible Navier-Stokes equations below), the code first computes the motion, and then solves the heat equation using the known motion.

<p align="center">
	<img src=eq_mass.png />
</p>

<p align="center">
	<img src=eq_momentum.png />
</p>

<p align="center">
	<img src=eq_energy.png />
</p>

The computation of the fluid motion is performed by explicitely integrating the velocity first, and then the pressure field is updated by solving a linear system AX = B with a LU decomposition LUX = B thanks to the Doolittle Algorithm.

Then, the computed motion is used in a $\theta$-schema for the temporal integration, and the Gauss-Seidel algorithm to solve the resulting linear system.

This project is coded as a C++ library interfaced with Python scripts, to benefit from the performance of the C++, the flexibility of Python and the graphical rendering of the Matplotlib library. The interface is realized with the PyBind11 library.

## Requirements
To compile the library :
- CMake  
- A C++ compiler  
- The PyBind11 library  

To use it with Python :
- Python 3+  
- The Numpy and Matplotlib libraries  

## How to install
Go to the code directory and execute the following commands to create a build directory to compile the library, a data directory to store the raw simulation outputs, and a results directory to save the postprocessed results.
```console
user@linux:~/path/to/code$ mkdir build  
user@linux:~/path/to/code$ mkdir data  
user@linux:~/path/to/code$ mkdir results  
```

Execute the following commands to compile the C++ library.
```console
user@linux:~/path/to/code$ cd build  
user@linux:~/path/to/code/build$ cmake ../src -DCMAKE_INSTALL_PREFIX=$PWD/install    
user@linux:~/path/to/code/build$ make install    
```

## How to use
Execute these commands to get to the python directory and add the library to the PYTHONPATH.
```console
user@linux:~/path/to/code$ cd python
user@linux:~/path/to/code/python$ export PYTHONPATH=$PWD/../build/install/lib/python
```

For each simulation, follow these steps :
- edit the parameters.py file
- execute the simulation.py file
- execute the plot.py file
