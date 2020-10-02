# 2D INC. NAVIER-STOKES SOLVER

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

The computation of the fluid motion is performed by iteratively solving a linear system AX = B at each time step with the Gauss-Seidel algorithm. The Navier-Stokes model can be reduced to this linear system by discretizing all the field operators and the temporal derivatives : the field operators are approximated by centered 2nd order finite differences, and the temporal integration is performed by a θ-scheme.

Then, the computed motion is used in another θ-scheme associated with another Gauss-Seidel algorithm to solve the resulting linear system for the temporal integration of the energy equation.

This project is coded as a C++ library interfaced with Python scripts, to benefit from the performance of the C++, the flexibility of Python and the graphical rendering of the Matplotlib library. The interface is realized with the PyBind11 library.

## Requirements
To compile the library :
- CMake  
- A C++ compiler  
- The PyBind11 library  

To generate the reference documentation :
- Doxygen  

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
user@linux:~/path/to/code/build$ make reference_doc
```

## How to use
Execute these commands to get to the python directory and add the library to the PYTHONPATH.
```console
user@linux:~/path/to/code$ cd python
user@linux:~/path/to/code/python$ export PYTHONPATH=$PWD/../build/install/lib/python
```

For each simulation, follow these steps :
- edit the parameters.py file to define your simulation
- execute the simulation.py file
- execute the plot.py file
