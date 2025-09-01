# drop_pool_coalescence
This repository contain the codes (written in C++) for computation of coalescence dynamics of a drop on a horizontal liquid pool using the geometric PLIC-CLSVOF method. The flow is considered to be axisymmetric in nature. The governing incompressible Navier-Stokes equation is solved in a mass-momentum consistent framework. Further, the finite volume method is used to discretize the governing equations in a co-located structured mesh. An additional volume fraction field is used to track the original drop liquid during the coalescence process; it acts as a passive tracer field which do not participate in the solution process.

The simulation results are qualitatively presented as an animation video:

PROVIDE LINK

# Software requirements
This solver needs:

- gcc

# How to install the required packages (on a Linux system)

To install gcc (compiler for C++ code)

```bash
sudo apt install build-essential
```

# How to compile and run the code

To compile the code

```bash
g++ driver.cpp -o output
```
To run this code

```bash
./output
```
