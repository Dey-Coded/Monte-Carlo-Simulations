# Monte-Carlo-Simulations
Monte Carlo simulation of 2D-magnetic materials
Read the project report for the detailed description and results of the simulations

# Compilation and execution
The simulation codes: ising.cpp and xy.cpp are parallelized using mpi while using the c++ wrappers of the python library matplotlib. So in order to compile them, use the following commands:

```shell
mpic++ ising.cpp -std=c++11 -I/usr/include/python2.7 -lpython2.7 -o ising
```
```shell
mpic++ xy.cpp -std=c++11 -I/usr/include/python2.7 -lpython2.7 -o xy
```
To execute the output files on 4 processes use the following commands:

```shell
mpiexec -n 4 ./ising
```
```shell
mpiexec -n 4 ./xy
```


