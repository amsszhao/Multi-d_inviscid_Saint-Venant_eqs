# Multidimensional Stability and Transverse Bifurcation of Hydraulic Shocks and Roll Waves in Open Channel Flow

Matlab codes for Evans function solvers and Python codes for numerical time evolutions of 2d- hydraulic shocks and roll waves of inviscid Saint-Venant equations 


_Authors:_ Zhao Yang and Kevin Zumbrun 

For questions/comments please contact either the first author via yangzhao@amss.ac.cn or the second author via kzumbrun@indiana.edu

## Prerequisites

### Programs

*  Python 
*  Matlab 

### Python libraries

* clawpack https://www.clawpack.org/

## Description 

### Evans solvers

Inside the folder Evans solvers and its subfolders, one can find codes for Evans/Evans-Lopatinsky solvers for smooth hydraulic shocks and discontinuous hydraulic shocks and a periodic Evans-Lopatinsky solver for roll waves, respectively.

### time-evolution solvers

Inside the folder time-evolution solvers and its subfolders, one can find sample python codes for generating time-evolution raw data files for the 2d inviscid Saint-Venant equations with either dam-break or roll waves initial data. One can also find matlab codes that make use of the raw data files to create movies.

Movies are downloadable from https://doi.org/10.5281/zenodo.8199967.

### data

Inside the folder data, one can find boundaries.mat which contains raw datum of various stability boundaries of the 2d inviscid roll waves.







