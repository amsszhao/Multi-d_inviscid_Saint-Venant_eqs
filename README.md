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

### movies

Inside the folder movies, one can find sample numerical time-evolutions with perturbed dam-break initial data in the sub-folder hydraulic shocks and perturb roll waves initial data in the sub-folder roll waves, respectively.

### Evans solvers

Inside the folder Evans solvers, one can find codes for Evans/Evans-Lopatinsky solvers in the sub-folder hydraulic shocks and a periodic Evans-Lopatinsky solver in the sub-folder roll waves, respectively.

### time-evolution solvers

Inside the folder time-evolution solvers, one can find sample python codes for generating time-evolutions of the 2d inviscid Saint-Venant equations with either dam-break or roll waves initial data. One can also find matlab codes that make use of the raw data files to create movies.

### data

Inside the folder data, one can find boundaries.mat which contains raw datum of various stability boundaries of the 2d inviscid roll waves.







