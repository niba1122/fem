# FEM Solver 

## Description

This solver is used in Takano Laboratory in Keio University.
This solver is written in fortran 2003 and assumed to be compiled by gfortran.

## How to use
### Installation
1. Download or `git clone ~` from GitHub.
2. Execute `make` or `make all` at 'fem/'.
3. Add 'export PATH="$HOME/fem/bin:$PATH"' to ~/.bash_profile

### Compile solver
Execute `fem compile SOLVER_NAME`

### Start calculation
Execute `fem exec SOLVER_NAME`

