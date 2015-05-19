OBJGROUP = ../../modules/fem_module.f90 linear_solver.f90
EXE = linear_solver.exe

all: clean femsolver

femsolver: $(OBJGROUP)
	gfortran -o $(EXE)  $(OBJGROUP)
	del *.mod

clean:
	del $(EXE)