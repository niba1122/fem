OBJGROUP = ../../modules/fem_module.f90 ../../modules/rotateD_module.f90 linear_solver_rot.f90
EXE = linear_solver_rot.exe

all: clean femsolver

femsolver: $(OBJGROUP)
	gfortran -o $(EXE)  $(OBJGROUP)
	del *.mod

clean:
	del $(EXE)