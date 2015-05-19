OBJGROUP = ../../modules/fem_module.f90 ../../modules/ebe_scg_module.f90 ebe_scg.f90
EXE = ebe_scg.exe

all: clean femsolver

femsolver: $(OBJGROUP)
	gfortran -fopenmp -o $(EXE)  $(OBJGROUP)
	del *.mod

clean:
	del $(EXE)