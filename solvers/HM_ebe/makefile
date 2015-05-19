OBJGROUP = ../../modules/fem_module.f90 ../../modules/hm_module.f90 ../../modules/ebe_scg_module.f90 hm_ebe.f90
EXE = hm_ebe.exe

all: clean femsolver

femsolver: $(OBJGROUP)
	gfortran -fopenmp -o $(EXE)  $(OBJGROUP)
	del *.mod

clean:
	del $(EXE)