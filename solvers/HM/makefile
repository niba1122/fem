OBJGROUP = ../../modules/fem_module.f90 ../../modules/hm_module.f90 hm.f90
EXE = hm.exe

all: clean femsolver

femsolver: $(OBJGROUP)
	gfortran -o $(EXE)  $(OBJGROUP)
	del *.mod

clean:
	del $(EXE)