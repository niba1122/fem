OBJGROUP = ../../modules/fem_module.f90 ../../modules/rotateD_module.f90 ../../modules/FRPGenerate_module.f90 mc_frp.f90 mt19937.f
EXE = mc_frp.exe

all: clean femsolver

femsolver: $(OBJGROUP)
	gfortran -o $(EXE)  $(OBJGROUP)
	del *.mod

clean:
	del $(EXE)