OBJGROUP = ../../modules/fem_module.f90 FRPGenerate_module.f90 mt19937.f 2d_frp_damage.f90
EXE = 2d_frp_damage.exe

all: clean femsolver

femsolver: $(OBJGROUP)
	gfortran -o $(EXE)  $(OBJGROUP)
	del *.mod

clean:
	del $(EXE)