OBJGROUP = ../../modules/fem_module.f90 MeshTester.f90
EXE = MeshTester.exe

all: clean meshtester

meshtester: $(OBJGROUP)
	gfortran -o $(EXE)  $(OBJGROUP)
	del *.mod

clean:
	del $(EXE)