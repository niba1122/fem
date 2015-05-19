OBJGROUP = GFRPMeshGenerator2D_1_1.f90 PeriodicMeshGenerator2D_1_3.f90
EXE = PeriodicMeshGenerator2D.exe

all: clean meshcopytest

meshcopytest: $(OBJGROUP)
	gfortran -o $(EXE)  $(OBJGROUP)
	del *.mod

clean:
	del $(EXE)