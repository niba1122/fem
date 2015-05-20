all: clean femsolver

femsolver: ./src/compile_mac.f90
	gfortran -o ./bin/compile.out ./src/compile_mac.f90

clean:
	find ./bin/ -name "*.out" -exec rm {} \;
