ciwi: interp.o params.0 ciwi.f90
	gfortran -o ciwi interp.o params.o ciwi.f90

interp.o: interp.f90
	gfortran -c interp.f90

params.0: params.f90
	gfortran -c params.f90

clean: 
	rm *.o *.mod ciwi