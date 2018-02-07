ciwi: interp.o ciwi.f90
	gfortran -o ciwi interp.o ciwi.f90

interp.o: interp.f90
	gfortran -c interp.f90

clean: 
	rm *.o *.mod ciwi