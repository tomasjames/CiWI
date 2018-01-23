ciwi: interpolate.o ciwi.f90
	gfortran -o ciwi interpolate.o ciwi.f90

interpolate.o: interpolate.f90
	gfortran -c interpolate.f90

clean: 
	rm *.o *.mod ciwi