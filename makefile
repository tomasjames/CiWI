wardle: wardle.f collfile.o
	gfortran -o wardle wardle.f collfile.o

collfile.o: collfile.f
	gfortran -c collfile.f

clean:
	rm ordout.d output.d wardle

ciwi: interp.o params.o ciwi.f90
	gfortran -o ciwi interp.o params.o ciwi.f90

interp.o: interp.f90
	gfortran -c interp.f90

params.o: params.f90
	gfortran -c params.f90

clean-ciwi: 
	rm *.o *.mod ciwi