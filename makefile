wardle: wardle.f collfile.o
	gfortran -o wardle wardle.f90 collfile.o

collfile.o: collfile.f
	gfortran -c collfile.f

clean:
	rm ordout.d output.d wardle