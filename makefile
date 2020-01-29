wardle: wardle.f collfile.o network.o 
	gfortran -o wardle wardle.f90 collfile.o network.o 

network.o: network.f90
	gfortran -c network.f90

collfile.o: collfile.f
	gfortran -c collfile.f

clean:
	rm ordout.d output.d wardle