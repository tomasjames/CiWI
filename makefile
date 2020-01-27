SHELL = /bin/sh
FC =gfortran
FFLAGS = -O2 -fPIC

wardle: wardle.f collfile.o 
	${FC} ${FFLAGS} -o wardle wardle.f90 collfile.o 

collfile.o: collfile.f
	${FC} ${FFLAGS} -c collfile.f

clean:
	rm *.o *.mod ordout.d output.d wardle