#FC = gfortran
FC=ifort
FFLAGS = -O -Wall -fcheck=all -g -fbacktrace -g -ffpe-trap=zero,invalid,overflow,underflow
# FFLAGS = -g -traceback -check all -fp-stack-check -debug all -warn -CB

wardle: wardle.f90 dlsode.o
	${FC} ${FFLAGS} -o wardle wardle.f90 dlsode.o

collfile.o: collfile.f
	${FC} ${FFLAGS} collfile.f

dlsode.o: dlsode.f90
	${FC} ${FFLAGS} -c dlsode.f90

clean:
	rm ordout.d output.d wardle *.o *.mod

ciwi: interp.o params.o ciwi.f90
	${FC} ${FFLAGS} -o ciwi interp.o params.o ciwi.f90

interp.o: interp.f90
	${FC} ${FFLAGS} -c interp.f90

params.o: params.f90
	${FC} ${FFLAGS} -c params.f90

clean-ciwi: 
	rm *.o *.mod ciwi

network: network.f90
	gfortran -o network network.f90 -ffpe-trap=invalid -fbacktrace

clean-network:
	rm network