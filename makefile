FC = gfortran
#FC = ifort
FFLAGS = -O2 -fPIC -g
#FFLAGS = -g -traceback -check all -fp-stack-check -debug all -warn -CB

wardle: wardle.f90 dvode.o
	${FC} ${FFLAGS} -o wardle wardle.f90 dvode.o

collfile.o: collfile.f
	${FC} ${FFLAGS} collfile.f

dlsode.o: dlsode.f90
	${FC} ${FFLAGS} -c dlsode.f90

dvode.o: dvode.f90
	${FC} ${FFLAGS} -c dvode.f90

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

network: networkConversion.f90
	gfortran -o network networkConversion.f90 -ffpe-trap=invalid -fbacktrace

clean-network:
	rm network
