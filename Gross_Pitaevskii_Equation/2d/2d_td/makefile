compiler=gfortran
options=-O3 -fbounds-check
all:mathf.mod io.mod setting.mod gpe.mod a.out constants.mod

constants.mod: constants.f90
	gfortran -c constants.f90 $(options)
	
mathf.mod mathf.o: mathf.f90 constants.mod
	gfortran -c mathf.f90 $(options)

io.mod io.o: io.f90 constants.mod
	gfortran -c io.f90 $(options)

setting.mod setting.o: setting.f90 constants.mod
	gfortran -c setting.f90 $(options)

gpe.mod gpe.o: gpe.f90 constants.mod
	gfortran -c gpe.f90 $(options)

a.out: gpe.o io.o setting.o mathf.o
	gfortran -o a.out *.o $(options)
	rm *.o *.mod 
clean:
	rm *.o *.mod a.out