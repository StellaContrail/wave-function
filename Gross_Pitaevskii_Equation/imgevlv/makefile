compiler=gfortran
options=-O3 -fbounds-check -llapack95 -llapack -lblas
all:mathf.mod io.mod setting.mod gpe.mod a.out

mathf.mod mathf.o: mathf.f90
	gfortran -c mathf.f90 $(options)

io.mod io.o: io.f90
	gfortran -c io.f90 $(options)

setting.mod setting.o: setting.f90
	gfortran -c setting.f90 $(options)

gpe.mod gpe.o: gpe.f90
	gfortran -c gpe.f90 $(options)

a.out: gpe.o io.o setting.o mathf.o
	gfortran -o a.out *.o $(options)

clean:
	rm *.o *.mod a.out