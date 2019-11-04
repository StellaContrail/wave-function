compile=gfortran -O3 -fcheck=all
gpe: gpe.o
	$(compile) -o gpe *.o
setting.o: setting.f90
	$(compile) -c setting.f90
mathf.o: mathf.f90
	$(compile) -c mathf.f90
gpe.o: gpe.f90 mathf.o setting.o
	$(compile) -c gpe.f90
clean:
	rm -f *.o *.mod