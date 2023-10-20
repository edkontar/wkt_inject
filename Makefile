#absoft Fortran
# CC = af95 
# CFLAGS = -Ofast -m64
# gfortran
# CC = f95
# Intel Fortran works faster
CC = ifort
CFLAGS = -O3
COPTIONS = 
FILE = wkt
OBJS = constant.o params.o nonlins.o reader.o writer.o initbeam.o 1beam.o

prog: $(OBJS)
	@echo "linking..."
	$(CC) $(CFLAGS) -o$(FILE) $(OBJS) 
	make clean


$(OBJS): 
	$(CC) -I $(CFLAGS) $(OPTIONS) -c $*.f90

clean:
	@echo "cleaning up..."
	rm -f *.o
	rm -f *~
	rm -f *.mod
