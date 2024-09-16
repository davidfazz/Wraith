# makefile - elapid - david gallina - university of kassel - 2022
#...............................................................................

# compiler
FC = gfortran
#FC = ifort
# flags
#FFLAGS = -O2 -C -fbounds-check -fopenmp
FFLAGS = -O3 -C -fopenmp -framework Accelerate
VPATH=pots:algs

# objects in source folder

OBJS = constants.o functions.o random.o main.o potential.o \
       lbfgs.o explore.o dneb.o efol.o rates.o heatcap.o mc.o \
	   dpls1d.o dpls2d.o unsurf.o chains.o dneb.o follow.o hysteresis.o \
	   open1d.o

LIBS = #liblapack.a libblas.a

EXEC = ../bin
 
# Rules and Targets
.SUFFIXES: 
.SUFFIXES: .o .f .f90

.f90.o:
	$(FC) $(FFLAGS) -c  $<
.f.o:
	$(FC) $(FFLAGS) -c  $<
	
$(EXEC)/chara: $(OBJS)
	$(FC) $(FFLAGS) -o $(EXEC)/xc *.o $(LIBS)

clean: 
	rm -f $(OBJS) *.o *.mod
	
# dependencies
$(OBJS):     	    $(LIBS)
main.o:         constants.o functions.o random.o potential.o explore.o follow.o \
	            hysteresis.o rates.o heatcap.o mc.o
explore.o:		constants.o functions.o random.o potential.o lbfgs.o efol.o dneb.o
dneb.o:  	  	constants.o functions.o potential.o  
efol.o:      	constants.o functions.o potential.o
lbfgs.o:     	constants.o functions.o potential.o
follow.o:		constants.o functions.o potential.o dneb.o lbfgs.o
hysteresis.o:   constants.o functions.o potential.o
potential.o:    constants.o functions.o dpls1d.o dpls2d.o unsurf.o chains.o open1d.o
dpls1d.o:    	constants.o functions.o random.o
dpls2d.o:    	constants.o functions.o random.o
open1d.o:		constants.o functions.o random.o
unsurf.o:	  	constants.o functions.o random.o
chains.o:	  	constants.o functions.o random.o
rates.o:		constants.o functions.o potential.o
mc.o:			constants.o functions.o potential.o random.o
heatcap.o: 		constants.o functions.o
random.o:      	constants.o
functions.o:    constants.o



