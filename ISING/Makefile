.SUFFIXES : .f90 .o

OBJS=modules.o main.o genseed.o function.o mcstep.o twop.o eloc.o

FC=gfortran
LK=gfortran
FFLAGS= -c -O3
vmc.x:$(OBJS)
	$(LK) -o vmc.x $(OBJS)
.f90.o:
	$(FC) $(FFLAGS) $<
clean:
	rm -f $(OBJS)
