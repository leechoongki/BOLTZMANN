.SUFFIXES : .f90 .o

#MODULE_LIB=modules.o minresqlpDataModule.o minresqlpBlasModule.o minresqlpModule.o
MODULE_LIB=modules.o minresqlpDataModule.o  zminresqlpDataModule.o minresqlpBlasModule.o zminresqlpBlasModule.o zminresqlpModule.o
OBJS= $(MODULE_LIB) main.o genseed.o function.o mcstep.o twop.o eloc.o sr.o

FC=gfortran
LK=gfortran
FFLAGS= -c -O3
vmc.x:$(OBJS)
	$(LK) -o vmc.x $(OBJS)
.f90.o:
	$(FC) $(FFLAGS) $<
clean:
	rm -f $(OBJS) *.mod
