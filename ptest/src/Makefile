FC=ifort
FFLAGS=-132 -mkl -mcmodel=medium

.SUFFIXES:    .f .for .o

all:gm2gid

GOBJS=gmesh2gid.o bound.o

.for.o:
	$(FC) $(FFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

gm2gid:$(GOBJS)
	$(FC) $(FFLAGS) -o $@ $(GOBJS)

clean:
	rm -f *.o gm2gid
