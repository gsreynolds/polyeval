MF=	Makefile

FC=	gfortran

FFLAGS= -O3 -fopenmp

LFLAGS=

EXE= PolyEval

SRC= \
	modPolyEval.f90 \
	PolyEval_test.f90

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)
MOD=	$(SRC:.f90=.mod)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF)

clean:
	rm -f $(OBJ) $(EXE) $(MOD) core
