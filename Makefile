#The path to ILOG/CPLEX
ILOGBASE = /opt/ibm/ILOG/CPLEX_Studio1263

#ILOG/CPLEX subdirectory and filename for statically linked library
ILOGLIBFORMAT = x86-64_linux/static_pic
#BLAS implementation to use with Armadillo
BLASLIB = blas

ILOGINCLUDES = -DILOGLUE -DIL_STD -I$(ILOGBASE)/cpoptimizer/include/ -I$(ILOGBASE)/concert/include/ -I$(ILOGBASE)/cplex/include/
ILOGLIBS  = -L$(ILOGBASE)/cpoptimizer/lib/$(ILOGLIBFORMAT) -L$(ILOGBASE)/concert/lib/$(ILOGLIBFORMAT) -L$(ILOGBASE)/cplex/lib/$(ILOGLIBFORMAT) -lcp -lilocplex -lcplex -lconcert -lm

ARMAINCLUDES = -I~/usr/include
ARMALIBS = -L~/usr/lib -larmadillo -l$(BLASLIB)

IDIR=.

CC=gcc
CXX=g++
CFLAGS=-fopenmp -g -Wall -w -O3 -I$(IDIR) $(ILOGINCLUDES) $(ARMAINCLUDES)

ODIR=.

LIBS=-L. $(ILOGLIBS) $(ARMALIBS)

_DEPS = AgileCombiFD.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = AgileCombiFD.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all:     main

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

main: $(OBJ) $(ODIR)/main.o
	$(CXX) -o agilefd $(OBJ) $(ODIR)/main.o $(CFLAGS) $(LIBS)


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ agilefd
