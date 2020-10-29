OOQP=./OOQP-0.99.26

OOQPINCLUDEDIR=$(OOQP)/include
OOQPLIBDIR=$(OOQP)/lib
GFORTRANLIB=/usr/local/gfortran/lib

CFLAGS = $(shell gsl-config --cflags)
CFLAGS += -O3 -I include -I /include
CFLAGS += -I$(OOQPINCLUDEDIR)

#DIRECTORY TO SEARCH FOR LIBS
LDFLAGS  =-L$(OOQPLIBDIR) -L$(GFORTRANLIB)
MA27LIB  = -L$(OOQP)/

#For MAC users, specify directories for quadmath, stdc++, libm, and gfortran--not needed!
#GCCLIB = /usr/local/gfortran/lib/gcc/x86_64-apple-darwin14/4.9.2
#GFLIB = /usr/local/gfortran/lib

#For LINUX users, specify directories for blas (or g2c), stdc++, libm, and gfortran. The numbered versions here will need to be updated.
LINGCC = /usr/lib/gcc/x86_x64-redhat-linux/4.1.2
GFLIB =
BLAS = 

#Replace -lblas with lg2c if necessary
FLIBSLIN    =  -L$(LINGCC) -L$(BLAS) -L$(GFLIB) -lgfortran  -lm -lstdc++ -lblas -lg2c

FLIBSMAC    =  -L/opt/local/lib/gcc10/gcc/x86_64-apple-darwin18/10.1.0 -L/opt/local/lib/gcc10/gcc/x86_64-apple-darwin18/10.1.0/../../.. -lgfortran -lquadmath -lm -lstdc++ -framework accelerate
#-lgfortran -lquadmath -lm -lstdc++ -framework accelerate


LIBS   = $(shell gsl-config --libs)
LIBS  += $(MA27LIB) 

LIBS2  = -looqpgensparse -looqpsparse  -looqpgondzio -looqpbase -lma27 $(LIBS) 

CC     = g++
BDIR   = bin
ODIR   = obj
SDIR   = subroutines
EDIR   = EXES
OS    := $(shell uname)

ifeq ($(OS),Linux)
	_OBJS = $(shell ls $(SDIR)/*.cpp | xargs -n 1 basename | sed -r 's/(\.cc|.cpp)/.o/')
	LIBS  +=$(FLIBSLIN)
else
	_OBJS = $(shell ls $(SDIR)/*.cpp | xargs -n 1 basename | sed -E 's/(\.cc|.cpp)/.o/')
	LIBS  +=$(FLIBSMAC)
endif

_EXECUTABLES = 	stoich_balance_knockdowns.exe \
##		stoich_balance.exe \



##look to find all subroutines
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))



EXECUTABLES = $(patsubst %,$(BDIR)/%,$(_EXECUTABLES))



_SOURCES = ${_EXECUTABLES:=.cpp}
SOURCES = $(patsubst %,$(EDIR)/%,$(_SOURCES))

all: dirs $(EXECUTABLES)

$(ODIR)/%.o: $(SDIR)/%.cpp
	@echo "Compiling $<"
	$(CC) $(CFLAGS) -c $< -o $@  
	@echo "------------"


$(EXECUTABLES): $(OBJS)
	@echo "Compiling $(EDIR)/$(@F:.exe=.cpp)"
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(EDIR)/$(@F:.exe=.cpp) $(OBJS) $(LIBS2)
	@echo "------------"

dirs: 
	mkdir -p bin
	mkdir -p obj

clean:
	rm -f $(ODIR)/*.o $(EXECUTABLES)


