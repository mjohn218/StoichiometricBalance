all:
all-depend: depend all

SHELL =    /bin/sh

CXX      = CC
F77      = gfortran
FC       = $(F77)
CXXFLAGS = -O
CC       = cc
CFLAGS   = -O
FFLAGS   = -O
INSTALL  = /usr/bin/install -c
CPPFLAGS =  $(INCLUDES) $(DEFS)
DEFS     = -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE_URL=\"\" -DSTDC_HEADERS=1 -DHAVE_GETRUSAGE=1
LDFLAGS  = 
FLIBS    =  -L/opt/local/lib/gcc10/gcc/x86_64-apple-darwin18/10.1.0 -L/opt/local/lib/gcc10/gcc/x86_64-apple-darwin18/10.1.0/../../.. -lgfortran -lquadmath -lm -lstdc++
BLAS     = -lblas -llapack
LIBS     = 
RANLIB   = ranlib
MAKEDEPEND = makedepend

SUPERLU     = /sandbox/SuperLU
SUPERLULIB  = $(SUPERLU)/superluLinux-i386.a
OBLIOLIB    = /sandbox/oblio/lib/oblio.a
MA57LIB     = libma57.a libfakemetis.a
MA27LIB     = libma27.a

srcdir      = src
prefix      = /usr/local
includedir  = ${prefix}/include
libdir      = ${exec_prefix}/lib
bindir      = ${exec_prefix}/bin
datadir     = ${prefix}/share
docdir      = $(datadir)/doc/packages/OOQP
exec_prefix = ${prefix}
locallibdir = ./lib
localincludedir = ./include

SOURCES =
PROGRAMS =
INCLUDES :=

#include all modules (subdirectories directly below this one that have
#a makefile.inc.)

MODULEMAKEFILES  := $(wildcard $(srcdir)/*/Makefile.inc)
MODULES          := $(dir $(MODULEMAKEFILES))
INCLUDES         :=  $(patsubst %, -I%, $(MODULES))
SOURCES   := $(foreach MODULE, $(MODULES), $(wildcard $(MODULE)*.[cC])) 

include $(MODULEMAKEFILES)
include $(wildcard $(srcdir)/*/MakefileTargets.inc)

all: $(PROGRAMS) all_headers all_libs
all_ma57:

# No user-serviceable parts below here
clean mostlyclean: 
	-find . -name "*.o" -exec rm \{\} \;

veryclean: clean headers_veryclean libs_veryclean
	-find . -name "*.exe" -exec rm \{\} \;

headers_veryclean:
	-rm -f $(localincludedir)/*.h

libs_veryclean:
	-rm -f $(locallibdir)/*.a

distclean: veryclean
	-rm -f config.status config.log config.cache TAGS
	-find . -name "*.out" -exec rm \{\} \;
	-find . -name "*.bak" -exec rm \{\} \;
	-find . -name "*~" -exec rm \{\} \;
	-if [ -d doc-src ]; then cd doc-src; $(MAKE) distclean; fi
	-rm  GNUmakefile PetscMakefile

maintainer-clean: distclean
	-if [ -d doc-src ]; then cd doc-src; $(MAKE) maintainer-clean; fi
	-rm -f configure

depend:
	$(MAKEDEPEND) -f GNUmakefile $(INCLUDES) $(CPPFLAGS) \
		$(DEFS) $(SOURCES) > /dev/null 2>&1
	@echo " *** Done making depend"

install: all install_headers install_libs install_docs install_exe
install_headers: all_headers
	$(INSTALL) -d $(includedir)/ooqp	
	$(INSTALL) $(localincludedir)/*.h $(includedir)/ooqp

install_libs: all_libs
	$(INSTALL) -d $(libdir)
	$(INSTALL) $(locallibdir)/*.a $(libdir)

install_exe:
	$(INSTALL) -d $(bindir)
	$(INSTALL) *.exe $(bindir)

install_docs:
	if [ -d doc-src ]; then cd doc-src; fi; make
	$(INSTALL) -d $(docdir)
	$(INSTALL) doc/*.pdf doc/*.html $(docdir)
	$(INSTALL) -d $(docdir)/reference-manual
	$(INSTALL) doc/reference-manual/*.html \
		doc/reference-manual/*.png \
		doc/reference-manual/*.css \
		$(docdir)/reference-manual
	$(INSTALL) -d $(docdir)/formulations
	$(INSTALL) doc/formulations/*.html $(docdir)/formulations
	$(INSTALL) -d $(docdir)/distribution-docs
	$(INSTALL) doc/distribution-docs/*.html $(docdir)/distribution-docs 

GNUmakefile: GNUmakefile.in config.status
	./config.status

config.status: configure
	./config.status --recheck

rebuild-manifest:
	cd doc-src; make
	make distclean
	find * .cvsignore \
	\( -path "*CVS" -prune \) -o \( -path doc-src -prune \) \
	-o -type f -print > MANIFEST

rebuild-VERSION:
	make qpgen-sparse-gondzio.exe
	./qpgen-sparse-gondzio.exe --version > VERSION

.SUFFIXES:
.SUFFIXES: .o .c .cc .C .f


