MkInclude = ${R_HOME}/etc${R_ARCH}/Makeconf

IEXTLIBS = -ISuiteSparse_config -ICAMD/Include -ICHOLMOD/Include
LEXTLIBS =	
E_LIBS = CHOLMOD.a CAMD.a AMD.a COLAMD.a SuiteSparse_config.a 
SUBDIRS = CHOLMOD COLAMD AMD CAMD SuiteSparse_config

PKG_CFLAGS = $(SHLIB_OPENMP_CFLAGS) $(IEXTLIBS)
PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS) $(IEXTLIBS)
PKG_LD = $(CXX)
PKG_LDFLAGS = $(ALL_CXXFLAGS) 
ALL_LIBS = $(E_LIBS) ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} -lm

all : excursions gaussint

EXCURSIONS_OBJ = excursions.o integration.o utils.o RngStream.o \
			gsl_fix.o xerbla.o
excursions : $(EXCURSIONS_OBJ) $(E_LIBS)
	$(PKG_LD) $(MAIN_LDFLAGS) $(LDFLAGS) $(PKG_LDFLAGS) $(LEXTLIBS) -o $@ $(EXCURSIONS_OBJ) $(ALL_LIBS)

GAUSSINT_OBJ = gaussint.o integration.o utils.o RngStream.o gsl_fix.o xerbla.o
gaussint : $(GAUSSINT_OBJ) $(E_LIBS)
	$(PKG_LD) $(MAIN_LDFLAGS) $(LDFLAGS) $(PKG_LDFLAGS) $(LEXTLIBS) -o $@ $(GAUSSINT_OBJ) $(ALL_LIBS) 

CHOLMOD.a:
	(cd CHOLMOD && CFLAGS="$(CFLAGS)" CXXFLAGS="$(CXXFLAGS)" \
		MkInclude="$(MkInclude)" $(MAKE) library)

CAMD.a:
	(cd CAMD && CFLAGS="$(CFLAGS)" CXXFLAGS="$(CXXFLAGS)" \
		MkInclude="$(MkInclude)" $(MAKE) library)

AMD.a:
	(cd AMD && CFLAGS="$(CFLAGS)" CXXFLAGS="$(CXXFLAGS)" \
		MkInclude="$(MkInclude)" $(MAKE) library)

COLAMD.a:
	(cd COLAMD && CFLAGS="$(CFLAGS)" CXXFLAGS="$(CXXFLAGS)" \
		MkInclude="$(MkInclude)" $(MAKE) library)

SuiteSparse_config.a:
	(cd SuiteSparse_config && CFLAGS="$(CFLAGS)" CXXFLAGS="$(CXXFLAGS)" \
		MkInclude="$(MkInclude)" $(MAKE) library)

clean: subclean
	@-rm -f *.o excursions gaussint

subclean:
	@-rm -f *.a
	@for d in $(SUBDIRS); do \
	  (cd $${d} && MkInclude="" $(MAKE) clean) || exit 0; \
	done

.PHONY: all clean subclean
