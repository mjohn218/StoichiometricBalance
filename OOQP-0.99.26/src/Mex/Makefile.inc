MEX          = 
MEXSFX       = mexmex
MEX_ENABLED  = no # MATLAB support has been disabled

OOQPMEXDIR = $(srcdir)/Mex

MEXUTILITYSRC = $(OOQPMEXDIR)/mexUtility.c
MEXUTILITYOBJ = $(OOQPMEXDIR)/mexUtility.o

MATLABFUNCTIONS = readmps.m ooqp_svm.m ooqp_huber.m ooqp.m ooqp_start.m
MEXFUNCTIONS = ooqp_mex.$(MEXSFX) readmps_mex.$(MEXSFX) \
	ooqp_svm.$(MEXSFX) ooqp_huber.$(MEXSFX) ooqp_start_mex.$(MEXSFX)

ooqp_mex_OBJ = $(CQPGENSPARSEOBJ) \
		$(OOQPMEXDIR)/mexUtility.o $(libooqpgondzio_STATIC) \
		$(libooqpgensparse_STATIC) \
		$(libooqpsparse_STATIC) $(libooqpbase_STATIC) \

ooqp_start_mex_OBJ = 	        $(CQPGENSPARSEOBJ) \
		$(OOQPMEXDIR)/mexUtility.o $(OOQPMEXDIR)/ooqp_start_cpp.o \
		$(libooqpgondzio_STATIC) \
		$(libooqpgensparse_STATIC) \
		$(libooqpsparse_STATIC) $(libooqpbase_STATIC) \

readmps_mex_OBJ = $(OOQPMEXDIR)/mexUtility.o $(MPSREADEROBJ) \
		 	$(libooqpbase_STATIC)

ooqp_svm_mex_OBJ = $(OOQPMEXDIR)/mexUtility.o $(SVMOBJECTS) \
	     	$(GONDZIOOBJECTS) $(ABOBJECTS) $(UTILOBJECTS) \
	     	$(DELINALGOBJECTS) $(VECTOROBJECTS) \
	     	$(OOQPMEXDIR)/cfunc_ooqp_svm.o

ooqp_huber_mex_OBJ =  $(OOQPMEXDIR)/mexUtility.o $(HUBEROBJECTS) \
	$(GONDZIOOBJECTS) $(ABOBJECTS) \
	$(VECTOROBJECTS) $(DELINALGOBJECTS) $(UTILOBJECTS) \
	$(OOQPMEXDIR)/cfunc_ooqp_huber.o
