#
# Makefile for XNet and associated programs
#
# Build search path
VPATH = .

# System configurations
# Architecture and installation specific MACROS
include Makefile_local

## Build configurations using variables defined in System Makefile
#CMODE = OPTIMIZE
CMODE = DEBUG
#CMODE = DEBUGOPT
#CMODE = PROFILE

FFLAGS   = $(R8) 
LDRFLAGS = $(R8)
CCFLAGS  = 
CXXFLAGS =

## Choose matrix solver package, "build" compiles LAPACK, others use existing libraries 
#MATRIX_SOLVER = PARDISO
#MATRIX_SOLVER = MA28
#MATRIX_SOLVER = MA48
#MATRIX_SOLVER = gpu
MATRIX_SOLVER = build
MATRIX_SOLVER ?= dense

## Choose EOS for screening; HELMHOLTZ (requires separate download) or BAHCALL
#EOS = HELMHOLTZ
EOS ?= BAHCALL

## Choose executable name
EXE = xnet

## SELECT MPI MODE, ON=MPI, OFF= noMPI 
#MPI_MODE = ON

## SELECT OPENMP MODE, Default=OFF
#OPENMP_MODE = ON

## SELECT GPU MODE, Default=OFF
#GPU_MODE = ON

## Defaults
MPI_MODE    ?= OFF
OPENMP_MODE ?= OFF
GPU_MODE    ?= OFF

# Optional features for XNet, uncomment to enable
#NSE_OBJ = nse.o 		# Also uncomment lines labeled !NSE or use sed -e 's/^\!\(.*\!NSE\)/ \1/' *.f90
			#        Comment lines labeled !NSE with   sed -e 's/^[^\!]\(.*\!NSE\)/!\1/' *.f90
#FFN_OBJ = ffn.o 		# Also uncomment lines labeled !FFN or use sed -e 's/^\!\(.*\!FFN\)/ \1/' *.f90
			#        Comment lines labeled !FFN with   sed -e 's/^[^\!]\(.*\!FFN\)/!\1/' *.f90
#NNU_OBJ = nu_nucleus.o 	# Also uncomment lines labeled !NNU or use sed -e 's/^\!\(.*\!NNU\)/ \1/' *.f90
			#        Comment lines labeled !NNU with   sed -e 's/^[^\!]\(.*\!NNU\)/!\1/' *.f90

#-------------------------------------------------------------------------------------------------
# The rules below are used by make to build a range of executables
# Users should not need to modify this portion of the Makefile
#

ifeq ($(CMODE),OPTIMIZE)
  FFLAGS   += $(OPT) $(EXTRA_FFLAGS)
  CCFLAGS  += $(OPT) $(EXTRA_CCFLAGS)
  CXXFLAGS += $(OPT) $(EXTRA_CXXFLAGS)
  LDRFLAGS += $(OPT) $(EXTRA_LDRFLAGS)
else ifeq ($(CMODE),DEBUG)
  FFLAGS   += $(DEBUG) $(EXTRA_FFLAGS)
  CCFLAGS  += $(DEBUG) $(EXTRA_CCFLAGS)
  CXXFLAGS += $(DEBUG) $(EXTRA_CXXFLAGS)
  LDRFLAGS += $(DEBUG) $(EXTRA_LDRFLAGS)
else ifeq ($(CMODE),DEBUGOPT)
  FFLAGS   += $(DEBUGOPT) $(EXTRA_FFLAGS)
  CCFLAGS  += $(DEBUGOPT) $(EXTRA_CCFLAGS)
  CXXFLAGS += $(DEBUGOPT) $(EXTRA_CXXFLAGS)
  LDRFLAGS += $(DEBUGOPT) $(EXTRA_LDRFLAGS)
else ifeq ($(CMODE),PROFILE)
  FFLAGS   += $(PROFILE) $(EXTRA_FFLAGS)
  CCFLAGS  += $(PROFILE) $(EXTRA_CCFLAGS)
  CXXFLAGS += $(PROFILE) $(EXTRA_CXXFLAGS)
  LDRFLAGS += $(PROFILE) $(EXTRA_LDRFLAGS)
endif

## Common XNet files
BASE = control.o data.o conditions.o $(FFN_OBJ) $(NNU_OBJ) common.o net_preprocess.o $(NSE_OBJ) full_net.o solve_be.o solve_bd.o match.o flux.o screening.o
MPI_FILES = data_distribute_mpi.o jacobian_bcast_dense.o jacobian_bcast_PARDISO.o jacobian_bcast_MA28.o net_mpi.o

## Configure Libraries and files to link appropriate EOS
ifeq ($(EOS),BAHCALL)
  EOS_OBJ = eos_bahcall.o
else ifeq ($(EOS),HELMHOLTZ)
  EOS_OBJ = eos_helm.o helmholtz.o
  VPATH += $(HELMHOLTZ_PATH)
endif
BASE += $(EOS_OBJ)


## Configure Libraries and files to link appropriate matrix solver
LIBS = 
ifeq ($(MATRIX_SOLVER),MA28)
  EXTRA_OBJ  += MA28.o
else ifeq ($(MATRIX_SOLVER),MA48)
  EXTRA_OBJ  += MA48.o
else ifeq ($(MATRIX_SOLVER),dense)
  LIBS       += $(LAPACK_LIBDIR) $(LAPACK_LIBS)
  EXTRA_OBJ  += $(LAPACK_OBJ)
else ifeq ($(MATRIX_SOLVER),PARDISO)
  LIBS       += $(PARDISO_LIBDIR) $(PARDISO_LIBS)
  EXTRA_OBJ  += $(PARDISO_OBJ)
else ifeq ($(MATRIX_SOLVER),gpu)
  LIBS       += $(LAPACK_LIBDIR) $(LAPACK_LIBS)
  EXTRA_OBJ  += $(LAPACK_OBJ)
  LIBS       += $(GPU_LAPACK_LIBDIR) $(GPU_LAPACK_LIBS)
  EXTRA_OBJ  += $(GPU_LAPACK_OBJ)
  LAPACK_INC += $(GPU_LAPACK_INC)
else ifeq ($(MATRIX_SOLVER),build)
  VPATH += ../tools/LAPACK
  BLAS_LVL3_SRC  = dgemm.o dtrsm.o
  BLAS_LVL2_SRC  = dger.o 
  BLAS_LVL1_SRC  = dlamch.o dscal.o dswap.o idamax.o lsame.o xerbla.o
  LAPACK_UTL_SRC = ieeeck.o iparmq.o
  LAPACK_AUX_SRC = dgetf2.o dlaswp.o ilaenv.o
  EXTRA_OBJ  +=  $(BLAS_LVL1_SRC) $(BLAS_LVL2_SRC) $(BLAS_LVL3_SRC)  $(LAPACK_UTL_SRC) $(LAPACK_AUX_SRC) dgetrf.o dgetrs.o dgesv.o
  MATRIX_SOLVER = dense
endif

# Define Driver
DRIVER_SM = net.o 
DRIVER_MPI = data_distribute_mpi.o net_mpi.o jacobian_bcast_$(MATRIX_SOLVER).o

# Configure flags for OpenMP
FFLAGS   += $(OMP_FLAG)
CCFLAGS  += $(OMP_FLAG)
CXXFLAGS += $(OMP_FLAG)
LDRFLAGS += $(OMP_FLAG)

# Configure flags for OpenACC
ifeq ($(GPU_MODE),ON)
  ifeq ($(MATRIX_SOLVER),dense)
    LIBS       += $(GPU_LAPACK_LIBDIR) $(GPU_LAPACK_LIBS)
    EXTRA_OBJ  += $(GPU_LAPACK_OBJ)
    LAPACK_INC += $(GPU_LAPACK_INC)
  endif
  ifeq ($(PE_ENV),PGI)
    ifeq ($(ACCMODE),OPENACC)
      CCFLAGS  += -acc -ta=nvidia
      CXXFLAGS += -acc -ta=nvidia
      FFLAGS   += -acc -ta=nvidia
      LDRFLAGS += -acc -ta=nvidia
    else
      CCFLAGS  += -ta=nvidia
      CXXFLAGS += -ta=nvidia
      FFLAGS   += -ta=nvidia
      LDRFLAGS += -ta=nvidia
    endif
  endif
else
  ifeq ($(PE_ENV),CRAY)
    CCFLAGS  += -h noacc
    CXXFLAGS += -h noacc
    FFLAGS   += -h noacc
    LDRFLAGS += -h noacc
  endif
endif

# Choose Driver and Compiler for MPI/nonMPI
ifeq ($(MPI_MODE),ON)
  DRIVER = $(DRIVER_MPI)
  FCL    = $(FC_MPI)
  CCL    = $(CC_MPI)
  CXXL   = $(CXX_MPI)
  LDRL   = $(LDR_MPI)
else
  DRIVER = $(DRIVER_SM)
  FCL    = $(FC)
  CCL    = $(CC)
  CXXL   = $(CXX)
  LDRL   = $(LDR)
endif

$(EXE): $(BASE) jacobian_$(MATRIX_SOLVER).o $(DRIVER) $(EXTRA_OBJ) 
	$(LDRL) $(LDRFLAGS) -o $(EXE) \
	    $(BASE) jacobian_$(MATRIX_SOLVER).o $(notdir $(EXTRA_OBJ)) $(DRIVER) \
	    $(EXTRA_LINK) $(LIBS)

#
# Rules to build specific configurations of XNet and supporting tools
#
xnet_dense: $(BASE) jacobian_dense.o $(DRIVER_SM) 
	$(LDR) $(LDRFLAGS) -o xnetd \
	    $(BASE) jacobian_dense.o $(notdir $(LAPACK_OBJ)) $(DRIVER_SM) \
	    $(LAPACK_LIBDIR) $(LAPACK_LIBS) $(EXTRA_LINK)

xnet_MA28: $(BASE) jacobian_MA28.o MA28.o $(EXTRA_OBJ)
	$(LDR) $(LDRFLAGS) -o xnetm \
	    $(BASE) jacobian_MA28.o $(DRIVER_SM) MA28.o \
	    $(EXTRA_LINK)

xnet_MA48: $(BASE) jacobian_MA48.o MA48.o $(EXTRA_OBJ)
	$(LDR) $(LDRFLAGS) -o xnetm \
	    $(BASE) jacobian_MA28.o $(DRIVER_SM) MA48.o \
	    $(EXTRA_LINK)

xnet_PARDISO: $(BASE) jacobian_PARDISO.o $(DRIVER) $(EXTRA_OBJ)
	$(LDR) $(LDRFLAGS) -o xnetp \
	    $(BASE) jacobian_PARDISO.o $(DRIVER_SM) \
	    $(PARDISO_INC) $(PARDISO_LIBDIR) $(PARDISO_LIBS) $(EXTRA_LINK)

xnet_nse: $(BASE) jacobian_$(MATRIX_SOLVER).o nse.o $(EXTRA_OBJ)
	sed -e 's/^\!\(.*\!NSE\)/ \1/' -e 's/^[^\!]\(.*\!NOTNSE\)/!\1/' net.f90 >| net_nse.f90
	$(LDR) $(LDRFLAGS) -o xnet_nse \
	    $(BASE) jacobian_$(MATRIX_SOLVER).o nse.o net_nse.f90 \
	    $(LIBS) $(EXTRA_LINK)

xnet_MPI: $(BASE) jacobian_$(MATRIX_SOLVER).o $(DRIVER_MPI) $(EXTRA_OBJ)
	$(LDR_MPI) $(LDRFLAGS) -o xnet_mpi \
	    $(BASE) jacobian_$(MATRIX_SOLVER).o $(DRIVER_MPI) \
	    $(LIBS) $(EXTRA_LINK)

xnet_MAGMA: $(BASE) jacobian_gpu.o $(DRIVER) $(MAGMA_OBJ) $(EXTRA_OBJ)
	$(LDR) $(LDRFLAGS) -o xnet_magma \
	    $(BASE) jacobian_gpu.o $(notdir $(MAGMA_OBJ)) $(DRIVER) \
	    $(MAGMA_LINK) $(EXTRA_LINK)

xnet_CDP: $(BASE) jacobian_gpu.o $(DRIVER) $(MAGMA_OBJ) $(CDP_LU_OBJ) $(EXTRA_OBJ)
	$(LDR) $(LDRFLAGS) -o xnet_cdp \
	    $(BASE) jacobian_gpu.o $(notdir $(MAGMA_OBJ)) $(notdir $(CDP_LU_OBJ)) cdpLink.o $(DRIVER) \
	    $(CDP_LU_LINK) $(MAGMA_LINK) $(EXTRA_LINK)

xnet_CULA: $(BASE) jacobian_gpu.o $(DRIVER) $(CULA_OBJ) $(EXTRA_OBJ)
	$(LDR) $(LDRFLAGS) -o xnet_cula \
	    $(BASE) jacobian_gpu.o $(notdir $(CULA_OBJ)) $(DRIVER) \
	    $(CULA_LINK) $(EXTRA_LINK)

xnet_LIBSCIACC: $(BASE) jacobian_gpu.o $(DRIVER) $(LIBSCIACC_OBJ) $(EXTRA_OBJ)
	$(LDR) $(LDRFLAGS) -o xnet_libsciacc \
	    $(BASE) jacobian_gpu.o $(notdir $(LIBSCIACC_OBJ)) $(DRIVER) \
	    $(LIBSCIACC_LINK) $(EXTRA_LINK)

net_setup: data.o conditions.o ffn.o net_preprocess.o net_setup.o $(NNU)
	$(LDR) $(LDRFLAGS) -o net_setup net_setup.o net_preprocess.o data.o conditions.o ffn.o $(NNU)

xinab: data.o init_abund.o 
	$(LDR)  $(LDRFLAGS) -o xinab \
	    init_abund.o data.o 

#
# Rules for compiling individual files.
#
MA28.o: MA28.f
	$(FC) $(FFLAGS) -c MA28.f -o $@
MA48.o: MA48.f
	$(FC) $(FFLAGS) -c MA48.f -o $@
$(MPI_FILES): %.o: %.f90
	$(FC_MPI) $(FFLAGS) -c $< -o $@
eosnom90.o: eosnom90.f
	$(FC) $(FFLAGS) -c eosnom90.f -o $@
eos_helm.o: eos_helm.f90
	$(FC) $(FFLAGS) -I$(HELMHOLTZ_PATH) -c eos_helm.f90 -o $@
press.o: press.f
	$(FC) $(FFLAGS) -c press.f -o $@
net.o: net.f90
	$(FC) $(FFLAGS) -c net.f90 -o $@
jacobian_$(MATRIX_SOLVER).o: jacobian_$(MATRIX_SOLVER).f90
	$(FC) $(FFLAGS) $(LAPACK_INC) -c jacobian_$(MATRIX_SOLVER).f90 -o $@
solve_be.o: solve_be.f90
	$(FC) $(FFLAGS) -c solve_be.f90 -o $@

ifeq ($(suffix LAPACK_SRC),.f90)
$(LAPACK_OBJ): %.o: %.f90
	$(FC) $(FFLAGS) $(LAPACK_INC) -c $< -o $(notdir $@)
else
$(LAPACK_OBJ): %.o: %.f
	$(FC) $(FFLAGS) $(LAPACK_INC) -c $< -o $(notdir $@)
endif

$(GPU_LAPACK_OBJ): %.o: %.f90
	$(FC) $(FFLAGS) $(GPU_LAPACK_INC) -c $< -o $(notdir $@)
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) $(CUDATK_INC)
%.o: %.f90
	$(FCL) $(FFLAGS) -c $< -o $@
%.o: %.f
	$(FCL) $(FFLAGS) -c $< -o $@
%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean:
	rm -f core *.o *.oo *.mod *.lst *.cub *.ptx *.i *.T *.diag xref.db
	rm -rf $(INLINE_DB)
tar:
	tar cvf net.tar ReadMe Changes Makefile *.f control th_const
#
#  Dependencies for object files
#
net.o:      full_net.o data.o flux.o match.o common.o control.o
flux.o:     full_net.o data.o match.o common.o
match.o:    data.o
full_net.o: data.o common.o
solve_be.o: jacobian_$(MATRIX_SOLVER).o
nu_nucleus.o: conditions.o
jacobian_$(MATRIX_SOLVER).o: $(EXTRA_OBJ) data.o full_net.o common.o
