
include ./Makefile.config

VER='\""0.1.0"\"'
COMPILE_TIME=$(shell date +"%d-%m-%Y %H:%M:%S")

LIBPATH_COMMON=
LINKERFLAG_COMMON=
INCPATH_COMMON=

ifeq ($(VENDOR), gnu)
	ifeq ($(LINEARALG), blas_lapack)
		LIBPATH_COMMON += "-L${BLAS_LIB}"
		LIBPATH_COMMON += "-L${LAPACK_LIB}"
		LINKERFLAG_COMMON += "-llapack -lblas"
	else ifeq ($(LINEARALG), openblas)
		LIBPATH_COMMON += "-L${OPENBLAS_LIB}"
		LINKERFLAG_COMMON += "-lopenblas"
	endif
endif

LIBPATH_COMMON += "-L${HDF5_LIB}"
LINKERFLAG_COMMON += "-lhdf5 -lhdf5_fortran -lm -ldl"

INCPATH_COMMON += "-I${HDF5_INC}"

ifeq ($(VENDOR), intel)
	MPIFC=${FC}
	LIBPATH_COMMON += "-L${MKL_LIB}"
endif

# ======================================= Compiler flag ======================================= #

COMP_FLAG_GNU_DEBUG=-Og -g -Wall -Wextra -Wunused-parameter -Wconversion $\
					-pedantic -fimplicit-none -fcheck=all -fbounds-check $\
					-fbacktrace -ftrapv -ffpe-trap=zero -finit-real=snan -cpp $\
					-D__COMPILE_TIME__='\""${COMPILE_TIME}"\"' -D__SPR_VERSION__=${VER}

COMP_FLAG_GNU_OPTIM=-O2 -flto -march=native -cpp -D__COMPILE_TIME__='\""${COMPILE_TIME}"\"' -D__SPR_VERSION__=${VER}

COMP_FLAG_INTEL_DEBUG=-O0 -debug -traceback -ftrapuv $\
					  -mp1 -fp-model precise $\
					  -W1 -warn all -warn noshape -warn uncalled -warn unused -warn usage $\
					  -check all -check bounds -check shape -check uninit -fpp $\
					  -D__COMPILE_TIME__='\""${COMPILE_TIME}"\"' -D__SPR_VERSION__=${VER}
#-diag-disable=10182

COMP_FLAG_INTEL_OPTIM=-O3 -xHost -fp-model source -ipo -flto -fpp -D__COMPILE_TIME__='\""${COMPILE_TIME}"\"' -D__SPR_VERSION__=${VER}

COMP_FLAG_COMMON=

ifeq ($(VENDOR), gnu)
	ifeq ($(DEBUG), 1)
		COMP_FLAG_COMMON=${COMP_FLAG_GNU_DEBUG}
	else 
		COMP_FLAG_COMMON=${COMP_FLAG_GNU_OPTIM}
	endif
endif

ifeq ($(VENDOR), intel)
	INTEL = yes
	ifeq ($(DEBUG), 1)
		COMP_FLAG_COMMON=${COMP_FLAG_INTEL_DEBUG}
	else 
		COMP_FLAG_COMMON=${COMP_FLAG_INTEL_OPTIM}
	endif
endif

# ======================================= Compiler flag ======================================= #

# ========================================== Python3 ========================================== #

ifeq ($(PY3_PATH),)
	PY3_PATH := $(shell which python3)
endif

ifeq ($(PY3_PATH),)
	PY3_PATH := $(shell which python)
endif

# ========================================== Python3 ========================================== #

PARENT_DIR=${PWD}

ifeq ($(INSTALL_PATH),)
    INSTALL_PATH := ${PARENT_DIR}/bin
endif

.PHONY: all omp FC2 FC3 FC4 CreateDM2 ReconstFC2 CreateDM3 ReconstFC3 CreateDM4 ReconstFC4 Renorm Disp Long_DD LSq TSS FC3omp FC4omp CreateDM2omp ReconstFC2omp CreateDM3omp ReconstFC3omp CreateDM4omp ReconstFC4omp FreeEnergy LineWidth kappaIter kappa clean

all: FC2 FC3 FC4 CreateDM2 ReconstFC2 CreateDM3 ReconstFC3 CreateDM4 ReconstFC4 Renorm Disp Long_DD LSq TSS FreeEnergy LineWidth kappaIter kappa

omp: FC3omp FC4omp CreateDM2omp ReconstFC2omp CreateDM3omp ReconstFC3omp CreateDM4omp ReconstFC4omp

COMP_FLAG=
LINKERFLAG=
LIBPATH=
EXEC_NAME=
INCPATH=

FC2DIR=${PARENT_DIR}/src/FC2/Src
FC2: 
	@if [ -f ${FC2DIR}/Makefile ]; then \
		rm ${FC2DIR}/Makefile; \
	fi
	@if [ -f ${FC2DIR}/variable.tmp ]; then \
		rm ${FC2DIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "FC2.x")
	@rm -f ${FC2DIR}/*.o ${FC2DIR}/*.mod ${FC2DIR}/${EXEC_NAME}
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qmkl=parallel")
	$(eval LINKERFLAG += "-liomp5")
endif
	$(eval LIBPATH += "-L${SPGLIB_LIB}")
	$(eval LINKERFLAG += "-lsymspg")
	@echo 
	@touch ${FC2DIR}/Makefile
	@touch ${FC2DIR}/variable.tmp
	@echo "FC=$(FC)" >> ${FC2DIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${FC2DIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${FC2DIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${FC2DIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${FC2DIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${FC2DIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${FC2DIR}/variable.tmp
	@cat ${FC2DIR}/variable.tmp >> ${FC2DIR}/Makefile
	@cat ${FC2DIR}/Makefile.template >> ${FC2DIR}/Makefile
	$(MAKE) -C ${FC2DIR}
	@install -d ${INSTALL_PATH}
	@install ${FC2DIR}/${EXEC_NAME} ${INSTALL_PATH}


FC3DIR=${PARENT_DIR}/src/FC3/Src
FC3: 
	@if [ -f ${FC3DIR}/Makefile ]; then \
		rm ${FC3DIR}/Makefile; \
	fi
	@if [ -f ${FC3DIR}/variable.tmp ]; then \
		rm ${FC3DIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "FC3.x")
	@rm -f ${FC3DIR}/*.o ${FC3DIR}/*.mod ${FC3DIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qmkl=parallel")
	$(eval LINKERFLAG += "-liomp5")
endif
	$(eval LIBPATH += "-L${SPGLIB_LIB}")
	$(eval LINKERFLAG += "-lsymspg")
	@echo 
	@touch ${FC3DIR}/Makefile
	@touch ${FC3DIR}/variable.tmp
	@echo "FC=$(FC)" >> ${FC3DIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${FC3DIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${FC3DIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${FC3DIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${FC3DIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${FC3DIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${FC3DIR}/variable.tmp
	@cat ${FC3DIR}/variable.tmp >> ${FC3DIR}/Makefile
	@cat ${FC3DIR}/Makefile.template >> ${FC3DIR}/Makefile
	$(MAKE) -C ${FC3DIR}
	@install -d ${INSTALL_PATH}
	@install ${FC3DIR}/${EXEC_NAME} ${INSTALL_PATH}


FC4DIR=${PARENT_DIR}/src/FC4/Src
FC4: 
	@if [ -f ${FC4DIR}/Makefile ]; then \
		rm ${FC4DIR}/Makefile; \
	fi
	@if [ -f ${FC4DIR}/variable.tmp ]; then \
		rm ${FC4DIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "FC4.x")
	@rm -f ${FC4DIR}/*.o ${FC4DIR}/*.mod ${FC4DIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qmkl=parallel")
	$(eval LINKERFLAG += "-liomp5")
endif
	$(eval LIBPATH += "-L${SPGLIB_LIB}")
	$(eval LINKERFLAG += "-lsymspg")
	@echo 
	@touch ${FC4DIR}/Makefile
	@touch ${FC4DIR}/variable.tmp
	@echo "FC=$(FC)" >> ${FC4DIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${FC4DIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${FC4DIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${FC4DIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${FC4DIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${FC4DIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${FC4DIR}/variable.tmp
	@cat ${FC4DIR}/variable.tmp >> ${FC4DIR}/Makefile
	@cat ${FC4DIR}/Makefile.template >> ${FC4DIR}/Makefile
	$(MAKE) -C ${FC4DIR}
	@install -d ${INSTALL_PATH}
	@install ${FC4DIR}/${EXEC_NAME} ${INSTALL_PATH}


CreateDM2DIR=${PARENT_DIR}/src/CreateDM2/Src_Dyn
CreateDM2: 
	@if [ -f ${CreateDM2DIR}/Makefile ]; then \
		rm ${CreateDM2DIR}/Makefile; \
	fi
	@if [ -f ${CreateDM2DIR}/variable.tmp ]; then \
		rm ${CreateDM2DIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "create_disp_mat2.x")
	@rm -f ${CreateDM2DIR}/*.o ${CreateDM2DIR}/*.mod ${CreateDM2DIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
ifdef CAF_SINGLE
	$(eval COMP_FLAG += "-coarray=single -qmkl=sequential")
else
ifdef CO_CONF_FILE_CREATEDM2
	$(eval COMP_FLAG += "-coarray=distributed -coarray-config-file=${CO_CONF_FILE_CREATEDM2} -qmkl=sequential")
else
	$(eval COMP_FLAG += "-coarray=shared -qmkl=sequential")
endif
endif
else
	$(eval COMP_FLAG += "-DGNU -fcoarray=lib")
	$(eval LIBPATH += "-L${COARR_LIB}")
	$(eval LINKERFLAG += "-lcaf_mpi")
endif
	@echo 
	@touch ${CreateDM2DIR}/Makefile
	@touch ${CreateDM2DIR}/variable.tmp
	@echo "FC=$(FC)" >> ${CreateDM2DIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${CreateDM2DIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${CreateDM2DIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${CreateDM2DIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${CreateDM2DIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${CreateDM2DIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${CreateDM2DIR}/variable.tmp
	@cat ${CreateDM2DIR}/variable.tmp >> ${CreateDM2DIR}/Makefile
	@cat ${CreateDM2DIR}/Makefile.template >> ${CreateDM2DIR}/Makefile
	$(MAKE) -C ${CreateDM2DIR}
	@install -d ${INSTALL_PATH}
	@install ${CreateDM2DIR}/${EXEC_NAME} ${INSTALL_PATH}


ReconstFC2DIR=${PARENT_DIR}/src/ReconstFC2/Src
ReconstFC2: 
	@if [ -f ${ReconstFC2DIR}/Makefile ]; then \
		rm ${ReconstFC2DIR}/Makefile; \
	fi
	@if [ -f ${ReconstFC2DIR}/variable.tmp ]; then \
		rm ${ReconstFC2DIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "reconst2.x")
	@rm -f ${ReconstFC2DIR}/*.o ${ReconstFC2DIR}/*.mod ${ReconstFC2DIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
ifdef CAF_SINGLE
	$(eval COMP_FLAG += "-coarray=single -qmkl=sequential")
else
ifdef CO_CONF_FILE_RECONST2
	$(eval COMP_FLAG += "-coarray=distributed -coarray-config-file=${CO_CONF_FILE_RECONST2} -qmkl=sequential")
else
	$(eval COMP_FLAG += "-coarray=shared -qmkl=sequential")
endif
endif
else
	$(eval COMP_FLAG += "-fcoarray=lib")
	$(eval LIBPATH += "-L${COARR_LIB}")
	$(eval LINKERFLAG += "-lcaf_mpi")
endif
	@echo 
	@touch ${ReconstFC2DIR}/Makefile
	@touch ${ReconstFC2DIR}/variable.tmp
	@echo "FC=$(FC)" >> ${ReconstFC2DIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${ReconstFC2DIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${ReconstFC2DIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${ReconstFC2DIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${ReconstFC2DIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${ReconstFC2DIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${ReconstFC2DIR}/variable.tmp
	@cat ${ReconstFC2DIR}/variable.tmp >> ${ReconstFC2DIR}/Makefile
	@cat ${ReconstFC2DIR}/Makefile.template >> ${ReconstFC2DIR}/Makefile
	$(MAKE) -C ${ReconstFC2DIR}
	@install -d ${INSTALL_PATH}
	@install ${ReconstFC2DIR}/${EXEC_NAME} ${INSTALL_PATH}


CreateDM3DIR=${PARENT_DIR}/src/CreateDM3/Src_Dyn
CreateDM3: 
	@if [ -f ${CreateDM3DIR}/Makefile ]; then \
		rm ${CreateDM3DIR}/Makefile; \
	fi
	@if [ -f ${CreateDM3DIR}/variable.tmp ]; then \
		rm ${CreateDM3DIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "create_disp_mat3.x")
	@rm -f ${CreateDM3DIR}/*.o ${CreateDM3DIR}/*.mod ${CreateDM3DIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
ifdef CAF_SINGLE
	$(eval COMP_FLAG += "-coarray=single -qmkl=sequential")
else
ifdef CO_CONF_FILE_CREATEDM3
	$(eval COMP_FLAG += "-coarray=distributed -coarray-config-file=${CO_CONF_FILE_CREATEDM3} -qmkl=sequential")
else
	$(eval COMP_FLAG += "-coarray=shared -qmkl=sequential")
endif
endif
else
	$(eval COMP_FLAG += "-DGNU -fcoarray=lib")
	$(eval LIBPATH += "-L${COARR_LIB}")
	$(eval LINKERFLAG += "-lcaf_mpi")
endif
	@echo 
	@touch ${CreateDM3DIR}/Makefile
	@touch ${CreateDM3DIR}/variable.tmp
	@echo "FC=$(FC)" >> ${CreateDM3DIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${CreateDM3DIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${CreateDM3DIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${CreateDM3DIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${CreateDM3DIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${CreateDM3DIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${CreateDM3DIR}/variable.tmp
	@cat ${CreateDM3DIR}/variable.tmp >> ${CreateDM3DIR}/Makefile
	@cat ${CreateDM3DIR}/Makefile.template >> ${CreateDM3DIR}/Makefile
	$(MAKE) -C ${CreateDM3DIR}
	@install -d ${INSTALL_PATH}
	@install ${CreateDM3DIR}/${EXEC_NAME} ${INSTALL_PATH}


ReconstFC3DIR=${PARENT_DIR}/src/ReconstFC3/Src
ReconstFC3: 
	@if [ -f ${ReconstFC3DIR}/Makefile ]; then \
		rm ${ReconstFC3DIR}/Makefile; \
	fi
	@if [ -f ${ReconstFC3DIR}/variable.tmp ]; then \
		rm ${ReconstFC3DIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "reconst3.x")
	@rm -f ${ReconstFC3DIR}/*.o ${ReconstFC3DIR}/*.mod ${ReconstFC3DIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
ifdef CAF_SINGLE
	$(eval COMP_FLAG += "-coarray=single -qmkl=sequential")
else
ifdef CO_CONF_FILE_RECONST3
	$(eval COMP_FLAG += "-coarray=distributed -coarray-config-file=${CO_CONF_FILE_RECONST3} -qmkl=sequential")
else
	$(eval COMP_FLAG += "-coarray=shared -qmkl=sequential")
endif
endif
else
	$(eval COMP_FLAG += "-fcoarray=lib")
	$(eval LIBPATH += "-L${COARR_LIB}")
	$(eval LINKERFLAG += "-lcaf_mpi")
endif
	@echo 
	@touch ${ReconstFC3DIR}/Makefile
	@touch ${ReconstFC3DIR}/variable.tmp
	@echo "FC=$(FC)" >> ${ReconstFC3DIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${ReconstFC3DIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${ReconstFC3DIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${ReconstFC3DIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${ReconstFC3DIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${ReconstFC3DIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${ReconstFC3DIR}/variable.tmp
	@cat ${ReconstFC3DIR}/variable.tmp >> ${ReconstFC3DIR}/Makefile
	@cat ${ReconstFC3DIR}/Makefile.template >> ${ReconstFC3DIR}/Makefile
	$(MAKE) -C ${ReconstFC3DIR}
	@install -d ${INSTALL_PATH}
	@install ${ReconstFC3DIR}/${EXEC_NAME} ${INSTALL_PATH}


CreateDM4DIR=${PARENT_DIR}/src/CreateDM4/Src_Dyn
CreateDM4: 
	@if [ -f ${CreateDM4DIR}/Makefile ]; then \
		rm ${CreateDM4DIR}/Makefile; \
	fi
	@if [ -f ${CreateDM4DIR}/variable.tmp ]; then \
		rm ${CreateDM4DIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "create_disp_mat4.x")
	@rm -f ${CreateDM4DIR}/*.o ${CreateDM4DIR}/*.mod ${CreateDM4DIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
ifdef CAF_SINGLE
	$(eval COMP_FLAG += "-coarray=single -qmkl=sequential")
else
ifdef CO_CONF_FILE_CREATEDM4
	$(eval COMP_FLAG += "-coarray=distributed -coarray-config-file=${CO_CONF_FILE_CREATEDM4} -qmkl=sequential")
else
	$(eval COMP_FLAG += "-coarray=shared -qmkl=sequential")
endif
endif
else
	$(eval COMP_FLAG += "-DGNU -fcoarray=lib")
	$(eval LIBPATH += "-L${COARR_LIB}")
	$(eval LINKERFLAG += "-lcaf_mpi")
endif
	@echo 
	@touch ${CreateDM4DIR}/Makefile
	@touch ${CreateDM4DIR}/variable.tmp
	@echo "FC=$(FC)" >> ${CreateDM4DIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${CreateDM4DIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${CreateDM4DIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${CreateDM4DIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${CreateDM4DIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${CreateDM4DIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${CreateDM4DIR}/variable.tmp
	@cat ${CreateDM4DIR}/variable.tmp >> ${CreateDM4DIR}/Makefile
	@cat ${CreateDM4DIR}/Makefile.template >> ${CreateDM4DIR}/Makefile
	$(MAKE) -C ${CreateDM4DIR}
	@install -d ${INSTALL_PATH}
	@install ${CreateDM4DIR}/${EXEC_NAME} ${INSTALL_PATH}


ReconstFC4DIR=${PARENT_DIR}/src/ReconstFC4/Src
ReconstFC4: 
	@if [ -f ${ReconstFC4DIR}/Makefile ]; then \
		rm ${ReconstFC4DIR}/Makefile; \
	fi
	@if [ -f ${ReconstFC4DIR}/variable.tmp ]; then \
		rm ${ReconstFC4DIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "reconst4.x")
	@rm -f ${ReconstFC4DIR}/*.o ${ReconstFC4DIR}/*.mod ${ReconstFC4DIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
ifdef CAF_SINGLE
	$(eval COMP_FLAG += "-coarray=single -qmkl=sequential")
else
ifdef CO_CONF_FILE_RECONST4
	$(eval COMP_FLAG += "-coarray=distributed -coarray-config-file=${CO_CONF_FILE_RECONST4} -qmkl=sequential")
else
	$(eval COMP_FLAG += "-coarray=shared -qmkl=sequential")
endif
endif
else
	$(eval COMP_FLAG += "-fcoarray=lib")
	$(eval LIBPATH += "-L${COARR_LIB}")
	$(eval LINKERFLAG += "-lcaf_mpi")
endif
	@echo 
	@touch ${ReconstFC4DIR}/Makefile
	@touch ${ReconstFC4DIR}/variable.tmp
	@echo "FC=$(FC)" >> ${ReconstFC4DIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${ReconstFC4DIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${ReconstFC4DIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${ReconstFC4DIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${ReconstFC4DIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${ReconstFC4DIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${ReconstFC4DIR}/variable.tmp
	@cat ${ReconstFC4DIR}/variable.tmp >> ${ReconstFC4DIR}/Makefile
	@cat ${ReconstFC4DIR}/Makefile.template >> ${ReconstFC4DIR}/Makefile
	$(MAKE) -C ${ReconstFC4DIR}
	@install -d ${INSTALL_PATH}
	@install ${ReconstFC4DIR}/${EXEC_NAME} ${INSTALL_PATH}


RenormDIR=${PARENT_DIR}/src/Renorm/Src
Renorm: 
	@if [ -f ${RenormDIR}/Makefile ]; then \
		rm ${RenormDIR}/Makefile; \
	fi
	@if [ -f ${RenormDIR}/variable.tmp ]; then \
		rm ${RenormDIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "renorm.x")
	@rm -f ${RenormDIR}/*.o ${RenormDIR}/*.mod ${RenormDIR}/${EXEC_NAME}
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qmkl=sequential")
endif
	@echo 
	@touch ${RenormDIR}/Makefile
	@touch ${RenormDIR}/variable.tmp
	@echo "FC=$(FC)" >> ${RenormDIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${RenormDIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${RenormDIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${RenormDIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${RenormDIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${RenormDIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${RenormDIR}/variable.tmp
	@cat ${RenormDIR}/variable.tmp >> ${RenormDIR}/Makefile
	@cat ${RenormDIR}/Makefile.template >> ${RenormDIR}/Makefile
	$(MAKE) -C ${RenormDIR}
	@install -d ${INSTALL_PATH}
	@install ${RenormDIR}/${EXEC_NAME} ${INSTALL_PATH}


DispDIR=${PARENT_DIR}/src/Disp/Src
Disp: 
	@if [ -f ${DispDIR}/Makefile ]; then \
		rm ${DispDIR}/Makefile; \
	fi
	@if [ -f ${DispDIR}/variable.tmp ]; then \
		rm ${DispDIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "disp.x")
	@rm -f ${DispDIR}/*.o ${DispDIR}/*.mod ${DispDIR}/${EXEC_NAME}
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qmkl=sequential")
endif
	@echo 
	@touch ${DispDIR}/Makefile
	@touch ${DispDIR}/variable.tmp
	@echo "FC=$(FC)" >> ${DispDIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${DispDIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${DispDIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${DispDIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${DispDIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${DispDIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${DispDIR}/variable.tmp
	@cat ${DispDIR}/variable.tmp >> ${DispDIR}/Makefile
	@cat ${DispDIR}/Makefile.template >> ${DispDIR}/Makefile
	$(MAKE) -C ${DispDIR}
	@install -d ${INSTALL_PATH}
	@install ${DispDIR}/${EXEC_NAME} ${INSTALL_PATH}


Long_DDDIR=${PARENT_DIR}/src/Long_DD/Src
Long_DD: 
	@if [ -f ${Long_DDDIR}/Makefile ]; then \
		rm ${Long_DDDIR}/Makefile; \
	fi
	@if [ -f ${Long_DDDIR}/variable.tmp ]; then \
		rm ${Long_DDDIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "longEw.x")
	@rm -f ${Long_DDDIR}/*.o ${Long_DDDIR}/*.mod ${Long_DDDIR}/${EXEC_NAME}
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qmkl=sequential")
endif
	@echo 
	@touch ${Long_DDDIR}/Makefile
	@touch ${Long_DDDIR}/variable.tmp
	@echo "FC=$(FC)" >> ${Long_DDDIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${Long_DDDIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${Long_DDDIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${Long_DDDIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${Long_DDDIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${Long_DDDIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${Long_DDDIR}/variable.tmp
	@cat ${Long_DDDIR}/variable.tmp >> ${Long_DDDIR}/Makefile
	@cat ${Long_DDDIR}/Makefile.template >> ${Long_DDDIR}/Makefile
	$(MAKE) -C ${Long_DDDIR}
	@install -d ${INSTALL_PATH}
	@install ${Long_DDDIR}/${EXEC_NAME} ${INSTALL_PATH}


LsqDIR=${PARENT_DIR}/src/LSq/Src
LSq: 
	@if [ -f ${LsqDIR}/Makefile ]; then \
		rm ${LsqDIR}/Makefile; \
	fi
	@if [ -f ${LsqDIR}/variable.tmp ]; then \
		rm ${LsqDIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "LSq.x")
	@rm -f ${LsqDIR}/*.o ${LsqDIR}/*.mod ${LsqDIR}/${EXEC_NAME}
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qmkl=parallel")
	$(eval LINKERFLAG += "-liomp5")
endif
	@echo 
	@touch ${LsqDIR}/Makefile
	@touch ${LsqDIR}/variable.tmp
	@echo "FC=$(FC)" >> ${LsqDIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${LsqDIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${LsqDIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${LsqDIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${LsqDIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${LsqDIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${LsqDIR}/variable.tmp
	@cat ${LsqDIR}/variable.tmp >> ${LsqDIR}/Makefile
	@cat ${LsqDIR}/Makefile.template >> ${LsqDIR}/Makefile
	$(MAKE) -C ${LsqDIR}
	@install -d ${INSTALL_PATH}
	@install ${LsqDIR}/${EXEC_NAME} ${INSTALL_PATH}


TSSDIR=${PARENT_DIR}/src/TSS/Src
TSS: 
	@if [ -f ${TSSDIR}/Makefile ]; then \
		rm ${TSSDIR}/Makefile; \
	fi
	@if [ -f ${TSSDIR}/variable.tmp ]; then \
		rm ${TSSDIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "TSS.x")
	@rm -f ${TSSDIR}/*.o ${TSSDIR}/*.mod ${TSSDIR}/${EXEC_NAME}
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qmkl=sequential")
else
	$(eval COMP_FLAG += "-DGNU")
endif
	@echo 
	@touch ${TSSDIR}/Makefile
	@touch ${TSSDIR}/variable.tmp
	@echo "FC=$(FC)" >> ${TSSDIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${TSSDIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${TSSDIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${TSSDIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${TSSDIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${TSSDIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${TSSDIR}/variable.tmp
	@cat ${TSSDIR}/variable.tmp >> ${TSSDIR}/Makefile
	@cat ${TSSDIR}/Makefile.template >> ${TSSDIR}/Makefile
	$(MAKE) -C ${TSSDIR}
	@install -d ${INSTALL_PATH}
	@install ${TSSDIR}/${EXEC_NAME} ${INSTALL_PATH}


FreeEnergyDIR=${PARENT_DIR}/src/FreeEnergy/Src
FreeEnergy: 
	@if [ -f ${FreeEnergyDIR}/Makefile ]; then \
		rm ${FreeEnergyDIR}/Makefile; \
	fi
	@if [ -f ${FreeEnergyDIR}/variable.tmp ]; then \
		rm ${FreeEnergyDIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "FreeEng.x")
	@rm -f ${FreeEnergyDIR}/*.o ${FreeEnergyDIR}/*.mod ${FreeEnergyDIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
ifdef CAF_SINGLE
	$(eval COMP_FLAG += "-coarray=single -qmkl=sequential")
else
ifdef CO_CONF_FILE_FREE_ENG
	$(eval COMP_FLAG += "-coarray=distributed -coarray-config-file=${CO_CONF_FILE_FREE_ENG} -qmkl=sequential")
else
	$(eval COMP_FLAG += "-coarray=shared -qmkl=sequential")
endif
endif
else
	$(eval COMP_FLAG += "-fcoarray=lib")
	$(eval LIBPATH += "-L${COARR_LIB}")
	$(eval LINKERFLAG += "-lcaf_mpi")
endif
	$(eval LIBPATH += "-L${SPGLIB_LIB}")
	$(eval LINKERFLAG += "-lsymspg")
	@echo 
	@touch ${FreeEnergyDIR}/Makefile
	@touch ${FreeEnergyDIR}/variable.tmp
	@echo "FC=$(FC)" >> ${FreeEnergyDIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${FreeEnergyDIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${FreeEnergyDIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${FreeEnergyDIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${FreeEnergyDIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${FreeEnergyDIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${FreeEnergyDIR}/variable.tmp
	@cat ${FreeEnergyDIR}/variable.tmp >> ${FreeEnergyDIR}/Makefile
	@cat ${FreeEnergyDIR}/Makefile.template >> ${FreeEnergyDIR}/Makefile
	$(MAKE) -C ${FreeEnergyDIR}
	@install -d ${INSTALL_PATH}
	@install ${FreeEnergyDIR}/${EXEC_NAME} ${INSTALL_PATH}


LineWidthDIR=${PARENT_DIR}/src/LineWidth/Src
LineWidth: 
	@if [ -f ${LineWidthDIR}/Makefile ]; then \
		rm ${LineWidthDIR}/Makefile; \
	fi
	@if [ -f ${LineWidthDIR}/variable.tmp ]; then \
		rm ${LineWidthDIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "LineWidth.x")
	@rm -f ${LineWidthDIR}/*.o ${LineWidthDIR}/*.mod ${LineWidthDIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
ifdef CAF_SINGLE
	$(eval COMP_FLAG += "-coarray=single -qmkl=sequential")
else
ifdef CO_CONF_FILE_LINEWIDTH
	$(eval COMP_FLAG += "-coarray=distributed -coarray-config-file=${CO_CONF_FILE_LINEWIDTH} -qmkl=sequential")
else
	$(eval COMP_FLAG += "-coarray=shared -qmkl=sequential")
endif
endif
else
	$(eval COMP_FLAG += "-fcoarray=lib")
	$(eval LIBPATH += "-L${COARR_LIB}")
	$(eval LINKERFLAG += "-lcaf_mpi")
endif
	@echo 
	@touch ${LineWidthDIR}/Makefile
	@touch ${LineWidthDIR}/variable.tmp
	@echo "FC=$(FC)" >> ${LineWidthDIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${LineWidthDIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${LineWidthDIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${LineWidthDIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${LineWidthDIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${LineWidthDIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${LineWidthDIR}/variable.tmp
	@cat ${LineWidthDIR}/variable.tmp >> ${LineWidthDIR}/Makefile
	@cat ${LineWidthDIR}/Makefile.template >> ${LineWidthDIR}/Makefile
	$(MAKE) -C ${LineWidthDIR}
	@install -d ${INSTALL_PATH}
	@install ${LineWidthDIR}/${EXEC_NAME} ${INSTALL_PATH}


kappaIterDIR=${PARENT_DIR}/src/Kappa/Src_iterative
kappaIter: 
	@if [ -f ${kappaIterDIR}/Makefile ]; then \
		rm ${kappaIterDIR}/Makefile; \
	fi
	@if [ -f ${kappaIterDIR}/variable.tmp ]; then \
		rm ${kappaIterDIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "kappa_iter.x")
	@rm -f ${kappaIterDIR}/*.o ${kappaIterDIR}/*.mod ${kappaIterDIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
ifdef CAF_SINGLE
	$(eval COMP_FLAG += "-coarray=single -qmkl=sequential")
else
ifdef CO_CONF_FILE_KAPPA_ITER
	$(eval COMP_FLAG += "-coarray=distributed -coarray-config-file=${CO_CONF_FILE_KAPPA_ITER} -qmkl=sequential")
else
	$(eval COMP_FLAG += "-coarray=shared -qmkl=sequential")
endif
endif
else
	$(eval COMP_FLAG += "-fcoarray=lib -DGNU")
	$(eval LIBPATH += "-L${COARR_LIB}")
	$(eval LINKERFLAG += "-lcaf_mpi")
endif
	@echo 
	@touch ${kappaIterDIR}/Makefile
	@touch ${kappaIterDIR}/variable.tmp
	@echo "FC=$(FC)" >> ${kappaIterDIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${kappaIterDIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${kappaIterDIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${kappaIterDIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${kappaIterDIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${kappaIterDIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${kappaIterDIR}/variable.tmp
	@cat ${kappaIterDIR}/variable.tmp >> ${kappaIterDIR}/Makefile
	@cat ${kappaIterDIR}/Makefile.template >> ${kappaIterDIR}/Makefile
	$(MAKE) -C ${kappaIterDIR}
	@install -d ${INSTALL_PATH}
	@install ${kappaIterDIR}/${EXEC_NAME} ${INSTALL_PATH}


kappaDIR=${PARENT_DIR}/src/Kappa/Src
kappa: 
	@if [ -f ${kappaDIR}/Makefile ]; then \
		rm ${kappaDIR}/Makefile; \
	fi
	@if [ -f ${kappaDIR}/variable.tmp ]; then \
		rm ${kappaDIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "kappa.x")
	@rm -f ${kappaDIR}/*.o ${kappaDIR}/*.mod ${kappaDIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
ifdef CAF_SINGLE
	$(eval COMP_FLAG += "-coarray=single -qmkl=sequential")
else
ifdef CO_CONF_FILE_KAPPA
	$(eval COMP_FLAG += "-coarray=distributed -coarray-config-file=${CO_CONF_FILE_KAPPA} -qmkl=sequential")
else
	$(eval COMP_FLAG += "-coarray=shared -qmkl=sequential")
endif
endif
else
	$(eval COMP_FLAG += "-fcoarray=lib -DGNU")
	$(eval LIBPATH += "-L${COARR_LIB}")
	$(eval LINKERFLAG += "-lcaf_mpi")
endif
	$(eval LIBPATH += "-L${SPGLIB_LIB}")
	$(eval LINKERFLAG += "-lsymspg")
	@echo 
	@touch ${kappaDIR}/Makefile
	@touch ${kappaDIR}/variable.tmp
	@echo "FC=$(FC)" >> ${kappaDIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${kappaDIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${kappaDIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${kappaDIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${kappaDIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${kappaDIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${kappaDIR}/variable.tmp
	@cat ${kappaDIR}/variable.tmp >> ${kappaDIR}/Makefile
	@cat ${kappaDIR}/Makefile.template >> ${kappaDIR}/Makefile
	$(MAKE) -C ${kappaDIR}
	@install -d ${INSTALL_PATH}
	@install ${kappaDIR}/${EXEC_NAME} ${INSTALL_PATH}


FC3DIR=${PARENT_DIR}/src/FC3/Src
FC3omp: 
	@if [ -f ${FC3DIR}/Makefile ]; then \
		rm ${FC3DIR}/Makefile; \
	fi
	@if [ -f ${FC3DIR}/variable.tmp ]; then \
		rm ${FC3DIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "FC3_omp.x")
	@rm -f ${FC3DIR}/*.o ${FC3DIR}/*.mod ${FC3DIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qmkl=parallel -qopenmp")
	$(eval COMP_FLAG += "-DOMP")
	$(eval LINKERFLAG += "-liomp5")
else
	$(eval COMP_FLAG += "-DOMP -fopenmp")
endif
	$(eval LIBPATH += "-L${SPGLIB_LIB}")
	$(eval LINKERFLAG += "-lsymspg")
	@echo 
	@touch ${FC3DIR}/Makefile
	@touch ${FC3DIR}/variable.tmp
	@echo "FC=$(FC)" >> ${FC3DIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${FC3DIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${FC3DIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${FC3DIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${FC3DIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${FC3DIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${FC3DIR}/variable.tmp
	@cat ${FC3DIR}/variable.tmp >> ${FC3DIR}/Makefile
	@cat ${FC3DIR}/Makefile.template >> ${FC3DIR}/Makefile
	$(MAKE) -C ${FC3DIR}
	@install -d ${INSTALL_PATH}
	@install ${FC3DIR}/${EXEC_NAME} ${INSTALL_PATH}


FC4DIR=${PARENT_DIR}/src/FC4/Src
FC4omp: 
	@if [ -f ${FC4DIR}/Makefile ]; then \
		rm ${FC4DIR}/Makefile; \
	fi
	@if [ -f ${FC4DIR}/variable.tmp ]; then \
		rm ${FC4DIR}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "FC4_omp.x")
	@rm -f ${FC4DIR}/*.o ${FC4DIR}/*.mod ${FC4DIR}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qmkl=parallel -qopenmp")
	$(eval COMP_FLAG += "-DOMP")
	$(eval LINKERFLAG += "-liomp5")
else
	$(eval COMP_FLAG += "-DOMP -fopenmp")
endif
	$(eval LIBPATH += "-L${SPGLIB_LIB}")
	$(eval LINKERFLAG += "-lsymspg")
	@echo 
	@touch ${FC4DIR}/Makefile
	@touch ${FC4DIR}/variable.tmp
	@echo "FC=$(FC)" >> ${FC4DIR}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${FC4DIR}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${FC4DIR}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${FC4DIR}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${FC4DIR}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${FC4DIR}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${FC4DIR}/variable.tmp
	@cat ${FC4DIR}/variable.tmp >> ${FC4DIR}/Makefile
	@cat ${FC4DIR}/Makefile.template >> ${FC4DIR}/Makefile
	$(MAKE) -C ${FC4DIR}
	@install -d ${INSTALL_PATH}
	@install ${FC4DIR}/${EXEC_NAME} ${INSTALL_PATH}


CreateDM2DIR_OMP=${PARENT_DIR}/src/CreateDM2/Src_OMP
CreateDM2omp: 
	@if [ -f ${CreateDM2DIR_OMP}/Makefile ]; then \
		rm ${CreateDM2DIR_OMP}/Makefile; \
	fi
	@if [ -f ${CreateDM2DIR_OMP}/variable.tmp ]; then \
		rm ${CreateDM2DIR_OMP}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "create_disp_mat2_omp.x")
	@rm -f ${CreateDM2DIR_OMP}/*.o ${CreateDM2DIR_OMP}/*.mod ${CreateDM2DIR_OMP}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qopenmp -qmkl=sequential")
else
	$(eval COMP_FLAG += "-fopenmp")
endif
	@echo 
	@touch ${CreateDM2DIR_OMP}/Makefile
	@touch ${CreateDM2DIR_OMP}/variable.tmp
	@echo "FC=$(FC)" >> ${CreateDM2DIR_OMP}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${CreateDM2DIR_OMP}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${CreateDM2DIR_OMP}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${CreateDM2DIR_OMP}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${CreateDM2DIR_OMP}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${CreateDM2DIR_OMP}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${CreateDM2DIR_OMP}/variable.tmp
	@cat ${CreateDM2DIR_OMP}/variable.tmp >> ${CreateDM2DIR_OMP}/Makefile
	@cat ${CreateDM2DIR_OMP}/Makefile.template >> ${CreateDM2DIR_OMP}/Makefile
	$(MAKE) -C ${CreateDM2DIR_OMP}
	@install -d ${INSTALL_PATH}
	@install ${CreateDM2DIR_OMP}/${EXEC_NAME} ${INSTALL_PATH}


ReconstFC2DIR_OMP=${PARENT_DIR}/src/ReconstFC2/Src_OMP
ReconstFC2omp: 
	@if [ -f ${ReconstFC2DIR_OMP}/Makefile ]; then \
		rm ${ReconstFC2DIR_OMP}/Makefile; \
	fi
	@if [ -f ${ReconstFC2DIR_OMP}/variable.tmp ]; then \
		rm ${ReconstFC2DIR_OMP}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "reconst2_omp.x")
	@rm -f ${ReconstFC2DIR_OMP}/*.o ${ReconstFC2DIR_OMP}/*.mod ${ReconstFC2DIR_OMP}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qopenmp -qmkl=sequential")
else
	$(eval COMP_FLAG += "-fopenmp")
endif
	@echo 
	@touch ${ReconstFC2DIR_OMP}/Makefile
	@touch ${ReconstFC2DIR_OMP}/variable.tmp
	@echo "FC=$(FC)" >> ${ReconstFC2DIR_OMP}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${ReconstFC2DIR_OMP}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${ReconstFC2DIR_OMP}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${ReconstFC2DIR_OMP}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${ReconstFC2DIR_OMP}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${ReconstFC2DIR_OMP}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${ReconstFC2DIR_OMP}/variable.tmp
	@cat ${ReconstFC2DIR_OMP}/variable.tmp >> ${ReconstFC2DIR_OMP}/Makefile
	@cat ${ReconstFC2DIR_OMP}/Makefile.template >> ${ReconstFC2DIR_OMP}/Makefile
	$(MAKE) -C ${ReconstFC2DIR_OMP}
	@install -d ${INSTALL_PATH}
	@install ${ReconstFC2DIR_OMP}/${EXEC_NAME} ${INSTALL_PATH}


CreateDM3DIR_OMP=${PARENT_DIR}/src/CreateDM3/Src_OMP
CreateDM3omp: 
	@if [ -f ${CreateDM3DIR_OMP}/Makefile ]; then \
		rm ${CreateDM3DIR_OMP}/Makefile; \
	fi
	@if [ -f ${CreateDM3DIR_OMP}/variable.tmp ]; then \
		rm ${CreateDM3DIR_OMP}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "create_disp_mat3_omp.x")
	@rm -f ${CreateDM3DIR_OMP}/*.o ${CreateDM3DIR_OMP}/*.mod ${CreateDM3DIR_OMP}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qopenmp -qmkl=sequential")
else
	$(eval COMP_FLAG += "-fopenmp")
endif
	@echo 
	@touch ${CreateDM3DIR_OMP}/Makefile
	@touch ${CreateDM3DIR_OMP}/variable.tmp
	@echo "FC=$(FC)" >> ${CreateDM3DIR_OMP}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${CreateDM3DIR_OMP}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${CreateDM3DIR_OMP}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${CreateDM3DIR_OMP}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${CreateDM3DIR_OMP}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${CreateDM3DIR_OMP}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${CreateDM3DIR_OMP}/variable.tmp
	@cat ${CreateDM3DIR_OMP}/variable.tmp >> ${CreateDM3DIR_OMP}/Makefile
	@cat ${CreateDM3DIR_OMP}/Makefile.template >> ${CreateDM3DIR_OMP}/Makefile
	$(MAKE) -C ${CreateDM3DIR_OMP}
	@install -d ${INSTALL_PATH}
	@install ${CreateDM3DIR_OMP}/${EXEC_NAME} ${INSTALL_PATH}


ReconstFC3DIR_OMP=${PARENT_DIR}/src/ReconstFC3/Src_OMP
ReconstFC3omp: 
	@if [ -f ${ReconstFC3DIR_OMP}/Makefile ]; then \
		rm ${ReconstFC3DIR_OMP}/Makefile; \
	fi
	@if [ -f ${ReconstFC3DIR_OMP}/variable.tmp ]; then \
		rm ${ReconstFC3DIR_OMP}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "reconst3_omp.x")
	@rm -f ${ReconstFC3DIR_OMP}/*.o ${ReconstFC3DIR_OMP}/*.mod ${ReconstFC3DIR_OMP}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qopenmp -qmkl=sequential")
else
	$(eval COMP_FLAG += "-fopenmp")
endif
	@echo 
	@touch ${ReconstFC3DIR_OMP}/Makefile
	@touch ${ReconstFC3DIR_OMP}/variable.tmp
	@echo "FC=$(FC)" >> ${ReconstFC3DIR_OMP}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${ReconstFC3DIR_OMP}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${ReconstFC3DIR_OMP}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${ReconstFC3DIR_OMP}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${ReconstFC3DIR_OMP}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${ReconstFC3DIR_OMP}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${ReconstFC3DIR_OMP}/variable.tmp
	@cat ${ReconstFC3DIR_OMP}/variable.tmp >> ${ReconstFC3DIR_OMP}/Makefile
	@cat ${ReconstFC3DIR_OMP}/Makefile.template >> ${ReconstFC3DIR_OMP}/Makefile
	$(MAKE) -C ${ReconstFC3DIR_OMP}
	@install -d ${INSTALL_PATH}
	@install ${ReconstFC3DIR_OMP}/${EXEC_NAME} ${INSTALL_PATH}


CreateDM4DIR_OMP=${PARENT_DIR}/src/CreateDM4/Src_OMP
CreateDM4omp: 
	@if [ -f ${CreateDM4DIR_OMP}/Makefile ]; then \
		rm ${CreateDM4DIR_OMP}/Makefile; \
	fi
	@if [ -f ${CreateDM4DIR_OMP}/variable.tmp ]; then \
		rm ${CreateDM4DIR_OMP}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "create_disp_mat4_omp.x")
	@rm -f ${CreateDM4DIR_OMP}/*.o ${CreateDM4DIR_OMP}/*.mod ${CreateDM4DIR_OMP}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qopenmp -qmkl=sequential")
else
	$(eval COMP_FLAG += "-fopenmp")
endif
	@echo 
	@touch ${CreateDM4DIR_OMP}/Makefile
	@touch ${CreateDM4DIR_OMP}/variable.tmp
	@echo "FC=$(FC)" >> ${CreateDM4DIR_OMP}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${CreateDM4DIR_OMP}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${CreateDM4DIR_OMP}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${CreateDM4DIR_OMP}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${CreateDM4DIR_OMP}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${CreateDM4DIR_OMP}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${CreateDM4DIR_OMP}/variable.tmp
	@cat ${CreateDM4DIR_OMP}/variable.tmp >> ${CreateDM4DIR_OMP}/Makefile
	@cat ${CreateDM4DIR_OMP}/Makefile.template >> ${CreateDM4DIR_OMP}/Makefile
	$(MAKE) -C ${CreateDM4DIR_OMP}
	@install -d ${INSTALL_PATH}
	@install ${CreateDM4DIR_OMP}/${EXEC_NAME} ${INSTALL_PATH}


ReconstFC4DIR_OMP=${PARENT_DIR}/src/ReconstFC4/Src_OMP
ReconstFC4omp: 
	@if [ -f ${ReconstFC4DIR_OMP}/Makefile ]; then \
		rm ${ReconstFC4DIR_OMP}/Makefile; \
	fi
	@if [ -f ${ReconstFC4DIR_OMP}/variable.tmp ]; then \
		rm ${ReconstFC4DIR_OMP}/variable.tmp; \
	fi
	$(eval EXEC_NAME = "reconst4_omp.x")
	@rm -f ${ReconstFC4DIR_OMP}/*.o ${ReconstFC4DIR_OMP}/*.mod ${ReconstFC4DIR_OMP}/${EXEC_NAME}; \
	$(eval LINKERFLAG = $(LINKERFLAG_COMMON))
	$(eval LIBPATH = $(LIBPATH_COMMON))
	$(eval INCPATH = $(INCPATH_COMMON))
	$(eval COMP_FLAG = $(COMP_FLAG_COMMON))
ifdef INTEL
	$(eval COMP_FLAG += "-qopenmp -qmkl=sequential")
else
	$(eval COMP_FLAG += "-fopenmp")
endif
	@echo 
	@touch ${ReconstFC4DIR_OMP}/Makefile
	@touch ${ReconstFC4DIR_OMP}/variable.tmp
	@echo "FC=$(FC)" >> ${ReconstFC4DIR_OMP}/variable.tmp
	@echo "MPIFC=$(MPIFC)" >> ${ReconstFC4DIR_OMP}/variable.tmp
	@echo "EXEC_NAME=$(EXEC_NAME)" >> ${ReconstFC4DIR_OMP}/variable.tmp
	@echo "LDFLAGS=$(LIBPATH)" >> ${ReconstFC4DIR_OMP}/variable.tmp
	@echo "LIBS=$(LINKERFLAG)" >> ${ReconstFC4DIR_OMP}/variable.tmp
	@echo "INCFLAGS=$(INCPATH)" >> ${ReconstFC4DIR_OMP}/variable.tmp
	@echo "FFLAGS=$(COMP_FLAG)" >> ${ReconstFC4DIR_OMP}/variable.tmp
	@cat ${ReconstFC4DIR_OMP}/variable.tmp >> ${ReconstFC4DIR_OMP}/Makefile
	@cat ${ReconstFC4DIR_OMP}/Makefile.template >> ${ReconstFC4DIR_OMP}/Makefile
	$(MAKE) -C ${ReconstFC4DIR_OMP}
	@install -d ${INSTALL_PATH}
	@install ${ReconstFC4DIR_OMP}/${EXEC_NAME} ${INSTALL_PATH}


clean:
	@rm -f ${FC2DIR}/*.o ${FC2DIR}/*.mod ${FC2DIR}/*.x ${FC2DIR}/Makefile ${FC2DIR}/variable.tmp
	@rm -f ${FC3DIR}/*.o ${FC3DIR}/*.mod ${FC3DIR}/*.x ${FC3DIR}/Makefile ${FC3DIR}/variable.tmp
	@rm -f ${FC4DIR}/*.o ${FC4DIR}/*.mod ${FC4DIR}/*.x ${FC4DIR}/Makefile ${FC4DIR}/variable.tmp
	@rm -f ${CreateDM2DIR}/*.o ${CreateDM2DIR}/*.mod ${CreateDM2DIR}/*.x ${CreateDM2DIR}/Makefile ${CreateDM2DIR}/variable.tmp
	@rm -f ${CreateDM3DIR}/*.o ${CreateDM3DIR}/*.mod ${CreateDM3DIR}/*.x ${CreateDM3DIR}/Makefile ${CreateDM3DIR}/variable.tmp
	@rm -f ${CreateDM4DIR}/*.o ${CreateDM4DIR}/*.mod ${CreateDM4DIR}/*.x ${CreateDM4DIR}/Makefile ${CreateDM4DIR}/variable.tmp
	@rm -f ${ReconstFC2DIR}/*.o ${ReconstFC2DIR}/*.mod ${ReconstFC2DIR}/*.x ${ReconstFC2DIR}/Makefile ${ReconstFC2DIR}/variable.tmp
	@rm -f ${ReconstFC3DIR}/*.o ${ReconstFC3DIR}/*.mod ${ReconstFC3DIR}/*.x ${ReconstFC3DIR}/Makefile ${ReconstFC3DIR}/variable.tmp
	@rm -f ${ReconstFC4DIR}/*.o ${ReconstFC4DIR}/*.mod ${ReconstFC4DIR}/*.x ${ReconstFC4DIR}/Makefile ${ReconstFC4DIR}/variable.tmp
	@rm -f ${RenormDIR}/*.o ${RenormDIR}/*.mod ${RenormDIR}/*.x ${RenormDIR}/Makefile ${RenormDIR}/variable.tmp
	@rm -f ${DispDIR}/*.o ${DispDIR}/*.mod ${DispDIR}/*.x ${DispDIR}/Makefile ${DispDIR}/variable.tmp
	@rm -f ${Long_DDDIR}/*.o ${Long_DDDIR}/*.mod ${Long_DDDIR}/*.x ${Long_DDDIR}/Makefile ${Long_DDDIR}/variable.tmp
	@rm -f ${LsqDIR}/*.o ${LsqDIR}/*.mod ${LsqDIR}/*.x ${LsqDIR}/Makefile ${LsqDIR}/variable.tmp
	@rm -f ${TSSDIR}/*.o ${TSSDIR}/*.mod ${TSSDIR}/*.x ${TSSDIR}/Makefile ${TSSDIR}/variable.tmp
	@rm -f ${FreeEnergyDIR}/*.o ${FreeEnergyDIR}/*.mod ${FreeEnergyDIR}/*.x ${FreeEnergyDIR}/Makefile ${FreeEnergyDIR}/variable.tmp
	@rm -f ${LineWidthDIR}/*.o ${LineWidthDIR}/*.mod ${LineWidthDIR}/*.x ${LineWidthDIR}/Makefile ${LineWidthDIR}/variable.tmp
	@rm -f ${kappaIterDIR}/*.o ${kappaIterDIR}/*.mod ${kappaIterDIR}/*.x ${kappaIterDIR}/Makefile ${kappaIterDIR}/variable.tmp
	@rm -f ${kappaDIR}/*.o ${kappaDIR}/*.mod ${kappaDIR}/*.x ${kappaDIR}/Makefile ${kappaDIR}/variable.tmp
	@rm -f ${CreateDM2DIR_OMP}/*.o ${CreateDM2DIR_OMP}/*.mod ${CreateDM2DIR_OMP}/*.x ${CreateDM2DIR_OMP}/Makefile ${CreateDM2DIR_OMP}/variable.tmp
	@rm -f ${CreateDM3DIR_OMP}/*.o ${CreateDM3DIR_OMP}/*.mod ${CreateDM3DIR_OMP}/*.x ${CreateDM3DIR_OMP}/Makefile ${CreateDM3DIR_OMP}/variable.tmp
	@rm -f ${CreateDM4DIR_OMP}/*.o ${CreateDM4DIR_OMP}/*.mod ${CreateDM4DIR_OMP}/*.x ${CreateDM4DIR_OMP}/Makefile ${CreateDM4DIR_OMP}/variable.tmp
	@rm -f ${ReconstFC2DIR_OMP}/*.o ${ReconstFC2DIR_OMP}/*.mod ${ReconstFC2DIR_OMP}/*.x ${ReconstFC2DIR_OMP}/Makefile ${ReconstFC2DIR_OMP}/variable.tmp
	@rm -f ${ReconstFC3DIR_OMP}/*.o ${ReconstFC3DIR_OMP}/*.mod ${ReconstFC3DIR_OMP}/*.x ${ReconstFC3DIR_OMP}/Makefile ${ReconstFC3DIR_OMP}/variable.tmp
	@rm -f ${ReconstFC4DIR_OMP}/*.o ${ReconstFC4DIR_OMP}/*.mod ${ReconstFC4DIR_OMP}/*.x ${ReconstFC4DIR_OMP}/Makefile ${ReconstFC4DIR_OMP}/variable.tmp

