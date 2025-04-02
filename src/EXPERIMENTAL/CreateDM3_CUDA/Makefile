FC=nvfortran
MPIFC=nvfortran

EXEC_NAME=create_dm3_CUDA.x

COMPILE_TIME=$(shell date +"%d-%m-%Y %H:%M:%S")
VER=\""0.1.0"\"

FFLAGS=
FFLINK=
LDFLAGS=
INCFLAGS=
LIBS =

HDF5_INC=/home/soham/Applications/Nvidia/hdf5_nvidia_new/lib/hdf5/include
HDF5_LIB=/home/soham/Applications/Nvidia/hdf5_nvidia_new/lib/hdf5/lib
#MKL_LIB=/home/soham/Applications/Nvidia/blas

#CUBLAS_INC=/home/soham/Applications/NVIDIA_HPC_SDK/Linux_x86_64/22.3/compilers/include
#CUBLAS_LIB=/home/soham/Applications/NVIDIA_HPC_SDK/Linux_x86_64/22.3/math_libs/11.6/targets/x86_64-linux/lib
#CUBLAS_LIB=/home/soham/Applications/NVIDIA_HPC_SDK/Linux_x86_64/22.3/math_libs/lib64/stubs

SRC= kinds.f90 constants.f90 parse_cmd_line.f90 hdf5_wrap_mod.f90 IndxChng_mod.f90 unit_mod.f90 FC3_data_mod.f90 NewFCtype.f90 DispMatFC3_mod.f90 main.f90
OBJ=$(SRC:.f90=.o)

FDEBUG=-O0 -debug -traceback -check -W1 -warn all -ftrapuv -qmkl=sequential -qopenmp

#NVFORTRAN specific flags
FOPTIM=-O3 -Wall -cpp -mp=multicore -cuda -gpu=cuda11.6 -gpu=cc86 -cudalib=cublas -D__COMPILE_TIME__=\""${COMPILE_TIME}"\" -D__SPR_VERSION__=${VER} #-Minfo

#FFLAGS+=$(FDEBUG)
FFLAGS+=$(FOPTIM)

INCFLAGS += -I${HDF5_INC} #-I${CUBLAS_INC}

LDFLAGS += -L${HDF5_LIB}
#LDFLAGS += -L${CUBLAS_LIB}

LIBS += -lhdf5 -lhdf5_fortran -lm -ldl #-lnvhpcwrapnvtx

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ $(INCFLAGS) -c $<

build: $(OBJ)
		$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXEC_NAME) $(OBJ) $(LIBS)

clean:
		rm -rf $(OBJ)
		rm -rf *.mod
		rm -rf $(EXEC_NAME)

