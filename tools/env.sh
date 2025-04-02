
#======================================= GCC Environments =======================================#
export LD_LIBRARY_PATH="/opt/hdf5_gcc/zlib_gcc/lib:$LD_LIBRARY_PATH"

export LD_LIBRARY_PATH="/opt/hdf5_gcc/hdf5_1.12.2/lib:$LD_LIBRARY_PATH"
export PATH="/opt/mpich/bin:$PATH"
export LD_LIBRARY_PATH="/opt/mpich/lib:$LD_LIBRARY_PATH"
export PATH="/opt/opencoarrays/bin:$PATH"
export LD_LIBRARY_PATH="/opt/opencoarrays/lib:$LD_LIBRARY_PATH"

export OPENBLAS_NUM_THREADS=4

#======================================= GCC Environments =======================================#

#====================================== Intel Environments ======================================#
export LD_LIBRARY_PATH="/opt/siesta-v4.1.5/DEPENDENCY/zlib_intel/lib:$LD_LIBRARY_PATH"

export LD_LIBRARY_PATH="/opt/siesta-v4.1.5/DEPENDENCY/hdf5_intel/lib:$LD_LIBRARY_PATH"
export MKL_NUM_THREADS=6
export FOR_COARRAY_NUM_IMAGES=6

export I_MPI_PIN=1
#====================================== Intel Environments ======================================#

#====================================== OpenMP Environment ======================================#
export OMP_NUM_THREADS=6
export OMP_STACKSIZE="20M"
export OMP_PROC_BIND=TRUE
#====================================== OpenMP Environment ======================================#

