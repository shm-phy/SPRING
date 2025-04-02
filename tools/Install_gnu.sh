#!/bin/bash

set -e

#======================================== Version =======================================#
GCC_VER=12.2.0

cmake_ver=3.23.3
mpich_ver=4.0.2
fftw_ver=3.3.10
opencoarray_ver=2.10.0
OpenBLAS_ver=0.3.23
blas_ver=3.10.0
lapack_ver=3.10.1
zlib_ver=1.2.13
hdf5_ver=1_12_2
spglib_ver=1.16.5
#======================================== Version =======================================#

#====================================== Source Dir ======================================#
DIR=$PWD
PARENT_DIR=$DIR/ALL_SOURCE_GNU
mkdir --verbose -p ${PARENT_DIR}
#====================================== Source Dir ======================================#

#====================================== Install Dir =====================================#
cmake_install_dir=${HOME}/Applications/GNU/cmake
mpich_install_dir=${HOME}/Applications/GNU/mpich
fftw_install_dir=${HOME}/Applications/GNU/fftw
opencoarray_dir=${HOME}/Applications/GNU/opencoarrays
OpenBLAS_install_dir=${HOME}/Applications/GNU/OpenBLAS
blas_install_dir=${HOME}/Applications/GNU/blas
lapack_install_dir=${HOME}/Applications/GNU/lapack
zlib_install_dir=${HOME}/Applications/GNU/zlib
hdf5_install_dir=${HOME}/Applications/GNU/hdf5
spglib_install_dir=${HOME}/Applications/GNU/spglib
#====================================== Install Dir =====================================#

#========================================== URL =========================================#
cmake_url=https://github.com/Kitware/CMake/releases/download
mpich_url=https://www.mpich.org/static/downloads
fftw_url=https://www.fftw.org
opencoarray_url=https://github.com/sourceryinstitute/OpenCoarrays/archive/refs/tags
               #https://github.com/sourceryinstitute/OpenCoarrays/releases/download/2.10.0/OpenCoarrays-2.10.0.tar.gz
OpenBLAS_url=https://github.com/xianyi/OpenBLAS/releases/download
blas_url=http://www.netlib.org/blas
lapack_url=https://github.com/Reference-LAPACK/lapack/archive/refs/tags
zlib_url=https://zlib.net
        #https://zlib.net/fossils
hdf5_url=https://github.com/HDFGroup/hdf5/archive/refs/tags
spglib_url=https://github.com/spglib/spglib/archive/refs/tags
#========================================== URL =========================================#

#============================ Manually Installed GCC, CMAKE =============================#
my_gfortran=gfortran #${HOME}/Applications/gcc-${GCC_VER}/bin/gfortran
my_gcc=gcc #${HOME}/Applications/gcc-${GCC_VER}/bin/gcc
my_cpp=g++ #${HOME}/Applications/gcc-${GCC_VER}/bin/g++
#============================ Manually Installed GCC, CMAKE =============================#

#========================================= cmake ========================================#
cd $PARENT_DIR
mkdir --verbose -p ${cmake_install_dir}

wget --verbose ${cmake_url}/v${cmake_ver}/cmake-${cmake_ver}.tar.gz
tar -zxvf cmake-${cmake_ver}.tar.gz

export CC=$my_gcc
export CXX=$my_cpp
cd $PARENT_DIR/cmake-${cmake_ver}

./bootstrap --prefix=${cmake_install_dir}

make
make install

#remove
cd $PARENT_DIR
rm -rf $PARENT_DIR/cmake-${cmake_ver}

my_cmake=${cmake_install_dir}/bin/cmake
#========================================= cmake ========================================#

#========================================= mpich ========================================#
cd $PARENT_DIR
mkdir --verbose -p ${mpich_install_dir}

wget --verbose ${mpich_url}/${mpich_ver}/mpich-${mpich_ver}.tar.gz
tar -zxvf mpich-${mpich_ver}.tar.gz

cd $PARENT_DIR/mpich-${mpich_ver}

./configure CC=$my_gcc \
    CXX=$my_cpp \
    FC=$my_gfortran \
    F77=$my_gfortran \
    FFLAGS=-fallow-argument-mismatch \
    FCFLAGS=-fallow-argument-mismatch \
    --prefix=${mpich_install_dir} \
    --enable-fortran=all \
    --enable-shared \
    --enable-static 

make
make install

#remove
cd $PARENT_DIR
rm -rf $PARENT_DIR/mpich-${mpich_ver}

export LD_LIBRARY_PATH=${mpich_install_dir}'/lib':$LD_LIBRARY_PATH
export PATH=${mpich_install_dir}'/bin':$PATH

#========================================= mpich ========================================#


#========================================== fftw ========================================#

cd $PARENT_DIR
mkdir --verbose -p ${fftw_install_dir}

wget --verbose ${fftw_url}/fftw-${fftw_ver}.tar.gz
tar -zxvf fftw-${fftw_ver}.tar.gz

cd $PARENT_DIR/fftw-${fftw_ver}

./configure CC=$my_gcc F77=$my_gfortran MPICC=${mpich_install_dir}/bin/mpicc \
    LDFLAGS=-L${mpich_install_dir}/libs CPPFLAGS=-I${mpich_install_dir}/include \
    --prefix=${fftw_install_dir} \
    --enable-float --enable-sse \
    --with-g77-wrappers \
    --enable-shared --enable-static --enable-mpi

make
make check
make install

#remove
cd $PARENT_DIR
rm -rf $PARENT_DIR/fftw-${fftw_ver}

#========================================== fftw ========================================#

#========================================= spglib =======================================#
cd $PARENT_DIR

mkdir --verbose -p ${spglib_install_dir}

wget --verbose ${spglib_url}/v${spglib_ver}.tar.gz
tar -zxvf v${spglib_ver}.tar.gz

cd $PARENT_DIR/spglib-${spglib_ver}
sed -i 's/option(USE_OMP "Build with OpenMP support" ON)/option(USE_OMP "Build with OpenMP support" OFF)/' CMakeLists.txt

mkdir --verbose -p $PARENT_DIR/spglib-${spglib_ver}/_build

cd $PARENT_DIR/spglib-${spglib_ver}/_build

${my_cmake} -DCMAKE_INSTALL_PREFIX=${spglib_install_dir} ..
${my_cmake} --build .
${my_cmake} --install . --prefix ${spglib_install_dir}

#remove
cd $PARENT_DIR
rm -rf $PARENT_DIR/spglib-${spglib_ver}
#========================================= spglib =======================================#

#===================================== Opencoarrays =====================================#
cd $PARENT_DIR
mkdir --verbose -p ${opencoarray_dir}

wget --verbose ${opencoarray_url}/${opencoarray_ver}.tar.gz
tar -zxvf ${opencoarray_ver}.tar.gz


cd $PARENT_DIR/OpenCoarrays-${opencoarray_ver}

PATH="${mpich_install_dir}/bin:$PATH" LD_LIBRARY_PATH="${mpich_install_dir}/lib:$LD_LIBRARY_PATH" ./install.sh \
    --install-prefix ${opencoarray_dir} \
    --with-fortran ${my_gfortran} \
    --with-c ${my_gcc} \
    --with-cxx ${my_cpp} \
    --with-mpi ${mpich_install_dir} \
    --with-cmake ${my_cmake}

#remove
cd $PARENT_DIR
rm -rf $PARENT_DIR/OpenCoarrays-${opencoarray_ver}
#===================================== Opencoarrays =====================================#

#======================================= OpenBLAS =======================================#
cd $PARENT_DIR 
mkdir --verbose -p ${OpenBLAS_install_dir}

wget --verbose ${OpenBLAS_url}/v${OpenBLAS_ver}/OpenBLAS-${OpenBLAS_ver}.tar.gz

tar -zxvf OpenBLAS-${OpenBLAS_ver}.tar.gz

cd $PARENT_DIR/OpenBLAS-${OpenBLAS_ver}

make CC=${my_gcc} FC=${my_gfortran} BINARY=64 USE_THREAD=1 BUFFERSIZE=25 BUILD_LAPACK_DEPRECATED=1 NO_PARALLEL_MAKE=1 PREFIX=${OpenBLAS_install_dir}

make CC=${my_gcc} FC=${my_gfortran} BINARY=64 USE_THREAD=1 BUFFERSIZE=25 BUILD_LAPACK_DEPRECATED=1 NO_PARALLEL_MAKE=1 PREFIX=${OpenBLAS_install_dir} install

#make PREFIX=${OpenBLAS_install_dir} install

#remove
cd $PARENT_DIR
rm -rf $PARENT_DIR/OpenBLAS-${OpenBLAS_ver}
#======================================= OpenBLAS =======================================#


#========================================== BLAS ========================================#
cd $PARENT_DIR

mkdir --verbose -p ${blas_install_dir}

wget --verbose --no-check-certificate ${blas_url}/blas-${blas_ver}.tgz
tar -xvzf blas-${blas_ver}.tgz

cd $PARENT_DIR/BLAS-${blas_ver}
make

mv blas_LINUX.a libblas.a
cp libblas.a ${blas_install_dir}

#remove
cd $PARENT_DIR
rm -rf $PARENT_DIR/BLAS-${blas_ver}
#========================================== BLAS ========================================#

#======================================== LAPACK ========================================#
cd $PARENT_DIR

mkdir --verbose --verbose -p ${lapack_install_dir}

wget ${lapack_url}/v${lapack_ver}.tar.gz
tar -zxvf v${lapack_ver}.tar.gz

cd $PARENT_DIR/lapack-${lapack_ver}
mv make.inc.example make.inc

make

cp liblapack.a librefblas.a libtmglib.a ${lapack_install_dir}

#remove
cd $PARENT_DIR
rm -rf $PARENT_DIR/lapack-${lapack_ver}
#======================================== LAPACK ========================================#

export CC="$my_gcc -fPIC"
export CXX="$my_cpp -fPIC"
export CFLAGS='-O2'
export CXXFLAGS='-O2'
export F77="$my_gfortran -fPIC"
export FC="$my_gfortran -fPIC"
export F90="$my_gfortran -fPIC"
export FFLAGS='-O2'
export CPP="$my_gcc -E -fPIC"
export CXXCPP="$my_cpp -E -fPIC"

#======================================== ZLIB ==========================================#
cd $PARENT_DIR
mkdir --verbose -p ${zlib_install_dir}

#wget https://zlib.net/fossils/zlib-${zlib_ver}.tar.gz
#tar -zxvf zlib-${zlib_ver}.tar.gz

#cd $PARENT_DIR/zlib-${zlib_ver}
#./configure --prefix=${zlib_install_dir} \

wget --verbose ${zlib_url}/zlib-${zlib_ver}.tar.gz
tar -zxvf zlib-${zlib_ver}.tar.gz

cd $PARENT_DIR/zlib-${zlib_ver}
./configure --prefix=${zlib_install_dir}

make
make check
make install

#remove
cd $PARENT_DIR
rm -rf $PARENT_DIR/zlib-${zlib_ver}

export LD_LIBRARY_PATH=${zlib_install_dir}'/lib':$LD_LIBRARY_PATH
#======================================== ZLIB ==========================================#

#export CC='gcc '
#export CXX='g++ '
#export CFLAGS='-O2 '
#export CXXFLAGS='-O2 '
#export F77='gfortran '
#export FC='gfortran '
#export F90='gfortran '
#export FFLAGS='-O2 '
#export CPP='gcc -E '
#export CXXCPP='g++ -E '

#======================================== HDF5 ==========================================#

cd $PARENT_DIR
mkdir --verbose -p ${hdf5_install_dir}

wget --verbose ${hdf5_url}/hdf5-${hdf5_ver}.tar.gz
tar -zxvf hdf5-${hdf5_ver}.tar.gz

cd $PARENT_DIR/hdf5-hdf5-${hdf5_ver}

./configure --prefix=${hdf5_install_dir} \
    --with-zlib=${zlib_install_dir}/include,${zlib_install_dir}/lib \
    --enable-shared \
    --enable-static \
    --enable-cxx \
    --enable-fortran \
    --enable-fortran2003

make
make check
make install
#env "PATH=$PATH" "LD_LIBRARY_PATH=$LD_LIBRARY_PATH" make install

#remove 
cd $PARENT_DIR
rm -rf $PARENT_DIR/hdf5-hdf5-${hdf5_ver}
#======================================== HDF5 ==========================================#


#============================ Install paths in activate.bash ============================#
cd $DIR
FILE=activate_gnu.sh

echo 'export PATH="'${cmake_install_dir}'/bin:$PATH"' > $FILE

echo 'export LD_LIBRARY_PATH="'${mpich_install_dir}'/lib:$LD_LIBRARY_PATH"' >> $FILE
echo 'export PATH="'${mpich_install_dir}'/bin:$PATH"' >> $FILE

echo 'export PATH="'${fftw_install_dir}'/bin:$PATH"' >> $FILE
echo 'export LD_LIBRARY_PATH="'${fftw_install_dir}'/lib:$LD_LIBRARY_PATH"' >> $FILE

echo 'export LD_LIBRARY_PATH="'${opencoarray_dir}'/lib64:$LD_LIBRARY_PATH"' >> $FILE
echo 'export PATH="'${opencoarray_dir}'/bin:$PATH"' >> $FILE

echo 'export LD_LIBRARY_PATH="'${OpenBLAS_install_dir}'/lib:$LD_LIBRARY_PATH"' >> $FILE

echo 'export LD_LIBRARY_PATH="'${zlib_install_dir}'/lib:$LD_LIBRARY_PATH"' >> $FILE
echo 'export LD_LIBRARY_PATH="'${hdf5_install_dir}'/lib:$LD_LIBRARY_PATH"' >> $FILE

echo 'export LD_LIBRARY_PATH="'${hdf5_install_dir}'/lib:$LD_LIBRARY_PATH"' >> $FILE
echo 'export PATH="'${hdf5_install_dir}'/bin:$PATH"' >> $FILE

echo 'BLAS_INC_DIR='${blas_install_dir} >> $FILE
echo 'LAPACK_INC_DIR='${lapack_install_dir} >> $FILE
echo 'SPGDIR='${spglib_install_dir}'/lib' >> $FILE
#============================ Install paths in activate.bash ============================#

