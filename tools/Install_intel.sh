#!/bin/bash

set -e

MY_CC=icx
MY_CXX=icpx
MY_FC=ifort

#======================================== Version =======================================#
cmake_ver=3.23.3
fftw_ver=3.3.10
zlib_ver=1.2.13
hdf5_ver=1_12_2
spglib_ver=1.16.5

NETCDF_C_VER=4.9.0
NETCDF_F_VER=4.5.4
#======================================== Version =======================================#

#====================================== Source Dir ======================================#
DIR=$PWD
PARENT_DIR=$DIR/ALL_SOURCE_INTEL
mkdir --verbose -p ${PARENT_DIR}
#====================================== Source Dir ======================================#

#====================================== Install Dir =====================================#
cmake_install_dir=${HOME}/Applications/Intel/cmake
fftw_install_dir=${HOME}/Applications/Intel/fftw
zlib_install_dir=${HOME}/Applications/Intel/zlib
hdf5_install_dir=${HOME}/Applications/Intel/hdf5
spglib_install_dir=${HOME}/Applications/Intel/spglib

netcdfC_install_dir=${HOME}/Applications/Intel/netcdf-c
netcdfF_install_dir=${HOME}/Applications/Intel/netcdf-fortran
#====================================== Install Dir =====================================#

#========================================== URL =========================================#
cmake_url=https://github.com/Kitware/CMake/releases/download
fftw_url=https://www.fftw.org
zlib_url=https://zlib.net
        #https://zlib.net/fossils
hdf5_url=https://github.com/HDFGroup/hdf5/archive/refs/tags
spglib_url=https://github.com/spglib/spglib/archive/refs/tags

NETCDF_C_URL=https://github.com/Unidata/netcdf-c/archive/refs/tags
NETCDF_F_URL=https://github.com/Unidata/netcdf-fortran/archive/refs/tags
#========================================== URL =========================================#

#========================================= cmake ========================================#
cd $PARENT_DIR
mkdir --verbose -p ${cmake_install_dir}

wget --verbose ${cmake_url}/v${cmake_ver}/cmake-${cmake_ver}.tar.gz
tar -zxvf cmake-${cmake_ver}.tar.gz

cd $PARENT_DIR/cmake-${cmake_ver}

export CC=$MY_CC
export CXX=$MY_CXX
./bootstrap --prefix=${cmake_install_dir}

make
make install

#remove
cd $PARENT_DIR
rm -rf $PARENT_DIR/cmake-${cmake_ver}

my_cmake=${cmake_install_dir}/bin/cmake

#========================================= cmake ========================================#

#========================================== fftw ========================================#

cd $PARENT_DIR
mkdir --verbose -p ${fftw_install_dir}

wget --verbose ${fftw_url}/fftw-${fftw_ver}.tar.gz
tar -zxvf fftw-${fftw_ver}.tar.gz

cd $PARENT_DIR/fftw-${fftw_ver}

./configure CC=$MY_CC F77=$MY_FC MPICC="mpiicc -CC=$MY_CC"\
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

${my_cmake} -DCMAKE_INSTALL_PREFIX=${spglib_install_dir} \
    -DCMAKE_C_COMPILER=$MY_CC -DCMAKE_CXX_COMPILER=$MY_CXX -DCMAKE_Fortran_COMPILER=$MY_FC ..
${my_cmake} --build .
${my_cmake} --install . --prefix ${spglib_install_dir}

#remove
cd $PARENT_DIR
rm -rf $PARENT_DIR/spglib-${spglib_ver}

#========================================= spglib =======================================#

MY_CC=icx
MY_CXX=icpx
MY_FC=ifort

if [ $MY_CC == "icc" ] && [ $MY_CXX == "icpc" ]; then

    echo "Using icc and icpc compiler ..."
    export CC="$MY_CC"
    export CXX="$MY_CXX"
    export CFLAGS='-O2 -xHost -ip -no-prec-div -static-intel'
    export CXXFLAGS='-O2 -xHost -ip -no-prec-div -static-intel'
    export CPP="$MY_CC -E"
    export CXXCPP="$MY_CXX -E"

elif [ $MY_CC == "icx" ] && [ $MY_CXX == "icpx" ]; then

    echo "Using icx and icpx compiler ..."
    export CC="$MY_CC -fPIC"
    export CXX="$MY_CXX -fPIC"
    export CFLAGS='-O2 -xHost'
    export CXXFLAGS='-O2 -xHost'
    export CPP="$MY_CC -E"
    export CXXCPP="$MY_CXX -E"

else

    echo "Unsupported combination of compilers"
    exit 1

fi

export F77="$MY_FC"
export FC="$MY_FC"
export F90="$MY_FC"
export FFLAGS='-O2 -xHost -ip -no-prec-div -static-intel'

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

#Set the LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${zlib_install_dir}'/lib':$LD_LIBRARY_PATH

#======================================== ZLIB ==========================================#

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

#Set the LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${hdf5_install_dir}'/lib':$LD_LIBRARY_PATH

#======================================== HDF5 ==========================================#

#NetCDF test fails for the new intel compiler. So go back to calssic compilers
MY_CC=icc
MY_CXX=icpc

if [ $MY_CC == "icc" ] && [ $MY_CXX == "icpc" ]; then

    echo "Using icc and icpc compiler ..."
    export CC="$MY_CC"
    export CXX="$MY_CXX"
    export CFLAGS='-O2 -xHost -ip -no-prec-div -static-intel'
    export CXXFLAGS='-O2 -xHost -ip -no-prec-div -static-intel'
    export CPP="$MY_CC -E"
    export CXXCPP="$MY_CXX -E"

elif [ $MY_CC == "icx" ] && [ $MY_CXX == "icpx" ]; then

    echo "Using icx and icpx compiler ..."
    export CC="$MY_CC -fPIC"
    export CXX="$MY_CXX -fPIC"
    export CFLAGS='-O2 -xHost'
    export CXXFLAGS='-O2 -xHost'
    export CPP="$MY_CC -E"
    export CXXCPP="$MY_CXX -E"

else

    echo "Unsupported combination of compilers"
    exit 1

fi

#======================================= NETCDF-C =======================================#

cd $PARENT_DIR
mkdir --verbose -p ${netcdfC_install_dir}

wget --verbose ${NETCDF_C_URL}/v${NETCDF_C_VER}.tar.gz

tar -xvf v${NETCDF_C_VER}.tar.gz

cd $PARENT_DIR/netcdf-c-${NETCDF_C_VER}

CPPFLAGS="-std=gnu99 -I${zlib_install_dir}/include -I${hdf5_install_dir}/include" LDFLAGS="-L${zlib_install_dir}/lib -L${hdf5_install_dir}/lib" ./configure --prefix=${netcdfC_install_dir} --enable-shared --enable-static --enable-netcdf-4 --disable-dap

make
make check
make install

#remove
cd $PARENT_DIR
rm -rf $PARENT_DIR/netcdf-c-${NETCDF_C_VER}

#Set the LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${netcdfC_install_dir}'/lib':$LD_LIBRARY_PATH

#======================================= NETCDF-C =======================================#

#==================================== NETCDF-FORTRAN ====================================#

cd $PARENT_DIR
mkdir --verbose -p ${netcdfF_install_dir}

wget --verbose ${NETCDF_F_URL}/v${NETCDF_F_VER}.tar.gz

tar -xvf v${NETCDF_F_VER}.tar.gz

cd $PARENT_DIR/netcdf-fortran-${NETCDF_F_VER}

NCDIR=${netcdfC_install_dir} NFDIR=${netcdfF_install_dir} CPPFLAGS=-I${netcdfC_install_dir}/include LDFLAGS=-L${netcdfC_install_dir}/lib ./configure --prefix=${netcdfF_install_dir} --enable-static --enable-shared 

make
make check
make install

#remove
cd $PARENT_DIR
rm -rf $PARENT_DIR/netcdf-fortran-${NETCDF_F_VER}

#==================================== NETCDF-FORTRAN ====================================#

#============================ Install paths in activate.bash ============================#
cd $DIR
FILE=activate_intel.sh
echo 'export PATH="'${cmake_install_dir}'/bin:$PATH"' > $FILE
echo 'export PATH="'${fftw_install_dir}'/bin:$PATH"' >> $FILE
echo 'export LD_LIBRARY_PATH="'${fftw_install_dir}'/lib:$LD_LIBRARY_PATH"' >> $FILE
echo 'export LD_LIBRARY_PATH="'${zlib_install_dir}'/lib:$LD_LIBRARY_PATH"' >> $FILE
echo 'export LD_LIBRARY_PATH="'${hdf5_install_dir}'/lib:$LD_LIBRARY_PATH"' >> $FILE
echo 'export PATH="'${hdf5_install_dir}'/bin:$PATH"' >> $FILE
echo 'export LD_LIBRARY_PATH="'${netcdfC_install_dir}'/lib:$LD_LIBRARY_PATH"' >> $FILE
echo 'export LD_LIBRARY_PATH="'${netcdfF_install_dir}'/lib:$LD_LIBRARY_PATH"' >> $FILE
echo 'SPGDIR='${spglib_install_dir}'/lib64' >> $FILE
#============================ Install paths in activate.bash ============================#

