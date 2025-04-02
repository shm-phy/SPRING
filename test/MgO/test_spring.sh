#!/bin/bash

set -e

# CLEAN: rm -r CreateDM2/ CreateDM3/ Disp/ Ewald/ IFC2/ IFC3/ Least_Sqr/ SecondOrderWith4th/

cd $PWD
DIR_CURRENT_CONFIG=$PWD
CURRENT_CONFIG=Configure8

#Material specific parameters
TEMP=300.0
prefix=MgO
NBASIS=2
 
INFILE='input.nml'
chunk=81

# Parameter to control processors of parallel execution
VENDOR=gnu
#VENDOR=intel

nthreads_2nd_fc=1
nthreads_3rd_fc=4

NPROCS_2nd=1
NPROCS_3rd=6

nthreads_lsq=4

NPROCS_reconst2nd=1
NPROCS_reconst3rd=4

# Executalble
PWD_EXEC1=${DIR_CURRENT_CONFIG}/../../bin
PWD_EXEC2=${PWD_EXEC1}

#Directory where the nml file and displacement-force data-set is stored
Dir_Nml_file=${DIR_CURRENT_CONFIG}/CreateSupCell
Dir_DispForce=${DIR_CURRENT_CONFIG}/ReadForce

# Directories where 2nd and 3rd order IFC information are stored
Dir_2nd=${DIR_CURRENT_CONFIG}/IFC2
Dir_3rd=${DIR_CURRENT_CONFIG}/IFC3
Dir_4th=


# Parent directory of FD 2nd FD 3rd
Dir_FD_2nd=${DIR_CURRENT_CONFIG}

Dir_FD_3rd=$Dir_FD_2nd
Dir_LongEw=$Dir_FD_2nd
Dir_LSq=$Dir_FD_2nd
Dir_Disp=$Dir_FD_2nd
Dir_SeconWith4th=$Dir_FD_2nd

#Directory where the python plot script is stored
Dir_plot_py=/home/soham/Dropbox/Linux@Dropbox/PROJECTS/PHONON/SPRING/tools

# ================================== Symmetry reduction of 2nd-Order IFCs ================================== #

EXEC_IFC2=${PWD_EXEC1}/FC2.x

mkdir -p --verbose ${Dir_2nd}

cp -v $Dir_Nml_file/$INFILE ${Dir_2nd}

cd ${Dir_2nd}

export MKL_NUM_THREADS=${nthreads_2nd_fc}
$EXEC_IFC2 -inp $INFILE -AdvncOut T  #2>&1 | tee IFC2_${prefix}.log

rm --verbose ${Dir_2nd}/InterMatChunk_2nd_${prefix}.h5

# ================================== Symmetry reduction of 2nd-Order IFCs ================================== #

# ================================== Symmetry reduction of 3rd-Order IFCs ================================== #

EXEC_IFC3=${PWD_EXEC1}/FC3.x

mkdir -p --verbose ${Dir_3rd}

cp -v $Dir_Nml_file/$INFILE ${Dir_3rd}

cd ${Dir_3rd}

export MKL_NUM_THREADS=${nthreads_3rd_fc}
$EXEC_IFC3 -inp $INFILE -AdvncOut T #2>&1 | tee IFC3_${prefix}.log

rm --verbose ${Dir_3rd}/InterMatChunk_3rd_${prefix}.h5

# ================================== Symmetry reduction of 3rd-Order IFCs ================================== #

# =============================== Second-Order Displacement Matrix creation ================================ #

EXEC_CDF2=${PWD_EXEC2}/create_disp_mat2.x


mkdir -p --verbose ${Dir_FD_2nd}/CreateDM2

cp -v $Dir_DispForce/disp_forc_${prefix}${TEMP}K.h5 ${Dir_FD_2nd}/CreateDM2

cp -v $Dir_Nml_file/$INFILE ${Dir_FD_2nd}/CreateDM2

cp -v $Dir_2nd/FC_2nd_common_${prefix}_F.h5 ${Dir_FD_2nd}/CreateDM2

cd ${Dir_FD_2nd}/CreateDM2
touch ${Dir_FD_2nd}/CreateDM2/CreateDM2_${prefix}.log


for ((bb=1; bb <= $NBASIS; bb++))
do
    for ((a=1; a <= 3; a++))
    do
        if [ ${VENDOR} = "gnu" ]; then
            mpirun -np $NPROCS_2nd $EXEC_CDF2 -inp $INFILE -T $TEMP -mu $bb -alpha $a -DynSchd T -Nchunk $chunk -IdleImg1 T #2>&1 | tee -a CreateDM2_${prefix}.log
        elif [ ${VENDOR} = "intel" ]; then
            export MKL_NUM_THREADS=1
            export FOR_COARRAY_NUM_IMAGES=$NPROCS_2nd
            $EXEC_CDF2 -inp $INFILE -T $TEMP -mu $bb -alpha $a -DynSchd T -Nchunk $chunk -IdleImg1 T #2>&1 | tee -a CreateDM2_${prefix}.log
        fi
    done
done

echo "          "
rm --verbose ${Dir_FD_2nd}/CreateDM2/FDpart_2nd_${prefix}${TEMP}K.h5
rm --verbose ${Dir_FD_2nd}/CreateDM2/FC_2nd_common_${prefix}_F.h5

# =============================== Second-Order Displacement Matrix creation ================================ #

# ============================== Second-Order Dipole-dipole IFCs calculation =============================== #

echo "          "

EXEC_longEw=${PWD_EXEC1}/longEw.x

mkdir -p --verbose ${Dir_LongEw}/Ewald

cp -v $Dir_DispForce/disp_forc_${prefix}${TEMP}K.h5 ${Dir_LongEw}/Ewald

cp -v $Dir_Nml_file/$INFILE ${Dir_LongEw}/Ewald

cd ${Dir_LongEw}/Ewald

export MKL_NUM_THREADS=1
$EXEC_longEw -inp $INFILE -T $TEMP #2>&1 | tee LongDD_${prefix}.log

echo "          "

# ============================== Second-Order Dipole-dipole IFCs calculation =============================== #


# ================================ Third-Order Displacement Matrix creation ================================ #

echo "          "

EXEC_CDF3=${PWD_EXEC2}/create_disp_mat3.x


mkdir -p --verbose ${Dir_FD_3rd}/CreateDM3

cp -v $Dir_DispForce/disp_forc_${prefix}${TEMP}K.h5 ${Dir_FD_3rd}/CreateDM3

cp -v $Dir_Nml_file/$INFILE ${Dir_FD_3rd}/CreateDM3

cp -v $Dir_3rd/FC_3rd_common_${prefix}_F.h5 ${Dir_FD_3rd}/CreateDM3

cd ${Dir_FD_3rd}/CreateDM3
touch ${Dir_FD_3rd}/CreateDM3/CreateDM3_${prefix}.log

for ((bb=1; bb <= $NBASIS; bb++))
do
    for ((a=1; a <= 3; a++))
    do
        if [ ${VENDOR} = "gnu" ]; then
            mpirun -np $NPROCS_3rd $EXEC_CDF3 -inp $INFILE -T $TEMP -mu $bb -alpha $a -DynSchd T -Nchunk $chunk -IdleImg1 T #2>&1 | tee -a CreateDM3_${prefix}.log
        elif [ ${VENDOR} = "intel" ]; then
            export MKL_NUM_THREADS=1
            export FOR_COARRAY_NUM_IMAGES=$NPROCS_3rd
            $EXEC_CDF3 -inp $INFILE -T $TEMP -mu $bb -alpha $a -DynSchd T -Nchunk $chunk -IdleImg1 T #2>&1 | tee -a CreateDM3_${prefix}.log
        fi
    done
done


rm --verbose ${Dir_FD_3rd}/CreateDM3/FDpart_3rd_${prefix}${TEMP}K.h5
rm --verbose ${Dir_FD_3rd}/CreateDM3/FC_3rd_common_${prefix}_F.h5

echo "          "

# ================================ Third-Order Displacement Matrix creation ================================ #

# ================================ Least Square and Reconstruction of IFC2 ================================= #

echo "          "

EXEC_LSq=${PWD_EXEC1}/LSq.x
EXEC_RECONST2=${PWD_EXEC2}/reconst2.x
EXEC_RECONST3=${PWD_EXEC2}/reconst3.x

mkdir -p --verbose ${Dir_LSq}/Least_Sqr

mv -v ${Dir_FD_2nd}/CreateDM2/FD_2nd_${prefix}${TEMP}K.h5 ${Dir_LSq}/Least_Sqr
mv -v ${Dir_FD_3rd}/CreateDM3/FD_3rd_${prefix}${TEMP}K.h5 ${Dir_LSq}/Least_Sqr
mv -v ${Dir_LongEw}/Ewald/FC2_dd_${prefix}${TEMP}K.h5 ${Dir_LSq}/Least_Sqr

cp -v $Dir_Nml_file/$INFILE ${Dir_LSq}/Least_Sqr

cd ${Dir_LSq}/Least_Sqr

export MKL_NUM_THREADS=${nthreads_lsq}
$EXEC_LSq -inp $INFILE -T $TEMP -CutOff -1 #2>&1 | tee LSq_${prefix}.log

# ------------------------ Reconstruct IFC2 ------------------------ #
echo "          "
cp -v $Dir_2nd/FC_2nd_common_${prefix}_F.h5 ${Dir_LSq}/Least_Sqr

if [ ${VENDOR} = "gnu" ]; then
    mpirun -np $NPROCS_reconst2nd $EXEC_RECONST2 -inp $INFILE -T $TEMP #2>&1 | tee Reconst2_${prefix}.log
elif [ ${VENDOR} = "intel" ]; then
    export MKL_NUM_THREADS=1
    export FOR_COARRAY_NUM_IMAGES=${NPROCS_reconst2nd}
    $EXEC_RECONST2 -inp $INFILE -T $TEMP #2>&1 | tee Reconst2_${prefix}.log
fi


echo "          "
rm --verbose ${Dir_LSq}/Least_Sqr/FC_2nd_common_${prefix}_F.h5
# ------------------------ Reconstruct IFC2 ------------------------ #

# ------------------------ Reconstruct IFC3 ------------------------ #
echo "          "
cp -v $Dir_3rd/FC_3rd_common_${prefix}_F.h5 ${Dir_LSq}/Least_Sqr


if [ ${VENDOR} = "gnu" ]; then
    mpirun -np $NPROCS_reconst3rd $EXEC_RECONST3 -inp $INFILE -T $TEMP #2>&1 | tee Reconst3_${prefix}.log
elif [ ${VENDOR} = "intel" ]; then
    export MKL_NUM_THREADS=1
    export FOR_COARRAY_NUM_IMAGES=${NPROCS_reconst3rd}
    export FOR_COARRAY_NUM_IMAGES=${NPROCS_reconst3rd}
    $EXEC_RECONST3 -inp $INFILE -T $TEMP #2>&1 | tee Reconst3_${prefix}.log
fi

echo "          "
rm --verbose ${Dir_LSq}/Least_Sqr/FC_3rd_common_${prefix}_F.h5
# ------------------------ Reconstruct IFC3 ------------------------ #

echo "          "

# ================================ Least Square and Reconstruction of IFC2 ================================= #


# =========================================== Second-Order with 4th ======================================== #

echo "          "

mkdir -p --verbose ${Dir_SeconWith4th}/SecondOrderWith4th

cp -v $Dir_Nml_file/$INFILE ${Dir_SeconWith4th}/SecondOrderWith4th
cp -v $Dir_plot_py/plot.py ${Dir_SeconWith4th}/SecondOrderWith4th

cp -v $Dir_2nd/FC_2nd_common_${prefix}_F.h5 ${Dir_SeconWith4th}/SecondOrderWith4th

echo "          "

# =========================================== Second-Order with 4th ======================================== #

# =============================================== Dispersion =============================================== #

echo "          "

EXEC_DISP=${PWD_EXEC1}/disp.x

mkdir -p --verbose ${Dir_Disp}/Disp

cp -v ${Dir_LSq}/Least_Sqr/FC2nd_${prefix}${TEMP}K_F.h5 ${Dir_Disp}/Disp

cp -v $Dir_plot_py/plot.py ${Dir_Disp}/Disp
cp -v $Dir_Nml_file/$INFILE ${Dir_Disp}/Disp

cd ${Dir_Disp}/Disp

export MKL_NUM_THREADS=1
$EXEC_DISP -inp $INFILE -T $TEMP #2>&1 | tee -a Disp_${prefix}.log

python3 plot.py -inp $INFILE -T $TEMP

echo "          "

# =============================================== Dispersion =============================================== #

