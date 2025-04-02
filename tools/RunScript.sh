#!/bin/bash

TEMP=300.0
prefix=NaCl
NBASIS=2
 
INFILE='input.nml'
chunk=81

VENDOR=gnu
#VENDOR=intel
NPROCS_2nd=1
NPROCS_3rd=6
#export FOR_COARRAY_NUM_IMAGES=6

#Directory where the nml file and displacement-force data-set is stored
Dir_Nml_file=/home/soham/Dropbox/Linux@Dropbox/PROJECTS/PHONON/Materials/NaCl/TSS_Type1/Pressure300K/Configure1/Iter2/IFC2
Dir_DispForce=/home/soham/Dropbox/Linux@Dropbox/PROJECTS/PHONON/Materials/NaCl/TSS_Type1/Pressure300K/Configure1/Iter2/IFC2

#Remote directories
CURRENT_CONFIG=Configure1/Iter3
REMOTE_DIR_PREFIX=/home/soham/PROJECTS/FourthOrder/NaCl/Pressure300K

# Directories where 2nd and 3rd order IFC information are stored
Dir_2nd=/home/soham/Dropbox/Linux@Dropbox/PROJECTS/PHONON/TEST/TEST_COMBILE/Dir_2nd
Dir_3rd=/home/soham/Dropbox/Linux@Dropbox/PROJECTS/PHONON/Materials/NaCl/TSS_Type1/Pressure300K/Configure1/Iter2/IFC3
Dir_4th=


# Parent directory of FD 2nd FD 3rd
Dir_FD_2nd=$PWD

Dir_FD_3rd=$Dir_FD_2nd
Dir_LongEw=$Dir_FD_2nd
Dir_LSq=$Dir_FD_2nd
Dir_Disp=$Dir_FD_2nd
Dir_SeconWith4th=$Dir_FD_2nd

#Directory where the python plot script is stored
Dir_plot_py=/home/soham/Dropbox/Linux@Dropbox/PROJECTS/PHONON/SPRING/tools

# =============================== Second-Order Displacement Matrix creation ================================ #

EXEC_CDF2=create_disp_mat2.x


mkdir -p --verbose ${Dir_FD_2nd}/IFC2

cp -v $Dir_DispForce/disp_forc_${prefix}${TEMP}K.h5 ${Dir_FD_2nd}/IFC2

cp -v $Dir_Nml_file/$INFILE ${Dir_FD_2nd}/IFC2

cp -v $Dir_2nd/FC_2nd_common_${prefix}_F.h5 ${Dir_FD_2nd}/IFC2

cd ${Dir_FD_2nd}/IFC2


for ((bb=1; bb <= $NBASIS; bb++))
do
    for ((a=1; a <= 3; a++))
    do
        if [ ${VENDOR} = "gnu" ]; then
            mpirun -np $NPROCS_2nd $EXEC_CDF2 -inp $INFILE -T $TEMP -mu $bb -alpha $a -DynSchd T -Nchunk $chunk -IdleImg1 T
        elif [ ${VENDOR} = "intel" ]; then
            export FOR_COARRAY_NUM_IMAGES=$NPROCS_2nd
            $EXEC_CDF2 -inp $INFILE -T $TEMP -mu $bb -alpha $a -DynSchd T -Nchunk $chunk -IdleImg1 T
        fi
    done
done

echo "          "
rm --verbose ${Dir_FD_2nd}/IFC2/FDpart_2nd_${prefix}${TEMP}K.h5
rm --verbose ${Dir_FD_2nd}/IFC2/FC_2nd_common_${prefix}_F.h5

# =============================== Second-Order Displacement Matrix creation ================================ #

# ============================== Second-Order Dipole-dipole IFCs calculation =============================== #

echo "          "

EXEC_longEw=longEw.x

mkdir -p --verbose ${Dir_LongEw}/Ewald

cp -v $Dir_DispForce/disp_forc_${prefix}${TEMP}K.h5 ${Dir_LongEw}/Ewald

cp -v $Dir_Nml_file/$INFILE ${Dir_LongEw}/Ewald

cd ${Dir_LongEw}/Ewald

$EXEC_longEw -inp $INFILE -T $TEMP

echo "          "

# ============================== Second-Order Dipole-dipole IFCs calculation =============================== #


# ================================ Third-Order Displacement Matrix creation ================================ #

echo "          "

EXEC_CDF3=create_disp_mat3.x


mkdir -p --verbose ${Dir_FD_3rd}/IFC3

cp -v $Dir_DispForce/disp_forc_${prefix}${TEMP}K.h5 ${Dir_FD_3rd}/IFC3

cp -v $Dir_Nml_file/$INFILE ${Dir_FD_3rd}/IFC3

cp -v $Dir_3rd/FC_3rd_common_${prefix}_F.h5 ${Dir_FD_3rd}/IFC3

cd ${Dir_FD_3rd}/IFC3

for ((bb=1; bb <= $NBASIS; bb++))
do
    for ((a=1; a <= 3; a++))
    do
        if [ ${VENDOR} = "gnu" ]; then
            mpirun -np $NPROCS_3rd $EXEC_CDF3 -inp $INFILE -T $TEMP -mu $bb -alpha $a -DynSchd T -Nchunk $chunk -IdleImg1 T
        elif [ ${VENDOR} = "intel" ]; then
            export FOR_COARRAY_NUM_IMAGES=$NPROCS_3rd
            $EXEC_CDF3 -inp $INFILE -T $TEMP -mu $bb -alpha $a -DynSchd T -Nchunk $chunk -IdleImg1 T
        fi
    done
done


rm --verbose ${Dir_FD_3rd}/IFC3/FDpart_3rd_${prefix}${TEMP}K.h5
rm --verbose ${Dir_FD_3rd}/IFC3/FC_3rd_common_${prefix}_F.h5

echo "          "

# ================================ Third-Order Displacement Matrix creation ================================ #

# ================================ Least Square and Reconstruction of IFC2 ================================= #

echo "          "

EXEC_LSq=LSq.x
EXEC_RECONST2=reconst2.x

mkdir -p --verbose ${Dir_LSq}/Least_Sqr

mv -v ${Dir_FD_2nd}/IFC2/FD_2nd_${prefix}${TEMP}K.h5 ${Dir_LSq}/Least_Sqr
mv -v ${Dir_FD_3rd}/IFC3/FD_3rd_${prefix}${TEMP}K.h5 ${Dir_LSq}/Least_Sqr
mv -v ${Dir_LongEw}/Ewald/FC2_dd_${prefix}${TEMP}K.h5 ${Dir_LSq}/Least_Sqr

cp -v $Dir_Nml_file/$INFILE ${Dir_LSq}/Least_Sqr

cd ${Dir_LSq}/Least_Sqr

$EXEC_LSq -inp $INFILE -T $TEMP

echo "          "
cp -v $Dir_2nd/FC_2nd_common_${prefix}_F.h5 ${Dir_LSq}/Least_Sqr

$EXEC_RECONST2 -inp $INFILE -T $TEMP
echo "          "

rm --verbose ${Dir_LSq}/Least_Sqr/FC_2nd_common_${prefix}_F.h5

echo "          "

# ================================ Least Square and Reconstruction of IFC2 ================================= #


# =================================== ssh the displacement matrix files ==================================== #

echo "          "

ssh soham@10.43.21.23 "mkdir -p --verbose ${REMOTE_DIR_PREFIX}/${CURRENT_CONFIG}/Least_Sqr"

scp FD_2nd_${prefix}${TEMP}K.h5 FD_3rd_${prefix}${TEMP}K.h5 FC2_dd_${prefix}${TEMP}K.h5 $INFILE soham@10.43.21.23:${REMOTE_DIR_PREFIX}/${CURRENT_CONFIG}/Least_Sqr

echo "          "

# =================================== ssh the displacement matrix files ==================================== #


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

EXEC_DISP=disp.x

mkdir -p --verbose ${Dir_Disp}/Disp

cp -v ${Dir_LSq}/Least_Sqr/FC2nd_${prefix}${TEMP}K_F.h5 ${Dir_Disp}/Disp

cp -v $Dir_plot_py/plot.py ${Dir_Disp}/Disp
cp -v $Dir_Nml_file/$INFILE ${Dir_Disp}/Disp

cd ${Dir_Disp}/Disp

$EXEC_DISP -inp $INFILE -T $TEMP

python3 plot.py -inp $INFILE -T $TEMP

echo "          "

# =============================================== Dispersion =============================================== #

