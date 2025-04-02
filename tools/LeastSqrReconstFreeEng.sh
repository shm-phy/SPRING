#!/bin/bash

TEMP=300.0
prefix=NaCl
NBASIS=2
 
INFILE='input.nml'
chunk=81

#VENDOR=gnu
VENDOR=intel
NPROCS_2nd=1
NPROCS_3rd=16
NPROCS_4th=16

export MKL_NUM_THREADS=16

PREFIX=/home/soham/Applications/Intel/spring
# ==================================== Least-Square up to fourth order ===================================== #

echo "          "

EXEC_LSq=$PREFIX/LSq.x
OUT_LSQ=LSq.log

$EXEC_LSq -inp $INFILE -T $TEMP 2>&1 | tee $OUT_LSQ

echo "          "

# ==================================== Least-Square up to fourth order ===================================== #

# ==================================== Second-Order IFCs reconstruction ==================================== #

echo "          "

EXEC_RECONST2=$PREFIX/reconst2.x
OUT_RECONST2=reconst2.log

if [ ${VENDOR} = "gnu" ]; then
    mpirun -np $NPROCS_2nd $EXEC_RECONST2 -inp $INFILE -T $TEMP 2>&1 | tee $OUT_RECONST2
elif [ ${VENDOR} = "intel" ]; then
    export FOR_COARRAY_NUM_IMAGES=$NPROCS_2nd
    $EXEC_RECONST2 -inp $INFILE -T $TEMP 2>&1 | tee $OUT_RECONST2
fi


echo "          "

# ==================================== Second-Order IFCs reconstruction ==================================== #


# ==================================== Third-Order IFCs reconstruction ===================================== #

echo "          "

EXEC_RECONST3=$PREFIX/reconst3.x
OUT_RECONST3=reconst3.log

if [ ${VENDOR} = "gnu" ]; then
    mpirun -np $NPROCS_3rd $EXEC_RECONST3 -inp $INFILE -T $TEMP -ShengBTE T 2>&1 | tee $OUT_RECONST3
elif [ ${VENDOR} = "intel" ]; then
    export FOR_COARRAY_NUM_IMAGES=$NPROCS_3rd
    $EXEC_RECONST3 -inp $INFILE -T $TEMP -ShengBTE T 2>&1 | tee $OUT_RECONST3
fi

echo "          "

# ==================================== Third-Order IFCs reconstruction ===================================== #

# =================================== Fourth-Order IFCs reconstruction ===================================== #

echo "          "

EXEC_RECONST4=$PREFIX/reconst4.x
OUT_RECONST4=reconst4.log

if [ ${VENDOR} = "gnu" ]; then
    mpirun -np $NPROCS_4th $EXEC_RECONST4 -inp $INFILE -T $TEMP -ShengBTE T 2>&1 | tee $OUT_RECONST4
elif [ ${VENDOR} = "intel" ]; then
    export FOR_COARRAY_NUM_IMAGES=$NPROCS_4th
    $EXEC_RECONST4 -inp $INFILE -T $TEMP -ShengBTE T 2>&1 | tee $OUT_RECONST4
fi

echo "          "

# =================================== Fourth-Order IFCs reconstruction ===================================== #

# ============================================= Free-Energy ================================================ #

echo "          "

EXEC_FREEENG=$PREFIX/FreeEng.x
OUT_FREEENG=FreeEng.log

$EXEC_FREEENG -inp $INFILE -T $TEMP 2>&1 | tee $OUT_FREEENG

echo "          "

# ============================================= Free-Energy ================================================ #

