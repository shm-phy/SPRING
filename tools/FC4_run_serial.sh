#!/bin/bash

EXEC=/home/soham/Applications/Intel/spring/create_disp_mat4.x
#VENDOR=gnu
VENDOR=intel

INFILE='input.nml'
TEMP=80
NBASIS=2
chunk=720

NPROCS=16
export FOR_COARRAY_NUM_IMAGES=16

for ((bb=1; bb <= $NBASIS; bb++))
do
    for ((a=1; a <= 3; a++))
    do
        if [ ${VENDOR} = "gnu" ]; then
            mpirun -np $NPROCS $EXEC -inp $INFILE -T $TEMP -mu $bb -alpha $a -DynSchd T -Nchunk $chunk -IdleImg1 T
        elif [ ${VENDOR} = "intel" ]; then
            $EXEC -inp $INFILE -T $TEMP -mu $bb -alpha $a -DynSchd T -Nchunk $chunk -IdleImg1 T
        fi 
    done
done

