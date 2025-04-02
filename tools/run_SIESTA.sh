#!/bin/bash

PBS_O_WORKDIR=${PWD}

#Job specific data and files
prefix=NaCl
TEMP=80.0

STR=0
END=60

SUP_DIR=SnapShots;
OUT_ALL=OUT;

BASE_FILE=input.fdf;

cd $PBS_O_WORKDIR
mkdir -p $PBS_O_WORKDIR/${OUT_ALL};

for i in $(eval echo "{$STR..$END}"); do 

    cd $PBS_O_WORKDIR

    cp $PBS_O_WORKDIR/${SUP_DIR}/${prefix}_T${TEMP}K_${i}.fdf $PBS_O_WORKDIR/system.fdf;
    outfile=${prefix}_sup_${i}.out;

    outdir=$PBS_O_WORKDIR/OUT-${i};
    mkdir -p ${outdir};

    cp ${BASE_FILE} $outdir;
    mv system.fdf $outdir;
    cp *.psf $outdir;

    #RUN
    cd $outdir;
    $INTELMPI -machinefile "$PBS_NODEFILE" -np $NUM_PROCS $SIESTA ${BASE_FILE} > ${outfile}
    
    ALL_SUP_DIR=$PBS_O_WORKDIR/${OUT_ALL}/SupCell-${i};
    mkdir -p ${ALL_SUP_DIR};
    mv system.fdf ${ALL_SUP_DIR};
    mv $outfile ${ALL_SUP_DIR};

    cd $PBS_O_WORKDIR

    rm -rf ${outdir} ;

done

