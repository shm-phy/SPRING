# create_disp_mat2.x 
This program creates the displacement matrix and the force row to find the independent Second Order IFCs through a least square method. 
This program can run in parallel through **Fortran Coarray**. _An OpenMP version of this program is also available for shared memory multiprocessing._ 
One must install [**OpenCoarrays**](http://www.opencoarrays.org/) to use the Fortran Coarray feature with the _gfortran_ compiler. For Intel Fortran Compiler (_ifort_ or _ifx_), no extra dependencies need to be installed (recommended). 
- If you use _gfortran_ with _OpenCoarrays_, run the executable like an ordinary MPI program:
```sh
mpirun -np NPROCS ./create_disp_mat2.x -inp [INPUT FILE] -T [TEMPERATURE] -mu [INTEGER] -alpha [INTEGER] -DynSchd [LOGICAL] -Nchunk [INTEGER] -IdleImg1 [LOGICAL]
```
- If you use Intel Fortran Compiler, you can set the number of processors to run the program in parallel through an environment variable:
```sh
export FOR_COARRAY_NUM_IMAGES=NPROCS
```
Then run the program like any ordinary executable:
```sh
./create_disp_mat2.x -inp [INPUT FILE] -T [TEMPERATURE] -mu [INTEGER] -alpha [INTEGER] -DynSchd [LOGICAL] -Nchunk [INTEGER] -IdleImg1 [LOGICAL]
```

##### Executable name: _create_disp_mat2.x_

#### Dependency
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [LAPACK and BLAS](https://netlib.org/lapack/lug/node11.html) or [Intel-mkl](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.mgzhln)

#### Input
Three input files are needed to run this program:
- The input namelist file. It uses three Namelist of the input Namelist file: **_FCInfo_**, **_CrystalInfo_**, and **_AtmInfo_**. Detailed descriptions of the variables in the namelist file can be found [here](src/FC2/README.md).
- The Force-displacement dataset file generated from the Molecular dynamics snapshots or from the Thermal Stochastic Snapshot (TSS) method. Filename format: disp_forc_`${prefix}${T}`K.h5
- The file containing information about the independent 2nd-Order IFCs. This file is generated from the program _FC2.x_. Filename format: FC_2nd_common_`${prefix}`_F.h5

 
**Description of Command line arguments:**
Seven variables can be optionally set through command line arguments:
- **_inp_:** The input filename is the name of the Namelist file. Default: _input.nml_
- **_T_:** Temperature T in K, at which the MD run was performed, or the temperature of the TSS method. The Force-displacement dataset file should have this temperature in the filename. Default: 300.0
- **_mu and alpha_:** Both are integer variables. _mu_ takes value 1 to _Nbasis_. _alpha_ always runs from 1 to 3 (Cartesian index). One can run the following simple script to seamlessly cover all the basis index (_mu_) and cartesian index (_alpha_):
```sh
#!/bin/bash

VENDOR=intel #gnu
INFILE='input.nml'
TEMP=80
NBASIS=2
NPROCS=6
chunk=81

export FOR_COARRAY_NUM_IMAGES=6

for ((bb=1; bb <= $NBASIS; bb++))
do
    for ((a=1; a <= 3; a++))
    do
        if [ ${VENDOR} = "gnu" ]; then
            mpirun -np $NPROCS ./create_disp_mat2.x -inp $INFILE -T $TEMP -mu $bb -alpha $a -DynSchd T -Nchunk $chunk -IdleImg1 T
        elif [ ${VENDOR} = "intel" ]; then
            ./create_disp_mat2.x -inp $INFILE -T $TEMP -mu $bb -alpha $a -DynSchd T -Nchunk $chunk -IdleImg1 T
        fi
    done
done
```
Between the first (1, 1) and last (_Nbasis_, 3) _mu_-_alpha_ set, all other (_mu_, _alpha_) can be run in any order or concurrently. 
- **_DynSchd_:** Logical variable. Determines whether Dynamic Scheduling or Static Scheduling is used for parallel processing. Sometimes _DynSchd_ `.true.` leads to faster execution. Default: `.false.`
- **_Nchunk_:** Integer variable. This variable is applicable only when _DynSchd_ is `.true.`. _Nchunk_ is the number of chunks the IFCs will be divided. It must be greater than _NPROCS_. Default: 60.
- **_IdleImg1_:** Logical variable. This variable is applicable only when _DynSchd_ is `.true.`. This variable determines whether Image 1 will be idle and coordinate the work of other images. Setting this variable to `.true.` sometimes leads to faster execution. Default: `.true.` 

#### Output
After successful execution, two files will be generated: 
- FDpart_2nd_`${prefix}${T}`K.h5
- FD_2nd_`${prefix}${T}`K.h5

FDpart_2nd_`${prefix}${T}`K.h5 can be deleted afterward. FD_2nd_`${prefix}${T}`K.h5 is required for the next steps. 

### OpenMP version: 
Set the number of threads through the environment variable `OMP_NUM_THREADS`. 
Only two variables needed to be set through command line arguments:
```sh
./create_disp_mat2.x -inp $INFILE -T $TEMP
```

#### License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)
