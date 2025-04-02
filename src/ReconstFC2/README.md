# reconst2.x 
This program reconstructs the whole set of Second-Order IFCs from the linear combination of independent IFCs. 
This program can run in parallel through **Fortran Coarray**. _An OpenMP version of this program is also available for shared memory multiprocessing._ 
One must install [**OpenCoarrays**](http://www.opencoarrays.org/) to use the Fortran Coarray feature with the _gfortran_ compiler. For Intel Fortran Compiler (_ifort_ or _ifx_), no extra dependencies need to be installed (recommended). 
- If you use _gfortran_ with _OpenCoarrays_, run the executable like an ordinary MPI program:
```sh
mpirun -np NPROCS ./reconst2.x -inp [INPUT FILE] -T [TEMPERATURE] 
```
- If you use Intel Fortran Compiler, you can set the number of processors to run the program in parallel through an environment variable:
```sh
export FOR_COARRAY_NUM_IMAGES=NPROCS
```
Then run the program like any ordinary executable:
```sh
./reconst2.x -inp [INPUT FILE] -T [TEMPERATURE] 
```

##### Executable name: _reconst2.x_

#### Dependency
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [LAPACK and BLAS](https://netlib.org/lapack/lug/node11.html) or [Intel-mkl](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.mgzhln)

#### Input
Three input files are needed to run this program:
- The input namelist file. It uses three Namelist of the input Namelist file: **_FCInfo_**, **_CrystalInfo_**, and **_AtmInfo_**. Detailed descriptions of the variables in the namelist file can be found [here](src/FC2/README.md).
- The file containing information about the independent 2nd-Order IFCs. This file is generated from the program _FC2.x_. Filename format: **FC_2nd_common_`${prefix}`_F.h5**
- The displacement-matrix file after the execution of the program _LSq.x_: **FD_2nd_`${prefix}${T}`K.h5**

 
**Description of Command line arguments:**
Two variables can be optionally set through command line arguments:
- **_inp_:** The input filename is the name of the Namelist file. Default: _input.nml_
- **_T_:** Temperature T in K, at which the MD run was performed, or the temperature of the TSS method. The Force-displacement dataset file should have this temperature in the filename. Default: 300.0

#### Output
After successful execution, the following file of Second-Order IFCs (only short-range part) of the material will be created: 
- **FC_2nd_`${prefix}${T}`K.h5**

### OpenMP version: 
Set the number of threads through the environment variable `OMP_NUM_THREADS`. 
Only two variables needed to be set through command line arguments:
```sh
./reconst2.x -inp $INFILE -T $TEMP
```

#### License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)

