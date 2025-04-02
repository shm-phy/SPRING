# FreeEng.x 
This program calculates Helmholtz's free energy considering anharmonicity up to fourth order. 
This program can run in parallel through **Fortran Coarray**.  
One must install [**OpenCoarrays**](http://www.opencoarrays.org/) to use the Fortran Coarray feature with the _gfortran_ compiler. For Intel Fortran Compiler (_ifort_ or _ifx_), no extra dependencies need to be installed (recommended). 
- If you use _gfortran_ with _OpenCoarrays_, run the executable like an ordinary MPI program:
```sh
mpirun -np NPROCS ./FreeEng.x -inp [INPUT FILE] -T [TEMPERATURE] -AnHar [T/F]
```
- If you use Intel Fortran Compiler, you can set the number of processors to run the program in parallel through an environment variable:
```sh
export FOR_COARRAY_NUM_IMAGES=NPROCS
```
Then run the program like any ordinary executable:
```sh
./FreeEng.x -inp [INPUT FILE] -T [TEMPERATURE] -AnHar [T/F]
```

##### Executable name: _LineWidth.x_ 

#### Dependency
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [LAPACK and BLAS](https://netlib.org/lapack/lug/node11.html) or [Intel-mkl](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.mgzhln)

#### Input
The following input files are needed to run this program:
- The input namelist file. It uses four Namelist of the input Namelist file: **_FCInfo_**, **_CrystalInfo_**, and **_AtmInfo_** . Detailed descriptions of the variables in these namelist can be found [here](src/FC2/README.md). 
- The Second-Order force-constant file (output of the program _reconst2.x_ or _renorm.x_): **FC_2nd_`${prefix}${T}`K.h5**
- **(Optional)** The Third-Order force-constant file (output of the program _reconst3.x_): **FC_3rd_`${prefix}${T}`K.h5**
- **(Optional)** The Fourth-Order force-constant file (output of the program _reconst4.x_): **FC_4th_`${prefix}${T}`K.h5**

**Description of Command line arguments:**
Three variables can be optionally set through command line arguments:
- **_inp_:** The input filename is the name of the Namelist file. Default: _input.nml_
- **_T_:** Temperature T in K, at which the MD run was performed, or the temperature of the TSS method. The Force-displacement dataset file should have this temperature in the filename. Default: 300.0
- **_AnHar_:** Logical variable. This variable determines whether to compute the anharmonic comtribution (3rd and 4th order) of free energy.

#### License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)
