# LineWidth.x 
The computation performed by this program involves determining the phonon linewidth along the high-symmetry path within the Brillouin zone. 
This program can run in parallel through **Fortran Coarray**.  
One must install [**OpenCoarrays**](http://www.opencoarrays.org/) to use the Fortran Coarray feature with the _gfortran_ compiler. For Intel Fortran Compiler (_ifort_ or _ifx_), no extra dependencies need to be installed (recommended). 
- If you use _gfortran_ with _OpenCoarrays_, run the executable like an ordinary MPI program:
```sh
mpirun -np NPROCS ./LineWidth.x -inp [INPUT FILE] -T [TEMPERATURE] 
```
- If you use Intel Fortran Compiler, you can set the number of processors to run the program in parallel through an environment variable:
```sh
export FOR_COARRAY_NUM_IMAGES=NPROCS
```
Then run the program like any ordinary executable:
```sh
./LineWidth.x -inp [INPUT FILE] -T [TEMPERATURE] 
```

##### Executable name: _LineWidth.x_ 

#### Dependency
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [LAPACK and BLAS](https://netlib.org/lapack/lug/node11.html) or [Intel-mkl](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.mgzhln)

#### Input
The following input files are needed to run this program:
- The input namelist file. It uses four Namelist of the input Namelist file: **_FCInfo_**, **_CrystalInfo_**, **_AtmInfo_**, **_LWInfo_**, and **_HighSymPathLW_**. Detailed descriptions of the variables in the first three namelist file can be found [here](src/FC2/README.md). Variables in **_LWInfo_**, and **_HighSymPathLW_** namelist are discussed below.
- The Second-Order force-constant file (output of the program _reconst2.x_ or _renorm.x_): **FC_2nd_`${prefix}${T}`K.h5**
- The Third-Order force-constant file (output of the program _reconst3.x_): **FC_3rd_`${prefix}${T}`K.h5**

**Description of Command line arguments:**
Two variables can be optionally set through command line arguments:
- **_inp_:** The input filename is the name of the Namelist file. Default: _input.nml_
- **_T_:** Temperature T in K, at which the MD run was performed, or the temperature of the TSS method. The Force-displacement dataset file should have this temperature in the filename. Default: 300.0

**Description of the variable in _LWInfo_, and _HighSymPathLW_ Namelist:**
```sh
&LWInfo
    qmesh           = 19    19  19
    num_points      =   2
/   
    
&HighSymPathLW
    q_high_sym(:, 1) = 0.0      0.0     0.0     0
    q_high_sym(:, 2) = 0.0      0.5     0.5     100
/   

```
- **_qmesh_:** Integer variable (1d array of length 3). This variable determines the q-point mesh for the phonon linewidth calculation. 
- **_num_points_:** Integer variable. This variable specifies the number of high-symmetry q-points along which phonon linewidth is to be calculated. These points must be specified in the _HighSymPathLW_ namelist.
- **_q_high_sym_:** Real variable (4 x _num_points_ array). The first three number in a column is the coordinate of the High symmetry points in the reciprocal space in units of reciprocal lattice vectors. The fourth number determines the number of points to be taken along the path connecting the current and previous points.


#### Output
After successful execution, the following file will be generated: 
- linewidth_`${prefix}${T}`K.txt

The phonon linewidth details for the q-points along the high-symmetry path in the Brillouin zone are stored in the file _linewidth_`${prefix}${T}`K.txt_.

#### License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)
