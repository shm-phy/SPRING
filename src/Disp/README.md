# disp.x 
This program calculates the phonon dispersion along the high-symmetry path in the reciprocal space.

##### Executable name: _disp.x_

#### Dependency
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [LAPACK and BLAS](https://netlib.org/lapack/lug/node11.html) or [Intel-mkl](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.mgzhln)

#### Input
Two input files are needed to run this program:
- The input namelist file. It uses six Namelist of the input Namelist file: **_FCInfo_**, **_CrystalInfo_**, **_AtmInfo_**, **_EwaldInfo_**, **_PhononDispInfo_**, and **_HighSymPath_**. Detailed descriptions of the variables in the namelist **_FCInfo_**, **_CrystalInfo_**, and **_AtmInfo_** can be found [here](src/FC2/README.md). Variables in the **_EwaldInfo_** namelist are discussed [here](src/Long_DD/README.md). Variables in **_PhononDispInfo_**, and **_HighSymPath_** namelist are discussed below.
- The Second-Order force-constant file (output of the program _reconst2.x_): **FC_2nd_`${prefix}${T}`K.h5**

**Description of the variable in _PhononDispInfo_ and _HighSymPath_ Namelist:**
```sh
&PhononDispInfo
    PhononDisp      =   .true.
    LongEW          =   .false.
    num_points      =   6
/
&HighSymPath
    q_high_sym(:, 1) = 0.0      0.0     0.0     0
    q_high_sym(:, 2) = 0.0      0.5     0.5     500
    q_high_sym(:, 3) = 0.25     0.625   0.625   500
    q_high_sym(:, 4) = 0.375     0.75   0.375   0
    q_high_sym(:, 5) = 0.0      0.0     0.0     500
    q_high_sym(:, 6) = 0.5      0.5     0.5     500
/
```

- **_PhononDisp_:** Logical variable. This variable determines whether to calculate phonon dispersion along the high-symmetry path. If it is `.true.`, then _HighSymPath_ namelist must be present.
- **_LongEW_:** Logical variable. This variable determines whether the long-range correction for the Second-Order IFCs is to be calculated through an Ewald-Summation technique. It should be `.true.` for polar materials. If it is `.true.`, then _EwaldInfo_ namelist must be present.
- **_num_points_:** Integer variable. This variable specifies the number of high-symmetry q-points along which phonon dispersion is to be calculated. These points must be specified in the _HighSymPath_ namelist.

- **_q_high_sym_:** Real variable (4 x _num_points_ array). The first three number in a column is the coordinate of the High symmetry points in the reciprocal space in units of reciprocal lattice vectors. The fourth number determines the number of points to be taken along the path connecting the current and previous points. 

**Description of Command line arguments:**
```sh
./disp.x -inp [INPUT FILE] -T [TEMPERATURE]
```
Two variables can be optionally set through command line arguments:
- **_inp_:** The input filename is the name of the Namelist file. Default: _input.nml_
- **_T_:** Temperature T in K, at which the MD run was performed, or the temperature of the TSS method. The Force-displacement dataset file should have this temperature in the filename. Default: 300.0

#### Output
If _PhononDisp_ is `.true.`, then the following file containing phonon frequencies along the high-symmetry path will be created:
- **dispersion_`${prefix}${T}`K.h5**

#### License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)
