# longEw.x 
This program calculates the long-range dipole-dipole interaction force through an Ewald-sum technique. This long-range force is subtracted from the actual MD force (or the force of the TSS method) to fit only the short-range part of the force constant in the least-square method. The long-range portion of the IFCs (in reciprocal space) will be later added to get the Dynamical Matrix. 

##### Executable name: _longEw.x_

#### Dependency
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [LAPACK and BLAS](https://netlib.org/lapack/lug/node11.html) or [Intel-mkl](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.mgzhln)

#### Input
Two input files are needed to run this program:
- The input namelist file. It uses three Namelist of the input Namelist file: **_FCInfo_**, **_CrystalInfo_**, **_AtmInfo_**, and **_EwaldInfo_**. Detailed descriptions of the variables in the namelist **_FCInfo_**, **_CrystalInfo_**, and **_AtmInfo_** can be found [here](src/FC2/README.md). Variables in the **_EwaldInfo_** namelist are discussed below. 
- The Force-displacement dataset file generated from the Molecular dynamics snapshots or from the Thermal Stochastic Snapshot (TSS) method. Filename format: **disp_forc_`${prefix}${T}`K.h5**

**Description of the variable in _EwaldInfo_ Namelist:**
```sh
&EwaldInfo
    Rmesh   =   5   5   5
    Gmesh   =   5   5   5
    Lmb     =   1.01
    decide_EwParam  =   .true.
    prec_Ew         =   1.0E-7

/
```
- **_Rmesh_:** Integer variable (1d array of length 3). This variable determines the Real-Space cutoff in the Ewald Summation technique.
- **_Gmesh_:** Integer variable (1d array of length 3). This variable determines the Reciprocal-Space cutoff in the Ewald Summation technique.
- **_Lmb_:** Real variable. _Lmb_ is a parameter of the Ewald Summation in the unit of Anstrom.
- **_decide_EwParam_:** Logical variable. This variable decides whether the parameters (_Rmesh_, _Gmesh_ and _Lmb_) will be predicted internally or use the user-given parameters. 
- **_prec_Ew_:** Real variable. This variable is applicable only when _decide_EwParam_ is `.true.`. It is used for the internal prediction of the Ewald Summation parameters.

**Description of Command line arguments:**
Two variables can be optionally set through command line arguments:
- **_inp_:** The input filename is the name of the Namelist file. Default: _input.nml_
- **_T_:** Temperature T in K, at which the MD run was performed, or the temperature of the TSS method. The Force-displacement dataset file should have this temperature in the filename. Default: 300.0

#### Output
After successful execution, the following file will be created: 
- **FC2_dd_`${prefix}${T}`K.h5**

#### License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)
