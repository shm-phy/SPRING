# TSS.x 
This program generates snapshots of displaced atoms of the supercell for the Thermal Stochastic Snapshot (TSS) method. 

##### Executable name: _TSS.x_

#### Dependency
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [LAPACK and BLAS](https://netlib.org/lapack/lug/node11.html) or [Intel-mkl](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.mgzhln)

#### Input
Two input files are needed to run this program:
- The input namelist file. It uses four Namelist of the input Namelist file: **_FCInfo_**, **_CrystalInfo_**, **_AtmInfo_**, and **_EwaldInfo_**. Detailed descriptions of the variables in the namelist **_FCInfo_**, **_CrystalInfo_**, and **_AtmInfo_** can be found [here](src/FC2/README.md). Variables in the **_EwaldInfo_** namelist are discussed [here](src/Long_DD/README.md).
- The Second-Order force-constant file (output of the program _reconst2.x_): **FC_2nd_`${prefix}${T}`K.h5**

**Description of Command line arguments:**
```sh
./TSS.x -inp [INPUT FILE] -T [TEMPERATURE] -Tsnap [TEMPERATURE] -Nsnap [INTEGER]
```
Four variables can be optionally set through command line arguments:
- **_inp_:** The input filename is the name of the Namelist file. Default: _input.nml_
- **_T_:** Temperature T in K, at which the MD run was performed, or the temperature of the TSS method. The Second-Order IFC file should have this temperature in the filename. Default: 300.0
- **_Tsnap_:** Real variable. _Tsnap_ is the temperature at which snapshots will be created. It can be different from _T_, the temperature of the Second-Order IFCs. For example, you can use Second-Order IFCs at 300K to generate snapshots for 320K. Then in a self-consistent manner, you can approach that particular temperature. Default: 300.0
- **_Nsnap_:** Integer variable. The number of snapshots that will be created. Default: 50

#### Output
- After successful execution, all the snapshots will be created inside the directory: _SnapShots/_. The files are in Quantum-Espresso input format.
- Another file containing the displacement dataset will also be created: **disp_forc_`${prefix}${Tsnap}`K.h5**

#### License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)
