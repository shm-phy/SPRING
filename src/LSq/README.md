# LSq.x 
This program finds the independent IFCs ( 2nd-Order, 3rd-Order(optionally), and 4th-Order(optionally) ) from a least-square fitting of the matrix equation: _D k = f_. Here _D_ is the displacement matrix ( 2nd-Order, 3rd_Order(optionally), and 4th-Order(optionally) combined ) found from the program create_disp_mat_`${Order}`.x. _f_ is a row of forces on atoms in the super-cell.

##### Executable name: _LSq.x_

#### Dependency
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [LAPACK and BLAS](https://netlib.org/lapack/lug/node11.html) or [Intel-mkl](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.mgzhln)

#### Input
Input files are needed to run this program:
- The input namelist file. It uses three Namelist of the input Namelist file: **_FCInfo_**, **_CrystalInfo_**, and **_AtmInfo_**. Detailed descriptions of the variables in the namelist file can be found [here](src/FC2/README.md).
- The displacement-matrix file generated from the program:  create_disp_mat_`${Order}`.x. **FD_2nd_`${prefix}${T}`K.h5** must present. To find  3rd-Order and 4th-Order IFCs (optionally, if _ThirdFC_ and _FourthFC_ are `.true.` in _FCInfo_ namelist), the respective displacement-matrix files are needed (**FD_3rd_`${prefix}${T}`K.h5** and **FD_4th_`${prefix}${T}`K.h5**)
- To separate the long-range dipole-dipole interaction part of the force (found from the Ewald-Sum technique), the following file must also be present: **FC2_dd_`${prefix}${T}`K.h5**. For this separation approach in the case of polar materials, _Force_dd_ must be `.true.` in _FCInfo_ namelist.

 
**Description of Command line arguments:**
```sh
./LSq.x -inp [INPUT FILE] -T [TEMPERATURE] -Renorm [T/F] -CutOff [REAL]
```
Four variables can be optionally set through command line arguments:
- **_inp_:** The input filename is the name of the Namelist file. Default: _input.nml_
- **_T_:** Temperature T in K, at which the MD run was performed, or the temperature of the TSS method. The displacement-matrix file should have this temperature in the filename. Default: 300.0
- **_Renorm_:** Logical variable. This variable determines whether to calculate renormalized 3rd-Oder and 4th-Order IFCs through the least-square method. Default: `.false.`
- **_CutOff_:** Real variable. Singular value cutoff for the least-square routine. Default: 1.0E-5.
#### Output
After successful execution, output will be written in the same files: 
- FD_2nd_`${prefix}${T}`K.h5 (optionally FD_3rd_`${prefix}${T}`K.h5 and FD_4th_`${prefix}${T}`K.h5 if 3rd-Order and 4th-Order are needed)

#### License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)
