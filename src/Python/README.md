# Python
This directory contains two python scripts:
_create_SupCell.py_: To create a supercell from the unit-cell of the crystalline solid.
_ReadForcDisp.py_: To read the output file of the molecular dynamics run (in _Siesta_) and create the force-displacement dataset. 

#### Dependency
- [NumPy](https://numpy.org/)
- [h5py](https://www.h5py.org/)
- [f90nml](https://f90nml.readthedocs.io/en/latest/)

#### Input

- Both the python script needs the input namelist file. . It uses three Namelist of the input Namelist file: **_FCInfo_**, **_CrystalInfo_**, and **_AtmInfo_**. Detailed descriptions of the variables in the namelist **_FCInfo_**, **_CrystalInfo_**, and **_AtmInfo_** can be found [here](src/FC2/README.md). 

ReadForcDisp.py needs the following files to run:
- The supercell information file that was created by _create_SupCell.py_. Filename format: **record_${prefix}.h5**
- The output log file of the MD run in _Siesta_. This file contains the displacements and force for each timestep of the MD run. 
- The MDE file created after the MD run in _Siesta_. This file contains thermodynamic information of the system at each timestep. 

**Description of Command line arguments:**
_create_SupCell.py_ can be run as: 
```sh
python3 create_SupCell.py -inp [INPUT FILE]
```
_ReadForcDisp.py_ can be run as: 
```sh
python3 ReadForcDisp.py -inp [INPUT FILE] -T [TEMPERATURE] -ts 5000 -DFfile out.log -MDEfile NaCl.MDE
```

- **_inp_:** The input filename is the name of the Namelist file. Default: _input.nml_
- **_T_:** Temperature T in K, at which the MD run was performed, or the temperature of the TSS method. Default: 300.0
- **_ts_:** Integer variable. This variable determines the total number of timesteps in the MD run. Default: 5000
- **_DFfile_:** The output log filename of the MD run in Siesta. This file contains the displacements and force for each timestep of the MD run. Default: _out.log_
- **_MDEfile_:** The MDE filename created after the MD run in Siesta. This file contains thermodynamic information of the system at each timestep.

#### Output
After successful execution, create_SupCell.py will create two files:
- **system_${prefix}.fdf:** This file contains the coordinate of the atoms in the supercell. It is in _Siesta_ input file format.
- **record_${prefix}.h5:** This file contains the supercell information which will be used by _ReadForcDisp.py_.

After successful execution, ReadForcDisp.py will create the following file containing the force-displacement dataset: 
- **disp_forc_`${prefix}${T}`K.h5**

#### License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)
