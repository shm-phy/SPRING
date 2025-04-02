# FC4.x
Find the independent Fourth-Order IFCs, combining point-group symmetry, permutation symmetry and accoustic sum rule (ASR) of IFCs.

##### Executable name: _FC4.x_

#### Dependency
- [Spglib](https://spglib.github.io/spglib/)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [LAPACK and BLAS](https://netlib.org/lapack/lug/node11.html) or [Intel-mkl](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.mgzhln)

#### Input
It uses three Namelist of the input Namelist file: **_FCInfo_**, **_CrystalInfo_**, and **_AtmInfo_**.
In the Namelist **_FCInfo_**, variables used by the executable _FC4.x_ are **_prefix, FourthFC_**, and **_nn_4th_**. 
_FC4.x_ requires all the variables (**_Nbasis_**, **_Ntyp_**, **_Latvec_**, **_Sup_Cell_**) of the **_CrystalInfo_** Namelist.
From the **_AtmInfo_** Namelist, it requires **_atm_pos, typ_indx_**, and **_typ_lbl_** variables as input.
Other Namelists may be present in the input file, but _FC4.x_ do not read them.

**Example Input Namelist file:**
```sh
&FCInfo
    prefix          =   'Si'
    SecondFC        =   .true.
    nn_2nd          =   11
    ThirdFC         =   .true.
    nn_3rd          =   6
    FourthFC        =   .true.
    nn_4th          =   2
    Force_dd        =   .false.
/
&CrystalInfo

    Nbasis          = 2
    Ntyp            = 1
    Latvec(:, 1)    = 3.840900274  0.0         0.0
    Latvec(:, 2)    = 1.920450137  3.326317210 0.0
    Latvec(:, 3)    = 1.920450137  1.108772404 3.136081944
    Sup_Cell        = 5     5     5

/
&AtmInfo

    atm_pos(:, 1)   = 0.50    0.50    0.50
    atm_pos(:, 2)   = 0.75    0.75    0.75
    typ_indx        =   1       1

    typ_lbl         =   'Si'
    atomic_no       =   14
    mass            = 28.0855

    dielec(:, 1)    = 14.581827722559         0.000000000000          0.000000000000
    dielec(:, 2)    =  0.000000000000         14.581827728166        -0.000000002887
    dielec(:, 3)    =  0.000000000000        -0.000000002887         14.581829202410

    BornZ(:, 1, 1)  =  0.0000000  0.0000000   0.0000000
    BornZ(:, 2, 1)  =  0.0000000  0.0000000  -0.0000000
    BornZ(:, 3, 1)  =  0.0000000 -0.0000000   0.0000000

    BornZ(:, 1, 2)  = -0.0000000   0.0000000    0.0000000
    BornZ(:, 2, 2)  =  0.0000000  -0.0000000    0.0000000
    BornZ(:, 3, 2)  =  0.0000000   0.0000000   -0.0000000

/
```
**Description:**
- **prefix:** Any string of characters. Generally, the name of the Material is a useful one.   This prefix is used in the output filename generated after the successful execution.  
- **SecondFC:** Logical variable. If `.false.`, the program will stop without any calculation for 2nd order IFCs.  
- **nn_2nd:** Integer variable. Nearest neighbor cut-off of 2nd order IFCs. Upto 11 or 12th nearest neighbor is a good choice for 2nd Order IFCs of Silicon. Do not increase this value arbitrarily, especially beyond the L/2 distance of the supercell. 
- **ThirdFC:** Same as _SecondFC_, but for the calculation of 3rd Order IFCs.
- **nn_3rd:** Same as _nn_2nd_, but for the 3rd Order IFCs.
- **FourthFC:** Same as _SecondFC_, but for the calculation of 4th Order IFCs.
- **nn_4th:** Same as _nn_2nd_, but for the 4th Order IFCs.
- **Force_dd:** Logical variable. Determines whether the Ewald-sum technique will be used for the long-range dipole-dipole interaction part of the 2nd Order IFCs. For polar materials (e.g., NaCl), it should be true to get the LO-TO splitting in the phonon dispersion. If it is `.true.`, one should provide the dielectric tensor (_dielec_) and Born-effective charge (_BornZ_) of the Material in the _AtmInfo_ Namelist.

- **Nbasis:** Integer variable. The number of basis atoms in a unit cell of the crystalline Material.
- **Ntyp:** Integer variable. The number of distinct elements present in the basis set. For Silicon, _Nbasis_ would be 2 and _Ntyp_ 1. For NaCl, both are 2.
- **Latvec:** Real variable (3x3 array). Three lattice vectors in Anstrom. _Latvec(:, 1)_ is the column containing one of the lattice vectors (Lx). _Latvec(:, 2)_ and _Latvec(:, 3)_ are other lattice vectors. 
- **Sup_Cell:** Integer variable (1d array of length 3). The extent of the supercell along three axes that will be formed from the unit cell. **Always provide an odd number**. 

- **atm_pos:** Real variable (3 x _Nbasis_ array). Positions of basis atoms in the unit of _Latvec_. For example, _atm_pos(:, 1)_ is a column containing the position of basis atom 1 in _Latvec_ unit. Similarly, for all other basis atoms up to _Nbasis_.
- **typ_indx:** Integer variable (1d array of length _Nbasis_). The type of a basis atom is an integer, starting from 1. It is a unique integer for each separate element in the basis set. For example, in the case of Silicon, there are 2 basis atoms per unit cell, each of which is the same element Si. So, _typ_indx_ = 1   1. In the case of NaCl _typ_indx_ = 1   2. 
- **typ_lbl:** Array of string (1d array of length _Ntyp_). For Silicon _typ_lbl_ = 'Si'. In the case of NaCl _typ_lbl_ = 'Na'  'Cl'.
- **atomic_no:** Integer variable (1d array of length _Ntyp_). Atomic-number of each unique element in the basis set.
- **mass:** Real variable (1d array of length _Ntyp_). The atomic mass of each unique element in the basis set in units of gm/mol.
- **dielec:** Real variable (3x3 array). Dielectric tensor of the Material. _dielec(:, 1)_ is the column for (e~xx~ e~xy~ e~xz~)
- **BornZ:** Real variable (3 x 3 x _Nbasis_ array). Born effective charge tensor for each atom in the basis set. _BornZ(:, :, 1)_ is the 3x3 array of Born effective charge tensor for basis atom 1. 

#### Command line arguments
Four variables can be optionally set through command line arguments. 
```sh
./FC4.x -inp [input filename] -ChnkSz [ChunkSize in MB] -AdvncOut [T/F] -eps [REAL]
```
Input filename is the name of the Namelist file. Default: _input.nml_
_ChnkSz_ is the Chunk Size of Inter-FC Matrix that will be written in file: InterMatChunk_4th_`${prefix}`.h5. Default: 512 MB. Reducing this may speed up the execution time. 
_AdvncOut_ is a logical variable that determines whether advancing output is required. Default: `.false.`
_eps_ is a Real variable for the point group symmetry precision. Default:  1.0E-5. 

#### Output
After successful execution, two files will be generated: 
- FC_4th_common_`${prefix}`_F.h5
- InterMatChunk_4th_`${prefix}`.h5

InterMatChunk_4th_`${prefix}`.h5 can be deleted afterward. FC_4th_common_`${prefix}`_F.h5 is required for the next steps. 

#### License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)
