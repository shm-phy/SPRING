[![Spring logo](tools/spring_logo2_bg.png)]()

[![doi](https://img.shields.io/badge/DOI-0.0.0-blue?logo=arxiv)](https://doi.org/10.48550/arXiv.2402.02787) [![doi](https://img.shields.io/badge/License-GPLv3-brightgreen?logo=gnu)](https://www.gnu.org/licenses/gpl-3.0.en.html) [![OS](https://img.shields.io/badge/Operating%20System-Linux-blue?logo=linux)]() [![version](https://img.shields.io/badge/Version-0.1.0-yellow?logo=github)]() 

SPRING stands for Self-consistent Phonon Renormalizing-package, representing an extensive numerical framework designed for evaluating the thermal and thermodynamic characteristics of crystalline insulators and semiconducting materials.
It is implemented according to the Fortran-2018 standard, with a substantial portion parallelized using Fortran Coarray and OpenMP.
For more comprehensive information, refer to our pre-print publication: https://arxiv.org/abs/2402.02787


## Features
- Calculation of interatomic force constants (IFCs) up to the fourth order using the temperature-dependent effective potential (TDEP) technique, incorporating permutation symmetry, the acoustic sum rule, and space group symmetry.
- Implementation of the Thermal Stochastic Snapshot (TSS) method for mitigating the computational cost of ab initio molecular dynamics.
- Incorporation of the precise Ewald summation method to handle long-range 2nd order IFCs in polar materials. The exact Ewald sum is applied in both the LO-TO splitting of phonon dispersion and the calculation of phonon group velocities.
- Self-consistent renormalization of harmonic IFCs, with the extension of this renormalization effect to anharmonic 3rd and 4th-order IFCs.
- Phonon linewidth calculation along the high-symmetry path within the Brillouin Zone. 
- Calculation of thermal conductivity involving three-phonon, four-phonon, and isotopic impurity scattering mechanisms.
- Computation of Helmholtz Free Energy considering anharmonicity up to fourth order.

## Dependency
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/): The Hierarchical Data Format version 5 (HDF5) is required as a dependency. HDF5 is a data model, library, and file format for storing and managing large and complex data. It provides a versatile and efficient way to store and organize various types of data, making it essential for handling data in the program.
- [LAPACK and BLAS](https://netlib.org/lapack/lug/node11.html) or [Intel-mkl](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.mgzhln) or [OpenBlas](https://www.openblas.net/): Linear Algebra PACKage (LAPACK) and Basic Linear Algebra Subprograms (BLAS) are libraries that provide a set of routines for solving linear equations, eigenvalue problems, singular value decomposition, and more. They are crucial for performing various mathematical computations required by the program. You can choose to use LAPACK and BLAS from various sources such as the official LAPACK/BLAS libraries, Intel Math Kernel Library (Intel MKL), or OpenBlas.
- [Spglib](https://spglib.github.io/spglib/): Spglib is a library for handling crystal symmetry in various crystallographic operations. It can generate symmetry operations, find primitive cells, obtain space group information, and more. Spglib is necessary for managing crystallographic aspects and symmetry operations within the program.
- Fortran Coarray: Fortran Coarray is utilized for the parallelization of various subprograms. One must install [**OpenCoarrays**](http://www.opencoarrays.org/) to use the Fortran Coarray feature with the _gfortran_ compiler. For Intel Fortran Compiler (_ifort_ or _ifx_), no extra dependencies need to be installed (recommended).

You have the flexibility to install these dependencies using your preferred method. Within the 'tools' directory, there are two bash scripts available—one for [GNU compilers](tools/Install_gnu.sh) and another for [Intel compilers](tools/Install_intel.sh). These scripts are designed to fetch and install all the required dependencies automatically. It's important to note that before running these scripts successfully, you must have various developer tools preinstalled, including compilers like GCC or Intel OneAPI, GNU make, wget, and others. 
## Installation
The building procedure follows the conventional `make` file-based approach. Once you've successfully installed all the necessary dependencies, proceed to set the library paths within the `Makefile.config` file. In this file, you also have the option to specify your preferred compiler vendor—either GNU or Intel. It's worth noting that, at present, only `gfortran` (GCC) and `ifort` (Intel OneAPI) compilers have been tested and are supported.
For additional flexibility in configuring various installation settings, you can modify the variables found in the `Makefile.config` file. This file has been extensively documented to ensure straightforward usability.
Once all the settings in `Makefile.config` file are set, run the following command to compile and install: 
```sh
make all
```
## Documentation
Comprehensive documentation for each subprogram is located within the subdirectories under the `src` directory.
## Workflow Diagram
[![Workflowdiag](tools/workflow_diagram.pdf)]()

## License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)
