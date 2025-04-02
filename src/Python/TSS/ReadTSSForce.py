
# ============================================================================= #
#    Copyright (C) 2022  Soham Mandal                                           #
#                                                                               #
#    This program is free software: you can redistribute it and/or modify       #
#    it under the terms of the GNU General Public License as published by       #
#    the Free Software Foundation, either version 3 of the License, or          #
#    (at your option) any later version.                                        #
#                                                                               #
#    This program is distributed in the hope that it will be useful,            #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#    GNU General Public License for more details.                               #
#                                                                               #
#    You should have received a copy of the GNU General Public License          #
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.     #
#                                                                               #
#    e-mail: phy.soham@gmail.com                                                #
# ============================================================================= #

import numpy as np
import h5py
import f90nml
import argparse
import sys


Ryau2eVAng = 25.71104309541616


class Unit:

    def __init__(self, nml_fname='input.nml') :

        print("\nReading input from file: '", nml_fname, "'\n")
        parser = f90nml.Parser()
        data_nml = parser.read(nml_fname)

        #............*****************************************************............#
        self.prefix = data_nml['FCInfo']['prefix']
        #............*****************************************************............#
        
        #............*****************************************************............#
        self.nbasis = data_nml['CrystalInfo']['Nbasis']
        self.ntyp = data_nml['CrystalInfo']['Ntyp']
        #............*****************************************************............#

        #............*****************************************************............#
        self.lattice = np.array(data_nml['CrystalInfo']['Latvec'], dtype=np.double)
        #............*****************************************************............#

        if ( len(self.lattice) != 3) :
            print("Error for input value in file: ", nml_fname) 
            print("All 3 lattice vectors are not properly provided")
            sys.exit()

        #............*****************************************************............#
        self.sup_dim = np.array(data_nml['CrystalInfo']['Sup_Cell'], dtype=np.intc)
        #............*****************************************************............#

        #............*****************************************************............#
        self.basis_frac = np.array(data_nml['AtmInfo']['atm_pos'], dtype=np.double)
        #............*****************************************************............#

        if ( len(self.basis_frac) != self.nbasis ) :
            print("Error for input value in file: ", nml_fname) 
            print("Number of basis atom position ('atm_pos') does not match with 'Nbasis'")
            sys.exit()

        #............*****************************************************............#
        self.typ_indx = np.array(data_nml['AtmInfo']['typ_indx'], dtype=np.intc)
        #............*****************************************************............#

        if ( self.ntyp != len(np.unique(self.typ_indx)) ) :
            print("Error for input value in file: ", nml_fname) 
            print("Number of unique index in 'typ_indx' does not matches with 'Ntyp'")
            sys.exit()

        #............*****************************************************............#
        self.mass_unq = np.zeros(self.ntyp, dtype=np.double)
        self.mass_unq = np.array(data_nml['AtmInfo']['mass'], dtype=np.double)
        self.atomic_no_unq = np.array(data_nml['AtmInfo']['atomic_no'], dtype=np.intc)
        self.typ_lbl_unq = data_nml['AtmInfo']['typ_lbl']
        #............*****************************************************............#

        #............*****************************************************............#
        self.mass = np.zeros(self.nbasis, dtype=np.double)
        
        for ii in range(self.nbasis) :

            atm_typ_indx = self.typ_indx[ii]
            if ( self.ntyp == 1 ) :
                self.mass[ii] = self.mass_unq
            else :
                self.mass[ii] = self.mass_unq[atm_typ_indx-1]
        #............*****************************************************............#

        #............*****************************************************............#
        self.basis_cart = np.matmul(self.basis_frac, self.lattice)
        #............*****************************************************............#


class ReadTSS :

    def __init__(self, nml_fname, T, calculator) :

        self.u = Unit(nml_fname)

        self.calculator = calculator

        self.N = np.copy(self.u.sup_dim)
        print("\tSuper cell dimension:", self.N)

        self.Natom = self.N[0]*self.N[1]*self.N[2]*len(self.u.basis_frac)

        self.recordfile = 'disp_forc_'+self.u.prefix+str(T)+'K.h5'
        self.cell_basis, self.cart_cord, self.Nstep = self.ReadRecord5()

        Nbasis = len(self.u.basis_frac)
        file0 = "./OUT/SupCell-0"+"/"+self.u.prefix+"_sup_0.out"
        print("\t*** Reading the forces at equilibrium positions. Filename: '", file0, "'***")
        self.initForc = self.InitForce()

        self.forc = np.zeros([self.Nstep, self.N[0], self.N[1], self.N[2], Nbasis, 3], dtype=np.double)

    
    def ReadRecord5(self) :

        cell_basis_record_dset = 'cell_basis'
        cart_cord_dset = 'cart_coord'
        disp_dset = 'disp_dataset'

        print("\t*** Reading record of Super-Cell information '", self.recordfile, "'***")
        h5f = h5py.File(self.recordfile, mode='r')

        cell_basis = h5f[cell_basis_record_dset][:]
        cart_cord = h5f[cart_cord_dset][:]

        Nstep = np.shape( h5f[disp_dset] )[0]

        h5f.close()

        return cell_basis, cart_cord, Nstep


    def FindLineNumber(self, lookup_word, SnapNum, filepath) :


        fl = open(filepath)

        num_lines = sum(1 for line in fl)
        read = fl.read()
        fl.seek(0)
        arr = []
        line_num = []
        for i in range(num_lines):
            arr.append(fl.readline())
        for i in range(len(arr)):
            if lookup_word in arr[i]:
                line_num.append(i+1)

        if ( (len(line_num) != 1) and (self.calculator == 'qe') ) :
            print("\t*** Error in reading file:", filepath, "***", len(line_num))

        elif ( (len(line_num) > 2) and (self.calculator == 'siesta') ) :
            print("\t*** Error in reading file:", filepath, "***", len(line_num))

        fl.close()

        return line_num[0]


    def InitForce(self) :

        SnapNum = 0

        if ( self.calculator == 'qe' ) :

            filepath = "./OUT/SupCell-"+str(SnapNum)+"/"+self.u.prefix+"_sup_"+str(SnapNum)+".out"
            lookup_word = "(cartesian axes, Ry/au):"
            line = self.FindLineNumber(lookup_word, SnapNum, filepath)

            initForc = np.loadtxt(filepath, dtype=np.double, skiprows=line+1, \
                                  usecols=(6, 7, 8), max_rows=self.Natom)

        elif ( self.calculator == 'siesta' ) :

            filepath = "./OUT/SupCell-"+str(SnapNum)+"/"+self.u.prefix+"_sup_"+str(SnapNum)+".out"
            lookup_word = "siesta: Atomic forces (eV/Ang):"
            line = self.FindLineNumber(lookup_word, SnapNum, filepath)

            initForc = np.loadtxt(filepath, dtype=np.double, skiprows=line, \
                                  usecols=(1, 2, 3), max_rows=self.Natom)

        else :

            print("\t*** Error: Only qe and siesta interface available: ", self.calulator)

        return initForc


    def ReadForce(self) :

        print()
        for SnapNum in range(1, self.Nstep+1) :

            if ( self.calculator == 'qe' ) :

                filepath = "./OUT/SupCell-"+str(SnapNum)+"/"+self.u.prefix+"_sup_"+str(SnapNum)+".out"
                lookup_word = "(cartesian axes, Ry/au):"
                line = self.FindLineNumber(lookup_word, SnapNum, filepath)

                force = np.loadtxt(filepath, dtype=np.double, skiprows=line+1, \
                                      usecols=(6, 7, 8), max_rows=self.Natom)
                force_scale = np.empty_like( force )
                force_scale = (force - self.initForc) * Ryau2eVAng

            elif ( self.calculator == 'siesta' ) :

                filepath = "./OUT/SupCell-"+str(SnapNum)+"/"+self.u.prefix+"_sup_"+str(SnapNum)+".out"
                lookup_word = "siesta: Atomic forces (eV/Ang):"
                line = self.FindLineNumber(lookup_word, SnapNum, filepath)

                force = np.loadtxt(filepath, dtype=np.double, skiprows=line, \
                                      usecols=(1, 2, 3), max_rows=self.Natom)
                force_scale = np.empty_like( force )
                force_scale = (force - self.initForc)

            else :

                print("\t*** Error: Only qe and siesta interface available: ", self.calulator)


            for atm in range(self.Natom) :
                nx, ny, nz, mu = self.cell_basis[atm, :]
                self.forc[SnapNum-1, nx, ny, nz, mu, :] = np.copy(force_scale[atm, :])

            print("\tRead Force: No. of TSS covered:", SnapNum, "/", self.Nstep, \
                    "/|\ Filename: '", filepath, "'", end='\r', flush=True)


    def Writeh5(self, T) :

        filename = 'disp_forc_'+self.u.prefix+str(T)+'K.h5'

        print()
        print("\t*** Writing force-dataset in file:'", filename, "'***")

        h5f = h5py.File(filename, mode='a')

        force_dset = 'force_dataset'

        if force_dset in h5f :
            del h5f[force_dset]

        h5f.create_dataset(force_dset, data=self.forc)
        
        h5f.close()

    
    #.........::::::::::::::: For Debug purpose :::::::::::::::............#
    def Readh5(self, T) :

        filename = 'disp_forc_'+self.u.prefix+str(T)+'K.h5'
        print("Debug: Reading file: {0}".format(filename))

        disp_dset = 'disp_dataset'
        force_dset = 'force_dataset'

        h5f = h5py.File(filename, mode='r')

        frc = h5f[force_dset][:]
        cord = h5f[disp_dset][:]

        h5f.close()

        return frc, cord
    #.........::::::::::::::: For Debug purpose :::::::::::::::............#

                    
if __name__ == "__main__" :

    parser = argparse.ArgumentParser()

    # Run: python3 ReadTSSForce.py -inp input.nml -T 80 
    parser.add_argument("-inp", type=str, default='input.nml', help="Enter the name of input file")
    parser.add_argument("-T", type=float, default=300, help="Temperature in Kelvin")
    parser.add_argument("-Fcalculator", type=str, default='siesta', help="Force calculator (qe or siesta)")

    args = parser.parse_args()
    nml_fname = args.inp
    T = args.T
    calculator = args.Fcalculator

    data = ReadTSS(nml_fname, T, calculator)

    data.ReadForce()
    data.Writeh5(T)

    print()

    #.........::::::::::::::: For Debug purpose :::::::::::::::............#
    debug = False
    if ( debug ) :

        frc, cord = data.Readh5(T)
        num = 52

        print("Displacement/Force data (Snap-shot num = ", num, "): ") 
        num_atm = 0
        for nx in range(data.N[0]) :
            for ny in range(data.N[1]) :
                for nz in range(data.N[0]) :
                    for bb in range(len(data.u.basis_frac)) :

                        #print(frc[(num-1), nx, ny, nz, bb, :])

                        print(cord[num-1, nx, ny, nz, bb, :] + \
                              data.cart_cord[num_atm, :])

                        num_atm += 1
        print()

        #print('cord[16, 0, 1, 2, 0, 0]:', cord[16, 0, 1, 2, 0, 0])
        #print('frc[16, 0, 1, 2, 0, 0]:', frc[16, 0, 1, 2, 0, 0])
    #.........::::::::::::::: For Debug purpose :::::::::::::::............#


