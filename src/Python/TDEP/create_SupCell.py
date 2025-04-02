
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
import argparse
import sys

import f90nml

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


class MakeSupCell :

    def __init__(self, fname) :

        self.u = Unit(fname)

        self.Nx = self.u.sup_dim[0]
        self.Ny = self.u.sup_dim[1]
        self.Nz = self.u.sup_dim[2]

        self.frac_coord = np.zeros([self.Nx*self.Ny*self.Nz*len(self.u.basis_frac), 3], dtype=np.double)
        self.cart_coord = np.zeros_like(self.frac_coord)
        self.cell_basis_record = np.zeros([self.Nx*self.Ny*self.Nz*len(self.u.basis_frac), 4], dtype=np.intc)
        
    def maketranslation(self) :

        num_basis = len(self.u.basis_frac)

        frac_basis = np.copy(self.u.basis_frac)
        cell_basis = np.zeros([num_basis, 4], dtype=np.intc)
        cell_basis[:, 3] = np.arange(num_basis)

        num_atm = 0
        for nx in range(self.Nx) :
            for ny in range(self.Ny) :
                for nz in range(self.Nz) :

                    n = np.array([nx, ny, nz])
                    self.frac_coord[num_atm:num_atm+num_basis, :] = n + frac_basis
                    self.cell_basis_record[num_atm:num_atm+num_basis, :] = cell_basis + np.pad(n, (0, 1), 'constant')
                    num_atm += num_basis

        print()
        print("\tTotal Number of atoms: ", num_atm)

        self.cart_coord = np.matmul(self.frac_coord, self.u.lattice)
        #print()
        #print(self.frac_coord)
        #print(self.cart_coord)
        #print(self.cell_basis_record)


    def write_h5(self) :

        filename = 'record_'+self.u.prefix+'.h5'

        cell_basis_record_dset = 'cell_basis'
        frac_cord_dset = 'frac_coord'
        cart_cord_dset = 'cart_coord'

        h5f = h5py.File(filename, mode='w')

        h5f.create_dataset(cell_basis_record_dset, data=self.cell_basis_record)
        h5f.create_dataset(frac_cord_dset, data=self.frac_coord)
        h5f.create_dataset(cart_cord_dset, data=self.cart_coord)

        h5f.close()
        print("\tCell-basis info of atoms is written in: {0}".format(filename))

    def write_siesta_data(self) :
    
        Nat = self.Nx*self.Ny*self.Nz*len(self.u.basis_frac)
        cord = self.cart_coord

        filname = 'system_'+self.u.prefix+'.fdf'
        f2 = open(filname, "w+")
        f2.write("NumberOfSpecies           %d\n" % (self.u.ntyp)) 
        f2.write("NumberOfAtoms             %d\n" % (Nat))

        f2.write("\n")

        f2.write("%block ChemicalSpeciesLabel\n")

        for l in range(self.u.ntyp) :
            if ( self.u.ntyp == 1 ) :
                atomic_no = self.u.atomic_no_unq
                label = self.u.typ_lbl_unq
            else :
                atomic_no = self.u.atomic_no_unq[l]
                label = self.u.typ_lbl_unq[l]
            pseudo_file = label+".psml"
            f2.write("  %d  %d  %s    %s\n" %(l+1, atomic_no, label, pseudo_file))

        f2.write("%endblock ChemicalSpeciesLabel\n")
        f2.write("\n")

        f2.write("%block AtomicMass\n")

        for l in range(self.u.ntyp) :
            if ( self.u.ntyp == 1 ) :
                mass = self.u.mass_unq
            else :
                mass = self.u.mass_unq[l]
            f2.write("  %d  %f\n" %(l+1, mass))

        f2.write("%endblock AtomicMass\n")
        f2.write("\n")

        ALAT = 1.0
        f2.write("LatticeConstant           %f Ang\n" %ALAT)
        f2.write("\n")

        lat_vec = np.copy(self.u.lattice)
        lat_vec[0, :] *= self.Nx
        lat_vec[1, :] *= self.Ny
        lat_vec[2, :] *= self.Nz
        
        f2.write("%block LatticeVectors\n")
        f2.write("  %s  %s  %s\n" % (format(lat_vec[0, 0], '.13g'), \
                                     format(lat_vec[0, 1], '.13g'), \
                                     format(lat_vec[0, 2], '.13g')))
        f2.write("  %s  %s  %s\n" % (format(lat_vec[1, 0], '.13g'), \
                                     format(lat_vec[1, 1], '.13g'), \
                                     format(lat_vec[1, 2], '.13g')))
        f2.write("  %s  %s  %s\n" % (format(lat_vec[2, 0], '.13g'), \
                                     format(lat_vec[2, 1], '.13g'), \
                                     format(lat_vec[2, 2], '.13g')))
        f2.write("%endblock LatticeVectors\n")
        f2.write("\n")

        f2.write("AtomicCoordinatesFormat Ang\n")
        f2.write("\n")

        f2.write("%block AtomicCoordinatesAndAtomicSpecies\n")
        accu_float = 3+1+13
        line_format = "{x: .16E}   {y: .16E}   {z: .16E}   {typ: <4s} {Num: <4s} {typlbl: <4s} {comment: <3s} {cellIndx1: <4s} {cellIndx2: <4s} {cellIndx3: <4s} {basisType: <4s}\n"
        for j in range(Nat) :

            basis_indx = self.cell_basis_record[j, 3] 
            species_indx = self.u.typ_indx[basis_indx] 

            if ( self.u.ntyp == 1 ) :
                species_lbl  = self.u.typ_lbl_unq
            else :
                species_lbl = self.u.typ_lbl_unq[species_indx-1]

            f2.write( line_format.format(x=cord[j, 0], \
                                         y=cord[j, 1], \
                                         z=cord[j, 2], \
                                         typ=str(species_indx), \
                                         Num=str(j+1), \
                                         typlbl=species_lbl, \
                                         comment="#", \
                                         cellIndx1=str(self.cell_basis_record[j, 0]+1), \
                                         cellIndx2=str(self.cell_basis_record[j, 1]+1), \
                                         cellIndx3=str(self.cell_basis_record[j, 2]+1), \
                                         basisType=str(self.cell_basis_record[j, 3]+1)) )
        f2.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")
        f2.write("\n")
        f2.close
        print("\tSIESTA data file '{0}' written".format(filname))

      
    def write_lammps_data(self) :
    
        print()

        Nat = self.Nx*self.Ny*self.Nz*len(self.u.basis_frac)
        cord = self.cart_coord

        filname = 'data_'+self.u.prefix+'.lmp'
        f = open(filname, "w+")
        f.write("LAMMPS Atom File\n\n")

        f.write("%d atoms\n" % Nat)
        f.write("%d bonds\n" % 0)
        f.write("%d angles\n" % 0)
        f.write("%d dihedrals\n" % 0)
        f.write("%d impropers\n\n" % 0)

        f.write("%d atom types\n" % (self.u.ntyp))
        f.write("%d bond types\n" % 0)
        f.write("%d angles types\n\n" % 0)

        lat_vec = np.copy(self.u.lattice)
        lat_vec[0, :] *= self.Nx
        lat_vec[1, :] *= self.Ny
        lat_vec[2, :] *= self.Nz

        f.write("%s  %s  xlo xhi\n" % (format(0.0, '.13g'), format(lat_vec[0, 0], '.13g')))
        f.write("%s  %s  ylo yhi\n" % (format(0.0, '.13g'), format(lat_vec[1, 1], '.13g')))
        f.write("%s  %s  zlo zhi\n" % (format(0.0, '.13g'), format(lat_vec[2, 2], '.13g')))
        f.write("%s  %s  %s xy xz yz\n\n" % (format(lat_vec[1, 0], '.13g'), \
                                             format(lat_vec[2, 0], '.13g'), \
                                             format(lat_vec[2, 1], '.13g')))
        f.write("Masses\n\n")

        for l in range(self.u.ntyp) :
            if ( self.u.ntyp == 1 ) :
                mass = self.u.mass_unq
            else :
                mass = self.u.mass_unq[l]
            f.write("   %d  %f\n" % (l+1, mass))

        f.write("\n")
        f.write("Atoms\n\n")

        line_format = "{atm_no: <6s}   {mol_ID: <6s}   {atm_typ: <6s} {charge: <6s}  {x: .16E} {y: .16E} {z: .16E}\n"
        accu_float = 3+1+13
    
        file2 = self.u.prefix+'.xyz'
        f2 = open(file2, "w+")
        f2.write("%d\n" % Nat)
        f2.write("xyz formatted file\n")

        line_format2 = "{atm_lbl: <6s} {x: .16E} {y: .16E} {z: .16E}\n"

        mol_ID = 1
        for j in range(Nat) :

            basis_indx = self.cell_basis_record[j, 3]
            species_indx = self.u.typ_indx[basis_indx]

            if ( self.u.ntyp == 1 ) :
                species_lbl  = self.u.typ_lbl_unq
            else :
                species_lbl = self.u.typ_lbl_unq[species_indx-1]

            f2.write( line_format2.format(atm_lbl=species_lbl, x=cord[j, 0], y=cord[j, 1], z=cord[j, 2]) )

            f.write(line_format.format(atm_no=str(j+1), mol_ID=str(mol_ID), atm_typ=str(species_indx), charge=str(0.0), \
                    x=cord[j, 0], y=cord[j, 1], z=cord[j, 2]))

    
        f.write("\n")
        f.close()
        f2.close()
        print("\tLAMMPS data file \'", filname, "\' written")
    
    

if __name__ == "__main__" :

    parser = argparse.ArgumentParser()

    parser.add_argument("-inp", type=str, default='input.nml', \
            help="Enter the name of input file")
    parser.add_argument("-MD", type=str, default='siesta', \
            help="Enter MD interface (siesta / lammps)")

    args = parser.parse_args()
    fname = args.inp
    MD = args.MD

    S = MakeSupCell(fname)
    S.maketranslation()
    S.write_h5()

    if ( MD == "siesta" ) :
        S.write_siesta_data()

    elif ( MD == "lammps" ) :
        S.write_lammps_data()

    else :

        print("ERROR: Unrecognized MD interface: ", MD)
        print("Only siesta and lammps interface available")

