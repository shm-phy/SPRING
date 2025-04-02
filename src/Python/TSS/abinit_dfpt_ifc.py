#!/usr/bin/python3

import argparse
import numpy as np
import f90nml
import h5py

import os
import shutil


BohrToAng = 0.5291772108 #1 Bohr = 0.5291772108 Angstroms
HaToeV = 27.2113845 #1 Ha = 27.2113845 eV
# There is no mention of the unit of IFCs generated from abinit in the documentation.
# Considering it is Ha/(Bohr**2)
# Ha/(Bohr**2) = 27.2113845/(0.5291772108**2) eV/(Ang**2) = 97.173618095 eV/(Ang**2)
Ha_Bh2ToeV_Ang2 = 97.173618095 # 1 Ha/(Bohr**2) = 97.173618095 eV/(Ang**2)

# Custom converter function to handle 'D' notation
def custom_converter(value):
    try:
        # Replace 'D' with 'e' and convert to float
        return float(value.decode('utf-8').replace('D', 'e'))
    except ValueError:
        # If conversion fails, return a default value (e.g., NaN)
        return np.nan


def insert_lines_after_line_number(file_path, line_number, new_lines):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        index_to_insert_after = line_number - 1

        with open(file_path, 'w') as file:
            file.writelines(lines[:index_to_insert_after + 1])

            file.writelines(new_lines)

            file.writelines(lines[index_to_insert_after + 1:])

    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"Error writing to the file: {e}")


def FindLineNumber(lookup_word, filepath, flag) :

    fl = open(filepath)

    num_lines = sum(1 for line in fl)
    read = fl.read()
    fl.seek(0)
    arr = []
    line_num = []
    for i in range(num_lines):
        arr.append(fl.readline())

    for i in range(len(arr)):

        if ( flag == "old" ) :

            if lookup_word in arr[i]:
                line_num.append(i+1)
        else:

            if f" {lookup_word} " in f" {arr[i]} ":
                line_num.append(i+1)

    if (len(line_num) != 1) :
        print("\t*** Error in reading file:", filepath, "***", len(line_num))

    fl.close()

    return line_num[0]


def Read_abinit_variable( file_path, line_num, column = 1, ncols=1, dtype="integer" ) :

    try:
        with open(file_path, 'r') as file:
            for i, line in enumerate(file, 1):
                if i == line_num:
                    columns = line.strip().split()
                    if len(columns) >= 2:
                        
                        if ( dtype == "integer" ) :
                            #var = int( columns[column] )
                            var = [int(item) for item in columns[column:column+ncols]]
                            if ( ncols == 1 ) :
                                var = var[0]
                            else :
                                var = np.array( var )

                        else :
                            #var = float( columns[column] )
                            var = [float(item) for item in columns[column:column+ncols]]
                            if ( ncols == 1 ) :
                                var = var[0]
                            else :
                                var = np.array( var )

                        return var

                    else:
                        print(f"Line {linenum} does not have at least two columns.")
                        return None
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None
    except Exception as e:
        print(f"Error reading the file: {e}")
        return None


def fold_indx( dimen, indx ) :

    fold = np.copy( indx )

    dimen_2 = dimen // 2

    for nn in range(3) :
        if ( fold[nn] > dimen_2[nn] ) :
            fold[nn] = fold[nn] - dimen[nn]

    return fold


class Unit:

    def __init__(self, nml_fname='input.nml') :

        #print("\nReading input from file: '", nml_fname, "'\n")
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


class AbinitDFPT_out :

    def __init__( self, prefix, nml_fname, T ) :

        self.u = Unit(nml_fname)
        self.nml_file = nml_fname
        self.prefix = prefix
        self.T = T

        self.B_Chrg = self.ReadBornCharge( )
        self.dielec = self.ReadDielectricTensor( )

        self.frac_coord, self.cart_coord, self.cell_basis_record = self.CreateSupCell( )

        self.FC2, self.atm_Indx, self.NumAtoms = self.ReadIFCs( )
        self.FC2 = self.FC2 * Ha_Bh2ToeV_Ang2 #Unit Conversion
        self.checkASR()
        self.write_h5()
        self.writeDielecBZ_nml()


    def writeDielecBZ_nml( self ) :

        file_path = self.prefix+"_spring.nml"
        lookup_string = "!Write dielectric and Born charge"
        line_number = FindLineNumber( lookup_string, file_path, flag="old" )

        for i in range(3) :
            string = "    dielec(:, "+str(i+1)+")    = "+"{:.13f}".format(self.dielec[i, 0])+"    {:.13f}".format(self.dielec[i, 1])+\
                                                       "    {:.13f}".format(self.dielec[i, 2])+"\n"
            insert_lines_after_line_number(file_path, line_number, string)
            line_number += 1

        string = "\n"
        insert_lines_after_line_number(file_path, line_number, string)
        line_number += 1

        for i in range(self.u.nbasis) :
            for alpha in range(3) :
                string = "    BornZ(:, "+str(alpha+1)+", "+str(i+1)+")  =  "+"{:.13f}".format(self.B_Chrg[i, alpha, 0])+\
                                                                              "    {:.13f}".format(self.B_Chrg[i, alpha, 1])+\
                                                                              "    {:.13f}".format(self.B_Chrg[i, alpha, 2])+"\n"
                insert_lines_after_line_number(file_path, line_number, string)
                line_number += 1

            string = "\n"
            insert_lines_after_line_number(file_path, line_number, string)
            line_number += 1


    def write_h5(self) :

        filename = "FC2nd_"+self.prefix+str(T)+"K_F.h5"

        fc_dset = 'FC2'
        cell_basis_record_dset = 'atm_Indx'
        atmNum_dset = 'NumAtoms'
        born_eff_chrd_dset = "BornEffectiveCharge"
        dielec_dset = "DielectricTensor"

        h5f = h5py.File(filename, mode='w')

        h5f.create_dataset(fc_dset, data=self.FC2)
        h5f.create_dataset(cell_basis_record_dset, data=self.atm_Indx)
        h5f.create_dataset(atmNum_dset, data=self.NumAtoms)
        h5f.create_dataset( born_eff_chrd_dset, data=self.B_Chrg )
        h5f.create_dataset( dielec_dset, data=self.dielec )

        h5f.close()
        print("\tIFCs from abinit converted to spring format and written in: {0}".format(filename))



    def CreateSupCell( self, fold=True ) :

        num_basis = len(self.u.basis_frac)

        file_path = self.prefix+"_anaddb.abo"
        lookup_word = "ngqpt"
        line_num = FindLineNumber( lookup_word, file_path, flag="new" )
        ngqpt = Read_abinit_variable( file_path, line_num, column=1, ncols=3, dtype="integer" )

        if ( not(np.array_equal(ngqpt, self.u.sup_dim)) ) :
            print(f"\tWARNING: ngqpt ({ngqpt}) and Sup_Cell ({self.u.sup_dim}) are not same")

        Nx = ngqpt[0] #self.u.sup_dim[0]
        Ny = ngqpt[1] #self.u.sup_dim[1]
        Nz = ngqpt[2] #self.u.sup_dim[2]

        frac_basis = np.copy(self.u.basis_frac)
        cell_basis = np.zeros([num_basis, 4], dtype=np.intc)
        cell_basis[:, 3] = np.arange(num_basis)

        frac_coord = np.zeros( [Nx*Ny*Nz*len(self.u.basis_frac), 3], dtype=np.double )
        cart_coord = np.zeros_like( frac_coord )
        cell_basis_record = np.zeros( [Nx*Ny*Nz*len(self.u.basis_frac), 4], dtype=np.intc )

        num_atm = 0
        for nx in range(Nx) :
            for ny in range(Ny) :
                for nz in range(Nz) :

                    n = np.array([nx, ny, nz])

                    if ( fold ) :
                        n_fold = fold_indx( self.u.sup_dim, n )

                    else :
                        n_fold = n

                    #print( n, " ==> ", n_fold )

                    frac_coord[num_atm:num_atm+num_basis, :] = n_fold + frac_basis
                    cell_basis_record[num_atm:num_atm+num_basis, :] = cell_basis + np.pad(n_fold, (0, 1), 'constant')
                    num_atm += num_basis

        #print()
        #print("\tTotal Number of atoms: ", num_atm)

        cart_coord = np.matmul(frac_coord, self.u.lattice) / BohrToAng

        #print( "frac_cord : " )
        #print(frac_coord)

        #print ( "\ncart_cord : " )
        #print(cart_coord)

        #print( "\ncell_basis_record: " )
        #print(cell_basis_record)

        #for ii in range( num_atm ) :
        #    print( cell_basis_record[ii, :], "  ==>  ", frac_coord[ii, :], "  ==>  ", cart_coord[ii, :] )

        return frac_coord, cart_coord, cell_basis_record


    def ReadBornCharge( self ) :

        file_path = self.prefix+"_anaddb.abo"
        lookup_string = "Effective charge tensors"
        line_number = FindLineNumber( lookup_string, file_path, flag="old" )
        line_number += 3

        atom_cart = np.zeros([self.u.nbasis*3, 2], dtype=np.intc)

        atom_cart = np.loadtxt(file_path, dtype=np.intc, skiprows=line_number, usecols=(0, 1), max_rows=self.u.nbasis*3)

        B_Z = np.loadtxt(file_path, dtype=np.double, skiprows=line_number, usecols=(2, 3, 4), max_rows=self.u.nbasis*3)

        Born_Eff_Chrg = np.zeros([self.u.nbasis, 3, 3], dtype=np.double)

        for i in range( self.u.nbasis*3 ) :

            atom = atom_cart[i, 0]
            cart = atom_cart[i, 1]

            Born_Eff_Chrg[atom-1, cart-1, :] = B_Z[i, :]

        return Born_Eff_Chrg


    def ReadDielectricTensor( self ) :

        file_path = self.prefix+"_dfpt.abo"
        lookup_string = "Dielectric tensor, in cartesian coordinates,"
        line_number = FindLineNumber( lookup_string, file_path, flag="old" )
        line_number += 3

        alpha_beta = np.zeros([9, 2], dtype=np.intc)
        alpha_beta = np.loadtxt(file_path, dtype=np.intc, skiprows=line_number, usecols=(0, 2), max_rows=11)

        dielc = np.zeros(9, dtype=np.double)
        dielc = np.loadtxt(file_path, dtype=np.double, skiprows=line_number, usecols=(4), max_rows=11)

        dielec_T = np.zeros([3, 3], dtype=np.double)

        for i in range( 9 ) :
            alpha = alpha_beta[i, 0]
            beta  = alpha_beta[i, 1]

            dielec_T[alpha-1, beta-1] = dielc[i]

        #Find the pressure prediction from abinit
        lookup_string = "pressure_GPa:"
        line_num = FindLineNumber( lookup_string, file_path, flag="old" )
        P = Read_abinit_variable( file_path, line_num, column=1, ncols=1, dtype="float" )
        print(f"\tPressure as predicted from abinit: {P} GPa")

        return dielec_T 


    def ReadIFCs( self ) :

        file_path = self.prefix+"_anaddb.abo"

        lookup_word = "natifc"
        line_num = FindLineNumber( lookup_word, file_path, flag="old" )
        natifc = Read_abinit_variable( file_path, line_num, column=1, ncols=1, dtype="integer" )

        lookup_word = "atifc"
        line_num = FindLineNumber( lookup_word, file_path, flag="new" )
        atifc = Read_abinit_variable( file_path, line_num, column=1, ncols=natifc, dtype="integer" )
        atifc_l = np.zeros( natifc, dtype=np.intc )
        atifc_l[:] = atifc

        lookup_word = "ifcout"
        line_num = FindLineNumber( lookup_word, file_path, flag="old" )
        numberOfatoms = Read_abinit_variable( file_path, line_num, column=1, ncols=1, dtype="integer" )

        nu_list = np.zeros( [natifc, numberOfatoms], dtype=np.intc )
        nu_atom_pos = np.zeros( [natifc, numberOfatoms, 3], dtype=np.double )
        FC2 = np.zeros( [natifc, numberOfatoms, 3, 3, 1], dtype=np.double )
        atm_Indx = np.zeros( [natifc, numberOfatoms, 4], dtype=np.intc )
        NumAtoms = np.zeros( natifc, dtype=np.intc )

        for atom in range( natifc ) :

            NumAtoms[atom] = numberOfatoms
            mu = atifc_l[atom]

            mu_indx = 0
            if (natifc != self.u.nbasis)  :
                print(f"WARNING: Number of basis atoms ( {self.u.nbasis} ) != natifc ( {natifc} ) ")
                print(f"         Please change the 'natifc' input variable in anaddb input file ")
                mu_indx += 1
            else :
                mu_indx = mu

            lookup_string = "generic atom number   "+str(mu)
            line_num = FindLineNumber( lookup_string, file_path, flag="old" )
            line_num += 6

            for nuatom in range( numberOfatoms ) :

                nu = Read_abinit_variable( file_path, line_num, column=4, ncols=1, dtype="integer" )
                nu_list[mu_indx-1, nuatom] = nu

                R_nu = np.loadtxt( file_path, dtype=np.double, skiprows=line_num, usecols=(2, 3, 4), max_rows=1)
                nu_atom_pos[mu_indx-1, nuatom, :] = R_nu

                #FC_3_3 = np.loadtxt( file_path, dtype=np.double, skiprows=line_num+2, usecols=(0, 1, 2), max_rows=3) #Just test
                FC_3_3 = np.loadtxt( file_path, dtype=np.double, skiprows=line_num+2, usecols=(6, 7, 8), max_rows=3)

                index = self.FindCellBasisIndex( R_nu, nu )
                atm_Indx[mu_indx-1, nuatom, :] = self.cell_basis_record[index, :]
                atm_Indx[mu_indx-1, nuatom, 3] = nu

                for alpha in range(3) :
                    for beta in range(3) :

                        fc = FC_3_3[beta, alpha]
                        FC2[mu_indx-1, nuatom, alpha, beta, 0] = fc

                line_num += 21

        #print( line_num )
        #print( nu_list )
        #print(nu_atom_pos[0, :, :])
        #print( FC2 )
        #print( atm_Indx )
        #print( NumAtoms )

        return FC2, atm_Indx, NumAtoms


    def checkASR( self, adhoc_asr=False ) :

        ZERO = 1.0E-8

        IFC_sum = np.sum( self.FC2 )

        if ( abs(IFC_sum) >= ZERO ) :
            print(f"\tWARNING:")
            print(f"\tAcoustic sum rule (ASR) is not maintained for the abinit generated IFCs : {IFC_sum} (eV/Ang**2)")
            print(f"\tWe don't have complete understanding of the format of abinit IFCs")
            print(f"\tPlease check the phonon dispresion from: disp.x -inp {self.nml_file} -T {self.T}")

            if ( adhoc_asr ) :
                print(f"\tTrying some crude ASR imposition")
                self.FC2 = self.FC2 - (IFC_sum / np.size( self.FC2 ))

        print()



    def FindCellBasisIndex( self, R_nu, nu, accuracy=1.0E-4 ) :

        matching_rows = np.all(np.isclose(self.cart_coord, R_nu, rtol=accuracy, atol=accuracy), axis=1)

        indices = np.where(matching_rows)[0]

        if (len(indices) == 1) :
            index = indices[0]
            #print(f"Row with specified accuracy found at index: {index} and nu: {nu-1}")

            if ( (nu-1) != self.cell_basis_record[index, 3] ) :
                print( f"ERROR: nu ( {nu-1} ) does not match with nu of self.cell_basis_record ( {self.cell_basis_record[index, 3]} ) " )

            return index

        elif (len(indices) > 1) :
            print(f"WARNING: Multiple matching rows found. Ambiguous result. {indices}")
            return None

        else:
            print("No matching row found.")
            return None



class AbinitDFPT_in :

    def __init__(self, prefix) :

        self.prefix = prefix

        self.nkpt = self.NumberofKPoints()
        print( f"Number of k points = {self.nkpt}" )

        self.kpt = self.KPoinTs( )
        #print( self.kpt )

        self.WrtieQpoints()
        self.MergeDDB_inputFile( )
        self.anaddb_input_file( )


    def anaddb_input_file( self ) :

        file_path = self.prefix+"_anaddb.files"

        try:
            # Open the file in write mode ('w')
            with open(file_path, 'w') as file:
                # Write the content to the file
                content = self.prefix+"_anaddb.abi\n"
                content += self.prefix+"_anaddb.abo\n"
                content += self.prefix+"_mrgddb.abo\n"
                content += self.prefix+"_band2eps\n"
                content += self.prefix+"_dummy1\n"
                content += self.prefix+"_dummy2\n"
                content += self.prefix+"_dummy3\n"

                file.write(content)

            print(f"anaddb input written to {file_path} successfully.")

        except Exception as e:
            print(f"Error writing to the file: {e}")



    def MergeDDB_inputFile( self ) :

        file_path = self.prefix+"_mrgddb.abi"

        try:
            # Open the file in write mode ('w')
            with open(file_path, 'w') as file:
                # Write the content to the file
                content = self.prefix+"_mrgddb.abo\n"
                content += self.prefix+" phonons on 5 5 5 mesh\n"
                content += str(self.nkpt)+"\n"

                file.write(content)

                #ZnS_dfpto_DS3_DDB
                for i in range( self.nkpt ) :

                    content = self.prefix+"_dfpto_DS"+str(i+3)+"_DDB\n"
                    file.write( content )

            print(f"mrgddb input written to {file_path} successfully.")

        except Exception as e:
            print(f"Error writing to the file: {e}")


    def WrtieQpoints( self ) :

        line_format = "  {qpt_str: <6s}       {qx: <19s} {qy: <19s} {qz: <19s}\n"

        accu_float = 3+1+13

        qpt_dset_str = "qpt2"
        zero = 0.0
        new_lines = line_format.format( qpt_str=qpt_dset_str, qx=str(zero)[:accu_float], \
                                                              qy=str(zero)[:accu_float], \
                                                              qz=str(zero)[:accu_float] )
        for i in range( self.nkpt ) :

            qpt_dset_str = "qpt"+str(i+3)
            new_lines += line_format.format( qpt_str=qpt_dset_str, qx=str(self.kpt[i, 0])[:accu_float], \
                                                                   qy=str(self.kpt[i, 1])[:accu_float], \
                                                                   qz=str(self.kpt[i, 2])[:accu_float] )
        #print( new_line )

        lookup_word = "#--WRITE_Qs--#"
        file_path = self.prefix+"_dfpt.abi"
        line_number = FindLineNumber( lookup_word, file_path, flag="old" )

        insert_lines_after_line_number(file_path, line_number, new_lines)

        lookup_word = "#--WRITE_NDSET--#"
        file_path = self.prefix+"_dfpt.abi"
        line_number = FindLineNumber( lookup_word, file_path, flag="old" )

        new_lines = "  ndtset  "+str(2+self.nkpt)+"\n"
        insert_lines_after_line_number(file_path, line_number, new_lines)

        print(f"qpoints written to {file_path} successfully.")



    def KPoinTs( self ) :

        lookup_word = "kpt"
        filepath = self.prefix+"_qpointso_DDB"

        linenum = FindLineNumber( lookup_word, filepath, flag="new" )

        kpt = np.zeros([self.nkpt, 3], dtype=np.double)
        kpt[0, :] = np.genfromtxt(filepath, dtype=float, \
                                  converters={1: custom_converter, 2: custom_converter, 3: custom_converter}, \
                                  skip_header=linenum-1, usecols=(1,2,3), max_rows=1, encoding='bytes')
        kpt[1:, :] = np.genfromtxt(filepath, dtype=float, \
                                  converters={0: custom_converter, 1: custom_converter, 2: custom_converter}, \
                                  skip_header=linenum, usecols=(0,1,2), max_rows=self.nkpt-1, encoding='bytes')

        return kpt

    def NumberofKPoints( self ) :

        lookup_word = "nkpt"
        filepath = self.prefix+"_qpointso_DDB"

        linenum = FindLineNumber( lookup_word, filepath, flag="new" )

        nkpt = Read_abinit_variable( filepath, linenum, column = 1, ncols=1, dtype="integer" )

        return nkpt


class LatVec :

    def __init__(self, prefix, nml_file, a, b, c, angl_bc, angl_ac, angl_ab) :

        self.prefix = prefix
        self.a = a
        self.b = b
        self.c = c
        self.angl_bc = angl_bc
        self.angl_ac = angl_ac
        self.angl_ab = angl_ab
        self.nml_file = nml_file
        self.pwd = os.getcwd()
        self.WriteLatVec()


    def WriteLatVec(self) :
    
        L1 = np.zeros(3, dtype=np.double)
        L2 = np.zeros(3, dtype=np.double)
        L3 = np.zeros(3, dtype=np.double)
    
        L1[0] = self.a
    
        L2[0] = self.b * np.cos( self.angl_ab )
        L2[1] = self.b * np.sin( self.angl_ab )
    
        L3[0] = self.c * np.cos( self.angl_ac )
        L3[1] = self.c * ( (np.cos(self.angl_bc) - np.cos(self.angl_ab) * np.cos(self.angl_ac))  / np.sin(self.angl_ab) )
        L3[2] = np.sqrt( self.c**2 - L3[0]**2 - L3[1]**2 )

        vec1 = "rprim    "+"   {:.13f}".format(L1[0])+"    {:.13f}".format(L1[1])+"    {:.13f}".format(L1[2])
        vec2 = "         "+"   {:.13f}".format(L2[0])+"    {:.13f}".format(L2[1])+"    {:.13f}".format(L2[2])
        vec3 = "         "+"   {:.13f}".format(L3[0])+"    {:.13f}".format(L3[1])+"    {:.13f}".format(L3[2])

        # ---------------------------------------------------------------------------- #
        abinit_file1 = self.pwd+"/"+self.prefix+"_qpoints.abi"
        fin = open(abinit_file1, "rt")
        data = fin.read()
        data = data.replace('REPLACE_STRING_WITH_Latvec1', vec1)
        data = data.replace('REPLACE_STRING_WITH_Latvec2', vec2)
        data = data.replace('REPLACE_STRING_WITH_Latvec3', vec3)
        fin.close()
        
        fin = open(abinit_file1, "wt")
        fin.write(data)
        fin.close()
        # ---------------------------------------------------------------------------- #

        # ---------------------------------------------------------------------------- #
        abinit_file2 = self.pwd+"/"+self.prefix+"_dfpt.abi"
        fin = open(abinit_file2, "rt")
        data = fin.read()
        data = data.replace('REPLACE_STRING_WITH_Latvec1', vec1)
        data = data.replace('REPLACE_STRING_WITH_Latvec2', vec2)
        data = data.replace('REPLACE_STRING_WITH_Latvec3', vec3)
        fin.close()
        
        fin = open(abinit_file2, "wt")
        fin.write(data)
        fin.close()
        # ---------------------------------------------------------------------------- #

        vec1 = "Latvec(:, 1)    = "+"   {:.13f}".format(L1[0])+"    {:.13f}".format(L1[1])+"    {:.13f}".format(L1[2])
        vec2 = "Latvec(:, 2)    = "+"   {:.13f}".format(L2[0])+"    {:.13f}".format(L2[1])+"    {:.13f}".format(L2[2])
        vec3 = "Latvec(:, 3)    = "+"   {:.13f}".format(L3[0])+"    {:.13f}".format(L3[1])+"    {:.13f}".format(L3[2])

        # ---------------------------------------------------------------------------- #
        spring_in_file = self.pwd+"/"+self.nml_file
        fin = open(spring_in_file, "rt")
        data = fin.read()
        data = data.replace('REPLACE_STRING_WITH_Latvec(:, 1)', vec1)
        data = data.replace('REPLACE_STRING_WITH_Latvec(:, 2)', vec2)
        data = data.replace('REPLACE_STRING_WITH_Latvec(:, 3)', vec3)
        fin.close()

        fin = open(spring_in_file, "wt")
        fin.write(data)
        fin.close()
        # ---------------------------------------------------------------------------- #



if __name__ == "__main__" :

    parser = argparse.ArgumentParser()

    parser.add_argument("-prefix", type=str, default='ZnS', help="Enter the prefix")
    parser.add_argument("-inp", type=str, default='input.nml', help="Enter the name of input file")
    parser.add_argument("-T", type=float, default=300.0, help="Temperature in Kelvin ")

    parser.add_argument("-a", type=float, default=1.0, help="Lattice constant a in Anstrom")
    parser.add_argument("-b", type=float, default=1.0, help="Lattice constant b in Anstrom")
    parser.add_argument("-c", type=float, default=1.0, help="Lattice constant c in Anstrom")

    parser.add_argument("-alpha", type=float, default=60.0, help="Angle alpha (between b and c) in degree")
    parser.add_argument("-beta", type=float, default=60.0, help="Angle beta (between a and c) in degree")
    parser.add_argument("-gama", type=float, default=60.0, help="Angle gama (between a and b) in degree")

    parser.add_argument("--abinit-in", dest='abinit_in_flag', action='store_true', \
            help="Enter boolean value to create input file for abinit dfpt")
    parser.add_argument("--abinit-out", dest='abinit_out_flag', action='store_true', \
            help="Enter boolean value to create IFCs in spring format from abinit dfpt output")
    parser.add_argument("--lat_vec_write", dest='lat_vec_flag', action='store_true', \
            help="Enter boolean value to write lattice vectors in input files")

    parser.set_defaults(abinit_in_flag=False)
    parser.set_defaults(abinit_out_flag=False)
    parser.set_defaults(lat_vec_flag=False)

    args = parser.parse_args()
    prefix = args.prefix
    nml_fname = args.inp
    T = args.T

    a = args.a
    b = args.b
    c = args.c
    alpha = args.alpha
    beta = args.beta
    gama = args.gama
    
    abinit_in_flag = args.abinit_in_flag
    abinit_out_flag = args.abinit_out_flag
    lat_vec_flag = args.lat_vec_flag

    if ( lat_vec_flag ) :
        angl_bc = alpha * np.pi / 180.0
        angl_ac = beta * np.pi / 180.0
        angl_ab = gama * np.pi / 180.0

        LatVec_write = LatVec(prefix, nml_fname, a, b, c, angl_bc, angl_ac, angl_ab)

    elif ( abinit_in_flag ) :
        DataDFPT_in = AbinitDFPT_in(prefix)

    elif ( abinit_out_flag ) :
        DataDFPT_out = AbinitDFPT_out( prefix, nml_fname, T )

