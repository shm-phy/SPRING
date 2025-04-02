
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



class ReadSIESTA :

    def __init__(self, nml_fname, siestafile, mde_file, ts_tot, \
                 Tmean, ts_cut, dT=50, d_time=20, init_force_zero=False) :

        self.u = Unit(nml_fname)

        self.N = np.copy(self.u.sup_dim)
        print("\tTemperature of MD run:", Tmean)
        print("\tNumber of timestep of MD run:", total_ts)
        print("\tSuper cell dimension:", self.N)

        self.Natom = self.N[0]*self.N[1]*self.N[2]*len(self.u.basis_frac)

        self.recordfile = 'record_'+self.u.prefix+'.h5'
        self.ts_tot = ts_tot
        self.siestafile = siestafile
        self.mdefile = mde_file
        self.cell_basis, self.cart_cord, self.frac_cord = self.ReadRecord5()

        self.Nstep, self.timestep = self.ReadMDE(Tmean, dT, ts_cut, d_time)

        Nbasis = len(self.u.basis_frac)

        print("\t*** Reading output of MD from SIESTA (displacement and force) '", self.siestafile, "'***")
        self.initForc = self.InitForce(init_force_zero)

        self.dU = np.zeros([self.Nstep, self.N[0], self.N[1], self.N[2], Nbasis, 3], dtype=np.double)

        self.forc = np.zeros([self.Nstep, self.N[0], self.N[1], self.N[2], Nbasis, 3], dtype=np.double)

        lat_vec = np.copy(self.u.lattice)
        self.lat_len = np.linalg.norm(lat_vec, axis=1)
        unit = 1.0 / self.lat_len
        lat_vec = lat_vec * unit[:, None]

        self.lat_unit_inv = np.linalg.inv(lat_vec)
        self.lat_len *= self.N

    
    def ReadMDE(self, mean_temp, dT, time_step_cut, dt) :

        print("\t*** Reading MDE file of SIESTA MD output '", self.mdefile, "'***")
        mde_dat = np.loadtxt(self.mdefile)
        mde_dat = mde_dat[np.where( mde_dat[:, 1] < mean_temp+dT )]
        mde_data = mde_dat[np.where( mde_dat[:, 1] > mean_temp-dT )]

        step = mde_data[:, 0]
        step = step[::-1].astype(int)

        timestep = []
        num = 0

        for n in step :

            if n > time_step_cut :
                if num == 0 :
                    timestep.append(n)
                    num += 1

                elif (num != 0) and abs(n - timestep[num-1]) > dt :
                    timestep.append(n)
                    num += 1

            elif n <= time_step_cut :
                break

        print()
        print("\tNumber of MD snapshots satisfy the given condition: ", len(timestep))

        return len(timestep), timestep


    def ReadRecord5(self) :

        cell_basis_record_dset = 'cell_basis'
        frac_cord_dset = 'frac_coord'
        cart_cord_dset = 'cart_coord'

        print("\t*** Reading record of Super-Cell information '", self.recordfile, "'***")
        h5f = h5py.File(self.recordfile, mode='r')

        cell_basis = h5f[cell_basis_record_dset][:]
        cart_cord = h5f[cart_cord_dset][:]
        frac_cord = h5f[frac_cord_dset][:]

        h5f.close()

        return cell_basis, cart_cord, frac_cord


    def FindLineNumber(self, lookup_word) :

        fl = open(self.siestafile)

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

        if ( len(line_num) > self.ts_tot ) :
            del line_num[-1]

        fl.close()

        return line_num


    def ReadCord(self) :

        line = self.FindLineNumber("outcoor")
        tstep = self.timestep

        num_step = 0
        print()
        for ts in tstep :

            ln = line[ts-1]

            cord = np.loadtxt(self.siestafile, dtype=np.double, skiprows=ln, \
                              usecols=(0, 1, 2), max_rows=self.Natom)

            dL = (cord - self.cart_cord)

            for atm in range(self.Natom) :

                disp_xyz = np.copy(dL[atm, :])

                dsp_xyz = self.CheckBound( atm, disp_xyz, np.copy(cord[atm, :]) )

                nx, ny, nz, mu = self.cell_basis[atm, :]
                self.dU[num_step, nx, ny, nz, mu, :] = dsp_xyz #np.copy(dL[atm, :])

            num_step += 1
            print("\tRead Coordinate: No. of MD snapshot covered:", num_step, "/", self.Nstep, \
                  "/|\ Current line number reading in file: ", ln, end='\r', flush=True)



    def CheckBound(self, atm, dsp_xyz, crd_xyz) :

        dsp_abc = np.matmul(dsp_xyz, self.lat_unit_inv)

        if np.any( abs(dsp_abc) >= (self.lat_len / 2.0) ) :
            print("*************** Atom moves out of box ***********************")
            print("PBC Wrap needed")
            lat_vec = np.copy(self.u.lattice)

            for xx in range(3) :

                if ( abs(dsp_abc[xx]) >= (self.lat_len[xx] / 2.0) ) :
                    if dsp_abc[xx] < 0.0 :
                        crd_xyz += lat_vec[xx, :]
                    elif dsp_abc[xx] > 0.0 :
                        crd_xyz -= lat_vec[xx, :]

            cord_zyz0 = np.copy(self.cart_cord[atm, :])
            dsp_xyz_new = crd_xyz - cord_xyz0

            return dsp_xyz_new

        else :
            return dsp_xyz



    def InitForce(self, init_force_zero) :

        if ( init_force_zero ) :
            initForc = np.zeros([self.Natom, 3], dtype=np.double)

        else :

            line = self.FindLineNumber("(eV/Ang):")

            init = line[0]

            initForc = np.loadtxt(self.siestafile, dtype=np.double, skiprows=init, \
                                  usecols=(1, 2, 3), max_rows=self.Natom)

        return initForc


    def ReadForce(self) :

        line = self.FindLineNumber("(eV/Ang):")

        tstep = self.timestep

        num_step = 0
        print()
        for ts in tstep :

            ln = line[ts-1]

            force = np.loadtxt(self.siestafile, dtype=np.double, skiprows=ln, \
                              usecols=(1, 2, 3), max_rows=self.Natom)

            force = (force-self.initForc)

            for atm in range(self.Natom) :
                nx, ny, nz, mu = self.cell_basis[atm, :]
                self.forc[num_step, nx, ny, nz, mu, :] = np.copy(force[atm, :])

            num_step += 1
            print("\tRead Force: No. of MD snapshot covered:", num_step, "/", self.Nstep, \
                  "/|\ Current line number reading in file: ", ln, end='\r', flush=True)

    
    def Writeh5(self, T) :

        filename = 'disp_forc_'+self.u.prefix+str(T)+'K.h5'

        print()
        print("\t*** Writing force-displacement dataset in file:'", filename, "'***")

        disp_dset = 'disp_dataset'
        force_dset = 'force_dataset'

        h5f = h5py.File(filename, mode='w')

        h5f.create_dataset(disp_dset, data=self.dU)
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


class Read_lammps :

    def __init__(self, nml_fname, lammps_init_force_file, lammps_trj_file, overwrite_initPos=False ) :

        self.u = Unit(nml_fname)

        self.N = np.copy(self.u.sup_dim)
        print("\tSuper cell dimension:", self.N)

        self.Natom = self.N[0]*self.N[1]*self.N[2]*len(self.u.basis_frac)

        self.recordfile = 'record_'+self.u.prefix+'.h5'

        if ( lammps_init_force_file == 'default' ) :
            self.init_frc_file = self.u.prefix+'_init.lammpstrj'
        else :
            self.init_frc_file = lammps_init_force_file

        if ( lammps_trj_file == 'default' ) :
            self.trj_file = self.u.prefix+'.lammpstrj'
        else :
            self.trj_file = lammps_trj_file

        print("\tInitial force of unperturbed system will be read from file: ", self.init_frc_file)
        print("\tLammps trajectory will be read from file: ", self.trj_file)

        self.cell_basis, self.cart_cord, self.frac_cord = self.ReadRecord5()

        self.Nstep, self.cord, self.force = self.ReadTRAJ()

        Nbasis = len(self.u.basis_frac)

        self.initPos, self.initForc = self.InitForce()

        self.dU = np.zeros([self.Nstep, self.N[0], self.N[1], self.N[2], Nbasis, 3], dtype=np.double)

        self.forc = np.zeros([self.Nstep, self.N[0], self.N[1], self.N[2], Nbasis, 3], dtype=np.double)

        #lat_vec = np.copy(self.u.lattice)
        #self.lat_len = np.linalg.norm(lat_vec, axis=1)
        #unit = 1.0 / self.lat_len
        #lat_vec = lat_vec * unit[:, None]

        #self.lat_unit_inv = np.linalg.inv(lat_vec)
        #self.lat_len *= self.N

        lat_vec = np.copy(self.u.lattice)
        #print( lat_vec )
        self.lat_len = np.linalg.norm(lat_vec, axis=1)
        unit = 1.0 / self.lat_len
        #print( unit )
        lat_vec = lat_vec * unit[:, None]
        #print( lat_vec )

        self.lat_unit_inv = np.linalg.inv(lat_vec)
        self.lat_len *= self.N
        #print( self.lat_len )

        lat_vec = np.copy(self.u.lattice)
        self.lat_vec = lat_vec * self.N[:, None]
        #print( lat_vec )


    def OverWriteInitPos_h5(self):

        print("Over-writing initial position of the system ...")

        cart_cord_dset = 'cart_coord'
        initial_force_dset = 'init_force'

        # Delete
        try:
            with h5py.File(self.recordfile, 'r+') as file:
                if cart_cord_dset in file:
                    cart_cord_save= h5f[cart_cord_dset][:]
                    del file[cart_cord_dset]
                    print(f"Dataset '{cart_cord_dset}' deleted successfully.")
                else:
                    print(f"Dataset '{cart_cord_dset}' not found in the file.")
                if initial_force_dset in file:
                    del file[initial_force_dset]
                    print(f"Dataset '{initial_force_dset}' deleted successfully.")
        except Exception as e:
            print(f"Error occurred: {e}")

        #Write
        h5f = h5py.File(self.recordfile, mode='r+')

        print( np.max( abs(self.initPos-cart_cord_save) ) )
        h5f.create_dataset(cart_cord_dset, data=self.initPos)
        h5f.create_dataset(initial_force_dset, data=self.initForc)

        h5f.close()

    
    def ReadTRAJ(self) :

        print("\t*** Reading trajectory file of lammps MD output '", self.trj_file, "'***")

        Nstep_guess = 40

        cord = np.zeros([Nstep_guess, self.Natom, 3], dtype=np.double)
        force = np.zeros([Nstep_guess, self.Natom, 3], dtype=np.double)

        Nstep = 0
        flag = 1

        while ( flag == 1 ) :

            skip_line = Nstep * ( 9 + self.Natom ) + 9

            try :
                data = np.loadtxt(self.trj_file, dtype=np.double, skiprows=skip_line, usecols=(1,2,3,4,5,6), max_rows=self.Natom)
            except :
                flag = 0

            if (flag == 1 ) :

                if ( (Nstep < Nstep_guess) ) :

                    cord[Nstep, :, :] = data[:, :3]
                    force[Nstep, :, :] = data[:, 3:]

                    Nstep += 1
                    #print("Nstep = ", Nstep)

                else :

                    UptoNstep = Nstep_guess

                    cord_cpy = np.copy( cord )
                    force_cpy = np.copy( force )

                    Nstep_guess += 40

                    cord = np.zeros([Nstep_guess, self.Natom, 3], dtype=np.double)
                    force = np.zeros([Nstep_guess, self.Natom, 3], dtype=np.double)

                    cord[:UptoNstep, :, :] = cord_cpy
                    force[:UptoNstep, :, :] = force_cpy

                    cord[Nstep, :, :] = data[:, :3]
                    force[Nstep, :, :] = data[:, 3:]

                    Nstep += 1

                    #print("Nstep = ", Nstep)

                    del( cord_cpy )
                    del( force_cpy )

        print()
        print("\tNumber of MD snapshots available: ", Nstep)

        cord_cpy = np.copy( cord )
        force_cpy = np.copy( force )

        cord = np.zeros([Nstep, self.Natom, 3], dtype=np.double)
        force = np.zeros([Nstep, self.Natom, 3], dtype=np.double)
        cord[:, :, :] = cord_cpy[:Nstep, :, :]
        force[:, :, :] = force_cpy[:Nstep, :, :]

        del( cord_cpy )
        del( force_cpy )

        return Nstep, cord, force


    def ReadRecord5(self) :

        cell_basis_record_dset = 'cell_basis'
        frac_cord_dset = 'frac_coord'
        cart_cord_dset = 'cart_coord'

        print("\t*** Reading record of Super-Cell information '", self.recordfile, "'***")
        h5f = h5py.File(self.recordfile, mode='r')

        cell_basis = h5f[cell_basis_record_dset][:]
        cart_cord = h5f[cart_cord_dset][:]
        frac_cord = h5f[frac_cord_dset][:]

        h5f.close()

        return cell_basis, cart_cord, frac_cord


    def ReadCord(self) :

        num_step = 0
        print()
        for ts in range(self.Nstep) :

            cord_ts = self.cord[ts, :, :]

            dL = (cord_ts - self.cart_cord)

            for atm in range(self.Natom) :

                disp_xyz = np.copy(dL[atm, :])

                dsp_xyz = self.CheckBound( atm, disp_xyz, np.copy(cord_ts[atm, :]) )

                nx, ny, nz, mu = self.cell_basis[atm, :]
                self.dU[num_step, nx, ny, nz, mu, :] = dsp_xyz #np.copy(dL[atm, :])

            num_step += 1
            print("\tRead Coordinate: No. of MD snapshot covered:", num_step, "/", self.Nstep, end='\r', flush=True)


    #-# def CheckBound(self, atm, dsp_xyz, crd_xyz) :

    #-#     dsp_abc = np.matmul(dsp_xyz, self.lat_unit_inv)

    #-#     if np.any( abs(dsp_abc) >= (self.lat_len / 2.0) ) :
    #-#         print("*************** Atom moves out of box ***********************")
    #-#         print("PBC Wrap needed")
    #-#         lat_vec = np.copy(self.u.lattice)

    #-#         for xx in range(3) :

    #-#             if ( abs(dsp_abc[xx]) >= (self.lat_len[xx] / 2.0) ) :
    #-#                 if dsp_abc[xx] < 0.0 :
    #-#                     crd_xyz += lat_vec[xx, :]
    #-#                 elif dsp_abc[xx] > 0.0 :
    #-#                     crd_xyz -= lat_vec[xx, :]

    #-#         cord_xyz0 = np.copy(self.cart_cord[atm, :])
    #-#         dsp_xyz_new = crd_xyz - cord_xyz0

    #-#         return dsp_xyz_new

    #-#     else :
    #-#         return dsp_xyz

    def CheckBound(self, atm, dsp_xyz, crd_xyz, ts, InitVaspData) :

        #print( self.lat_unit_inv )
        dsp_abc = np.matmul(dsp_xyz, self.lat_unit_inv)

        if np.any( abs(dsp_abc) >= (self.lat_len / 2.0) ) :
            print("*************** Atom moves out of box ***********************")
            print("PBC Wrap needed")
            ##print( self.cord[ts, :, :] )
            #print( ts, atm )
            #print( InitVaspData[atm, :], "--->", self.cord[ts, atm, :], "--->", dsp_abc )
            #print( )
            #lat_vec = np.copy(self.u.lattice)

            #print( lat_vec )

            for xx in range(3) :

                if ( abs(dsp_abc[xx]) >= (self.lat_len[xx] / 2.0) ) :
                    #print( xx )
                    if dsp_abc[xx] < 0.0 :
                        crd_xyz += self.lat_vec[xx, :]
                    elif dsp_abc[xx] > 0.0 :
                        crd_xyz -= self.lat_vec[xx, :]

            cord_xyz0 = np.copy(InitVaspData[atm, :])
            dsp_xyz_new = crd_xyz - cord_xyz0
            #print( cord_xyz0, "--->", crd_xyz, "--->", dsp_xyz_new )
            #print( dsp_xyz_new )

            return dsp_xyz_new

        else :
            return dsp_xyz


    def InitForce(self) :

        force_init = np.zeros([self.Natom, 3], dtype=np.double)
        pos_init = np.zeros([self.Natom, 3], dtype=np.double)

        skip_line = 9
        try :

            data = np.loadtxt(self.init_frc_file, dtype=np.double, skiprows=skip_line, usecols=(1,2,3,4,5,6), max_rows=self.Natom)
            pos_init = data[:, :3]
            force_init = data[:, 3:]

        except :

            print("\tProblem in reading lammps trajectory file for initial force: '", self.init_frc_file, " '")
            print("\tInitial force of unperterbed system is set to zero")

        return pos_init, force_init


    def ReadForce(self) :

        num_step = 0
        print()
        for ts in range(self.Nstep) :

            force_ts = self.force[ts, :, :]

            force = (force_ts-self.initForc)

            for atm in range(self.Natom) :
                nx, ny, nz, mu = self.cell_basis[atm, :]
                self.forc[num_step, nx, ny, nz, mu, :] = np.copy(force[atm, :])

            num_step += 1
            print("\tRead Force: No. of MD snapshot covered:", num_step, "/", self.Nstep, end='\r', flush=True)

    
    def Writeh5(self, T) :

        filename = 'disp_forc_'+self.u.prefix+str(T)+'K.h5'

        print()
        print("\t*** Writing force-displacement dataset in file:'", filename, "'***")

        disp_dset = 'disp_dataset'
        force_dset = 'force_dataset'

        h5f = h5py.File(filename, mode='w')

        h5f.create_dataset(disp_dset, data=self.dU)
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

    # Run: python3 ReadForcDisp.py -inp input_Config1.nml -T 300 -DFfile out.log -MDEfile SiC.MDE -ts 5000 -Taccu 15 -ts_cut 3500 -ts_skip 15
    parser.add_argument("-inp", type=str, default='input.nml', help="Enter the name of input file")
    parser.add_argument("-T", type=float, default=300, help="Temperature in Kelvin ")
    parser.add_argument("-MD", type=str, default="siesta", help="Provide the MD interface (only siesta and lammps available)")

    parser.add_argument("-ts", type=int, default=5000, help="Total number of timesteps in MD run (Effective for siesta only interface)")
    parser.add_argument("-DFfile", type=str, default='out.log', help="Displacement-Force output file of MD run (Effective only for siesta interface)")
    parser.add_argument("-MDEfile", type=str, default='NaCl.MDE', help="MD run thermodynamic information (Effective only for siesta interface)")
    parser.add_argument("-Taccu", type=float, default=10, help="Maximum shift from mean temperature (in K)")
    parser.add_argument("-ts_cut", type=int, default=4000, help="Trajectory to be read after this many timesteps (Effective only for siesta interface)")

    parser.add_argument("-ts_skip", type=int, default=15, help="Minimum timestep difference between 2 trajectory (Effective only for siesta interface)")
    parser.add_argument("--init_force_zero", dest='init_force_flag', action='store_true', \
                        help="Whether you want force on the atoms of starting configuration to be zero")

    parser.add_argument("-lammps_initforce", type=str, default='default', help="Trajectory containg intial for of the unperturbed system (Effective only for lammps interface)")
    parser.add_argument("-lammps_trj", type=str, default='default', help="Trajectory filename from lammps (Effective only for lammps interface)")

    parser.set_defaults(init_force_flag=False)

    args = parser.parse_args()
    nml_fname = args.inp
    MD = args.MD

    T = args.T
    siestafile = args.DFfile
    mde_file = args.MDEfile
    total_ts = args.ts
    Taccu = args.Taccu
    ts_cut = args.ts_cut
    t_skip = args.ts_skip
    ZeroIntialTF = args.init_force_flag

    lammps_init_force_file = args.lammps_initforce
    lammps_trj_file = args.lammps_trj

    #............::::::::::::::.............#
    Tmean = T
    #Taccu = 10
    #ts_cut = 4000 #ts - 1000
    #t_skip = 15
    #............::::::::::::::.............#

    if ( MD == 'siesta' ) :

        data = ReadSIESTA(nml_fname, siestafile, mde_file, total_ts, \
                          Tmean, ts_cut, dT=Taccu, d_time=t_skip, init_force_zero=ZeroIntialTF)
        data.ReadCord()
        data.ReadForce()
        data.Writeh5(T)

    elif ( MD == 'lammps' ) :

        data = Read_lammps(nml_fname, lammps_init_force_file, lammps_trj_file)
        data.ReadCord()
        data.ReadForce()
        data.Writeh5(T)

    else :

        print("\t *** ERROR: unrecognized MD interface", MD, "  (only siesta and lammps are available) ***")

    print()

    #.........::::::::::::::: For Debug purpose :::::::::::::::............#
    debug = False
    if ( debug ) :

        num = 50

        if ( MD == 'siesta' ) :
            lookup_word = "outcoor"
            line_cord = data.FindLineNumber(lookup_word)
            lookup_word = "(eV/Ang):"
            line_frc = data.FindLineNumber(lookup_word)

            tsteps = data.timestep
            ts = tsteps[num]
            #ln = line_frc[ts-1]
            ln = line_cord[ts-1]

            print("Line number in the file: ", ln)
            print()

        frc, cord = data.Readh5(T)

        print("Data: ") 
        num_atm = 0
        for nx in range(data.N[0]) :
            for ny in range(data.N[1]) :
                for nz in range(data.N[0]) :
                    for bb in range(len(data.u.basis_frac)) :

                        #print(frc[num, nx, ny, nz, bb, :])
                        print(cord[num, nx, ny, nz, bb, :] + \
                              data.cart_cord[num_atm, :])
                        num_atm += 1
        print()
        print('cord[16, 0, 1, 2, 0, 0]:', cord[16, 0, 1, 2, 0, 0])
        print('frc[16, 0, 1, 2, 0, 0]:', frc[16, 0, 1, 2, 0, 0])
    #.........::::::::::::::: For Debug purpose :::::::::::::::............#

