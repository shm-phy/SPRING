
! ============================================================================= !
!    Copyright (C) 2022  Soham Mandal                                           !
!                                                                               !
!    This program is free software: you can redistribute it and/or modify       !
!    it under the terms of the GNU General Public License as published by       !
!    the Free Software Foundation, either version 3 of the License, or          !
!    (at your option) any later version.                                        !
!                                                                               !
!    This program is distributed in the hope that it will be useful,            !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of             !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              !
!    GNU General Public License for more details.                               !
!                                                                               !
!    You should have received a copy of the GNU General Public License          !
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.     !
!                                                                               !
!    e-mail: phy.soham@gmail.com                                                !
! ============================================================================= !

program main

    use kinds,                      only : dp
    use parse_cmd_line,             only : get_command_line, ShowWelcomeBanner
    use unit_cell,                  only : cell
    use FC2_data_mod,               only : FC2_dat, RedInFCsp
    use Helper,                     only : ReadPointDefVar
    use hdf5_wrap,                  only : r_force_disp, wFD2_part, wFD_JoinAll
    use DispMatFC2,                 only : IterOver_N2, find_force

    implicit none

    type(cell)                                              :: sys
    type(FC2_dat)                                           :: FC2

    character(len=256)                                      :: filename
    character(len=128)                                      :: fc2file
    character(len=24)                                       :: Temp_char
    character(len=8)                                        :: frmt
    character(len=128)                                      :: fdfile
    character(len=64)                                       :: outfile, file_full

    real(dp), dimension(:,:,:,:,:,:), allocatable           :: disp
    real(dp), dimension(:,:,:,:,:,:), allocatable           :: force
    real(dp), allocatable, dimension(:,:)                   :: Mat
    real(dp), allocatable, dimension(:)                     :: force_mua, &
                                                             & Mass_def !** Variable for point-defect calculation **!
    real(dp)                                                :: T

    integer, allocatable, dimension(:)                      :: flags !** Variable for point-defect calculation **!
    integer, allocatable, dimension(:, :)                   :: cell_atm_def !** Variable for point-defect calculation **!
    integer                                                 :: mu, alpha
    integer                                                 :: NIndFC, Nbasis, &
                                                             & nstep, num_row, &
                                                             & NumDef_sup !** Variable for point-defect calculation **!
    integer, dimension(3)                                   :: sup_dim, &
                                                             & qMesh_pd !** Variable for point-defect calculation **!

    logical                                                 :: PointDef !** Variable for point-defect calculation **!



    call ShowWelcomeBanner()

    call get_command_line(filename, T, mu, alpha)
    call sys%init_data(filename)
    Nbasis = sys%natm !** Point defect related **!

    fc2file = 'FC_2nd_common_'//trim(sys%prefix)//'_F.h5'
    call FC2%set_FC2dat(fc2file)

    ! =================== Read Point-Defect realted input from namelist file ==================== !

    allocate( Mass_def(Nbasis) )
    call ReadPointDefVar( filename, Nbasis, sys, qMesh_pd, Mass_def, PointDef, NumDef_sup, cell_atm_def )

    ! =================== Read Point-Defect realted input from namelist file ==================== !

    frmt = '(F6.1)'
    write(Temp_char, frmt) T
    fdfile = 'disp_forc_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
    call r_force_disp(fdfile, disp, force)

    outfile = 'FDpart_2nd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

#ifdef _USE_DAXPY
    write(*, 75)
    75 FORMAT( 4X, 4('-'), '(dev) Computations to be done with daxpy', 4('-') )
#else
    write(*, 76)
    76 FORMAT( 4X, 4('-'), '(dev) Computations to be done with out daxpy', 4('-') )
#endif

    nstep = size(disp, 1)
    sup_dim = sys%sup_cell
    num_row = nstep*product(sup_dim)
    NIndFC = FC2%ind_fc

    allocate( Mat(num_row, NIndFC) )
    allocate( flags(num_row) )

    write(*, *)
    mu_loop: do mu = 1, Nbasis
        alpha_loop: do alpha = 1, 3

            write(*, 100) mu, alpha
            Mat = 0.0_dp
            flags = 0

            call IterOver_N2( mu, alpha, num_row, nstep, NIndFC, NumDef_sup, PointDef, &
                            & cell_atm_def, sup_dim, FC2, disp, flags, Mat )

            call find_force(sup_dim, nstep, num_row, mu, alpha, force, force_mua)

            call wFD2_part(Mat, force_mua, flags, outfile, mu, alpha)

            write(*, *)

            last_chk: if ( (mu == sys%natm) .and. (alpha == 3) ) then

                file_full = 'FD_2nd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
                call wFD_JoinAll(sys%natm, num_row, FC2%ind_fc, outfile, file_full)

            end if last_chk

            deallocate( force_mua )

        end do alpha_loop
    end do mu_loop

    100 FORMAT("Starting calculation of 2nd-order displacement matrix with mu = ", I3, " and alpha = ", I2, " ...")

end program main

