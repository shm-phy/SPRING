
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

    use hdf5,                       only : HID_T, H5F_ACC_RDONLY_F, H5P_DEFAULT_F, h5fopen_f, &
                                         & h5fclose_f

    use kinds,                      only : dp
    use parse_cmd_line,             only : get_command_line, ShowWelcomeBanner
    use unit_cell,                  only : cell
    use FC2_mod,                    only : FC2type
    use EwaldMod,                   only : EwaldParam
    use CreateMesh,                 only : mesh_points, q_points_highsymm
    use Irr_q_point,                only : q_points_data
    use phonon_m,                   only : Phon
    use ThermalStochasticSnap,      only : GenerateInitialPos, CreateSnapshotType1, &
                                         & CreateSnapshotType2
    use hdf5_wrap,                  only : wCellBasisDisp

    implicit none

    type(cell)                                      :: sys
    type(FC2type)                                   :: FC2

    type(Phon)                                      :: phonon_dat
    type(EwaldParam)                                :: EwaldConst
    logical                                         :: LongEW

    character(len=256)                              :: filename, outfile
    real(dp)                                        :: T, Tsnap
    character(len=24)                               :: Temp_char
    character(len=8)                                :: frmt

    character(len=64)                               :: fc2file

    integer(HID_T)                                  :: file_id
    integer                                         :: hdferr

    integer, dimension(3)                           :: mesh, shift
    type(q_points_data)                             :: Qpoints

    integer                                         :: num_irr_q

    real(dp), dimension(:,:), allocatable           :: InitPos
    real(dp), dimension(:,:), allocatable           :: R_Ncell
    integer, dimension(:,:), allocatable            :: cell_basis_record

    integer                                         :: Nsnap, seed
    real(dp), dimension(:,:,:,:,:,:), allocatable   :: dU

    integer                                         :: typ
    character(len=8)                                :: calculator
    logical                                         :: abinit_xred

    
    ! ------------------------------------------------------------------------------------------- !
    call ShowWelcomeBanner()

    call get_command_line(filename, T, Tsnap, Nsnap, seed, typ, calculator, abinit_xred)
    call sys%init_data(filename)

    frmt = '(F6.1)'
    write(Temp_char, frmt) T

    fc2file = 'FC2ndRenorm_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
    call h5fopen_f( fc2file, H5F_ACC_RDONLY_F, file_id, hdferr, &  
                  & H5P_DEFAULT_F )

        CheckIfFileExists: if ( hdferr < 0 ) then

            write(*, *)
            write(*, 10)
            write(*, 54) 
            write(*, 55) fc2file
            write(*, 20)
            write(*, 10)
            write(*, *)

            fc2file = 'FC2nd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
            call h5fopen_f( fc2file, H5F_ACC_RDONLY_F, file_id, hdferr, &  
                          & H5P_DEFAULT_F )

            NormalFileExists: if ( hdferr < 0 ) then

                write(*, 56) 
                write(*, 57) fc2file
                56 FORMAT( 'Trying to read unrenormalized 2nd-Order IFCs.' )
                57 FORMAT( 'FAILURE: file does not exists: ', A64, '.  STOPPING ...')
                STOP

            end if NormalFileExists

        end if CheckIfFileExists

    call h5fclose_f(file_id, hdferr)

    call FC2%set_FC2(fc2file)

    !mesh = (/5, 5, 5/)
    mesh = sys%sup_cell
    shift(:) = 0

    call Qpoints%make_tetrahedron(sys, mesh, shift)

    num_irr_q = size( Qpoints%irr_q_int, 2 )

    write(*, 12) num_irr_q

    12 FORMAT('Number of irreducible q-points: ', I5)

    ! ======================= find phonon freq., velocity and eigenvectors ====================== !

    LongEW = sys%Force_dd

    Ew_chk: if ( LongEW ) then

        call EwaldConst%set_EwaldParam(sys, filename)

        write(*, *)
        write(*, 100) EwaldConst%Lmb

        write(*, 200) EwaldConst%Gmesh
        write(*, *)

        100 FORMAT('The alpha parameter for Ewald Sum: ', ES10.3)
        200 FORMAT('The NKcut for Ewald Sum:', 3I3)

    end if Ew_chk

    call phonon_dat%set_Phon(sys, FC2, Qpoints%irr_q, EwaldConst, Tsnap, LongEW)

    ! ======================= find phonon freq., velocity and eigenvectors ====================== !

    call GenerateInitialPos(sys, InitPos, R_Ncell, cell_basis_record)

    typ_chk: if ( typ == 1 ) then

        call CreateSnapshotType1(sys, Tsnap, Nsnap, seed, mesh, Qpoints, phonon_dat, &
                               & InitPos, R_Ncell, cell_basis_record, calculator, abinit_xred, dU)

    else typ_chk

        call CreateSnapshotType2(sys, Tsnap, Nsnap, seed, mesh, Qpoints, phonon_dat, &
                               & InitPos, R_Ncell, cell_basis_record, calculator, abinit_xred, dU)

    end if typ_chk

    write(Temp_char, frmt) Tsnap
    outfile = 'disp_forc_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

    call wCellBasisDisp(cell_basis_record, InitPos, dU, outfile)

    write(*, '(/)')

    10 FORMAT("============----------------------------------------------------------============")
    54 FORMAT( 'Trying to read renormalized 2nd-Order IFCs.') 
    55 FORMAT( 'FAILURE: file does not exists: ', A64)
    20 FORMAT("Trying to read unrenormalized IFC2 (IGNORE THE ABOVE MESSAGES FROM HDF5)")

end program main

