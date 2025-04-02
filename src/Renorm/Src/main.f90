
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

    use kinds,          only : dp
    use parse_cmd_line, only : get_command_line, ShowWelcomeBanner
    use unit_cell,      only : cell
    use FC2_mod,        only : FC2type
    use FC4_mod,        only : FC4type
    use EwaldMod,       only : EwaldParam
    use hdf5_wrap,      only : w2dh5_no_grp, wFC2_Renorm
    use CreateMesh,     only : mesh_points, q_points_highsymm
    use DynaMat,        only : SetDynEw_allq, get_phonon_disp
    use Renorm,         only : Renormalize
    use FindForce,      only : CalculateForce

    implicit none

    type(cell)                                  :: sys
    type(FC2type)                               :: FC2_in
    type(FC2type)                               :: FC2_out
    type(FC4type)                               :: FC4
    type(EwaldParam)                            :: EwaldConst

    character(len=256)                          :: filename
    real(dp)                                    :: T
    character(len=24)                           :: Temp_char
    character(len=8)                            :: frmt

    character(len=64)                           :: fc2file, fc4file

    integer, dimension(3)                       :: Qmesh
    logical                                     :: shift
    real(dp)                                    :: accu
    real(dp), allocatable                       :: q_points(:, :)
    real(dp), dimension(3)                      :: qshift
    namelist        /RenormInfo/   Qmesh, shift, accu

    character(len=64)                           :: fdfile, renorm_out_file
    real(dp), dimension(:), allocatable         :: force_renorm, force_md

    logical                                     :: PhononDisp, LongEW
    integer                                     :: num_points 
    namelist        /PhononDispInfo/     PhononDisp, LongEW, num_points

    real(dp), allocatable                       :: q_high_sym(:, :)
    namelist        /HighSymPath/     q_high_sym
    real(dp), allocatable                       :: qpoints_hs(:, :)
    real(dp), allocatable                       :: freq_dat_hs(:, :)
    character(len=50)                           :: outfile

    integer                                     :: err, kk
    character(len=512)                          :: err_msg

    call ShowWelcomeBanner()

    call get_command_line(filename, T)
    call sys%init_data(filename)

    frmt = '(F6.1)'
    write(Temp_char, frmt) T

    fc2file = 'FC2nd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
    call FC2_in%set_FC2(fc2file)

    fc4file = 'FC4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
    call FC4%set_FC4(fc4file, FC2_in%atmIndx)

    open(unit=5, file=filename, status='OLD', iostat=err, iomsg=err_msg, &
         action='READ', delim='APOSTROPHE')

        open_chk: if ( err /= 0 ) then
            write(*, 300) err
            300 FORMAT('Input file OPEN failed: iostat =', I3) 
            write(*, 400) err_msg
            400 FORMAT('Error message = ', A)
        end if open_chk

        read(unit=5, nml=PhononDispInfo)

        disp_chk: if ( PhononDisp ) then

            allocate(q_high_sym(4, num_points))
            read(unit=5, nml=HighSymPath)

        end if disp_chk

        read(unit=5, nml=RenormInfo)

    close(unit=5, status='KEEP', iostat=err, iomsg=err_msg)

    call mesh_points(cell_vec = sys%G, mesh = Qmesh, div_cell = .true., mesh_cord = q_points)

    shift_chk: if ( shift ) then

        qshift = ( (sys%G(:, 1) / (Qmesh(1) *2.0_dp)) + &
                  &(sys%G(:, 2) / (Qmesh(2) *2.0_dp)) + &
                  &(sys%G(:, 3) / (Qmesh(3) *2.0_dp)) )

        loop_q: do kk = 1, size(q_points, 2)
            q_points(:, kk) = q_points(:, kk) + qshift
        end do loop_q

    end if shift_chk

    Ew_chk: if ( LongEW ) then

        call EwaldConst%set_EwaldParam(sys, filename)
        write(*, 100) EwaldConst%Lmb
        100 FORMAT('The alpha parameter for Ewald Sum: ', ES10.3)
        write(*, 200) EwaldConst%Gmesh
        200 FORMAT('The NKcut for Ewald Sum:', 3I3)

        call SetDynEw_allq(sys, q_points, EwaldConst)

    end if Ew_chk

    call Renormalize(sys, T, accu, EwaldConst, &
                    &LongEw, FC2_in, FC4, q_points, FC2_out)

    fdfile = 'disp_forc_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
    call CalculateForce(fdfile, FC2_out, force_renorm, force_md)

    renorm_out_file = 'FC2ndRenorm_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
    call wFC2_Renorm(renorm_out_file, FC2_out%fC, FC2_out%atmNum, FC2_out%atmIndx, &
                   & force_renorm, force_md)

    disp: if ( PhononDisp ) then

        call q_points_highsymm(q_high_sym, sys%G, qpoints_hs)

        if ( LongEW ) call SetDynEw_allq(sys, qpoints_hs(2:4, :), EwaldConst)

        call get_phonon_disp(sys, FC2_out, qpoints_hs, EwaldConst, LongEw, freq_dat_hs)

        outfile = 'dispersion_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
        call w2dh5_no_grp(freq_dat_hs, outfile)

    end if disp

    write(*, '(/)')

    deallocate(q_high_sym, qpoints_hs, freq_dat_hs, q_points)

end program main

