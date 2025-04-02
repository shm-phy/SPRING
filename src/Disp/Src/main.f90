
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
    use EwaldMod,       only : EwaldParam
    use hdf5_wrap,      only : r2ndFC, w2dh5_no_grp
    use CreateMesh,     only : q_points_highsymm
    use DynaMat,        only : get_phonon_disp

    implicit none

    type(cell)                  :: sys
    character(len=256)          :: filename
    real(dp)                    :: T
    character(len=24)           :: Temp_char
    character(len=8)            :: frmt

    type(EwaldParam)            :: EwaldConst

    character(len=50)           :: fc2file
    integer, allocatable        :: atmNum(:)
    integer, allocatable        :: atmIndx(:, :, :)
    real(dp), allocatable       :: FC2(:, :, :, :, :)


    logical                     :: PhononDisp, LongEW
    integer                     :: num_points 
    namelist        /PhononDispInfo/     PhononDisp, LongEW, num_points

    real(dp), allocatable       :: q_high_sym(:, :)
    namelist        /HighSymPath/     q_high_sym

    real(dp), allocatable       :: qpoints(:, :)

    real(dp), allocatable       :: freq_dat(:, :)
    character(len=50)           :: outfile

    integer                     :: err
    character(len=512)          :: err_msg

    call ShowWelcomeBanner()

    call get_command_line(filename, T)
    call sys%init_data(filename)

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

    close(unit=5, status='KEEP', iostat=err, iomsg=err_msg)

    Ew_chk: if ( LongEW ) then

        call EwaldConst%InitEwald(sys, filename)
        write(*, 100) EwaldConst%Lmb
        100 FORMAT('The parameter alpha for Ewald Sum: ', ES10.3)
        write(*, 200) EwaldConst%Gmesh
        200 FORMAT('The NKcut for Ewald Sum:', 3I3)

    end if Ew_chk

    frmt = '(F6.1)'
    write(Temp_char, frmt) T

    fc2file = 'FC2nd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
    call r2ndFC(fc2file, atmNum, atmIndx, FC2)

    call q_points_highsymm(q_high_sym, sys%G, qpoints)
    call get_phonon_disp(sys, atmNum, atmIndx, FC2, &
                       & qpoints, EwaldConst, LongEw, freq_dat)

    outfile = 'dispersion_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
    call w2dh5_no_grp(freq_dat, outfile)

    write(*, '(/)')

end program main

