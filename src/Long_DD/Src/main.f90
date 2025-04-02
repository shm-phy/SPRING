
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

    use kinds,              only : dp
    use parse_cmd_line,     only : get_command_line, ShowWelcomeBanner
    use unit_cell,          only : cell
    use EwaldSum,           only : DoEwaldSum, Inverse_FT
    use FindForce,          only : CalculateForce
    use hdf5_wrap,          only : wFC2_dd!, r_force_disp

    implicit none

    type(cell)                                      :: sys
    character(len=256)                              :: filename 
    real(dp)                                        :: T
    character(len=24)                               :: Temp_char
    character(len=8)                                :: frmt

    integer, dimension(3)                           :: qmesh, Gmesh, Rmesh
    logical                                         :: decide_EwParam=.true.
    real(dp)                                        :: Lmb = 1.01_dp, prec_Ew=1.0E-6_dp
    namelist        /EwaldInfo/     qmesh, Gmesh, Rmesh, Lmb, decide_EwParam, prec_Ew

    real(dp), allocatable                           :: qpoints(:, :)
    complex(dp), allocatable                        :: CEw(:, :, :, :, :)

    complex(dp), dimension(:,:,:,:), allocatable    :: IFC
    real(dp), dimension(:,:,:,:), allocatable       :: FC2_dd
    integer, dimension(:,:,:), allocatable          :: RCb

    character(len=128)                              :: fdfile
    real(dp), dimension(:), allocatable             :: force_dd, force_md

    integer                                         :: err
    character(len=512)                              :: err_msg

    character(len=128)                              :: outfile
    

    call ShowWelcomeBanner()

    call get_command_line(filename, T)
    call sys%init_data(filename)

    open(unit=5, file=filename, status='OLD', iostat=err, iomsg=err_msg, &
         action='READ', delim='APOSTROPHE')

        open_chk: if ( err /= 0 ) then
            write(*, 100) err
            100 FORMAT('Input file OPEN failed: iostat = ', I3)
            write(*, 200) err_msg
            200 FORMAT('Error message = ', A128)
        end if open_chk

        read(unit=5, nml=EwaldInfo)

    close(unit=5, status='KEEP', iostat=err, iomsg=err_msg)

    qmesh = sys%sup_cell

    call DoEwaldSum(sys, qmesh, Rmesh, Gmesh, Lmb, decide_EwParam, prec_Ew, qpoints, CEw)
    
    call Inverse_FT(sys, CEw, qpoints, IFC, RCb)

    allocate( FC2_dd(3, 3, product(sys%sup_cell) * sys%natm, sys%natm) )
    FC2_dd = dble(IFC)

    write(*, 300) maxval(dabs(dimag(IFC)))
    300 FORMAT('Maximum value of the imaginary part (ideally 0) of the IFC: ', E14.6)

    frmt = '(F6.1)'
    write(Temp_char, frmt) T

    fdfile = 'disp_forc_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
    call CalculateForce(fdfile, FC2_dd, RCb, force_dd, force_md)

    outfile = 'FC2_dd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
    call wFC2_dd(FC2_dd, RCb, force_dd, force_md, outfile)

    write(*, '(/)')

    deallocate(IFC)

end program main

