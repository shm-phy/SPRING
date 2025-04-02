
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
    use hdf5_wrap,                  only : wFC2_Reconst
    use ReconstFC2,                 only : IterOver_N2

    implicit none

    type(cell)                                                          :: sys
    character(len=256)                                                  :: filename
    real(dp)                                                            :: T

    type(FC2_dat)                                                       :: FC2
    character(len=128)                                                  :: fc2file, &
                                                                         & fc2indfc_file
    character(len=24)                                                   :: Temp_char
    character(len=8)                                                    :: frmt

    real(dp), allocatable, dimension(:,:,:,:)                           :: FC2_Reconst

    integer                                                             :: Nbasis, Natm_Max

    character(len=64)                                                   :: outfile


    call ShowWelcomeBanner()

    call get_command_line(filename, T)
    call sys%init_data(filename)

    frmt = '(F6.1)'
    write(Temp_char, frmt) T

    fc2file = 'FC_2nd_common_'//trim(sys%prefix)//'_F.h5'
    fc2indfc_file = 'FD_2nd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

    call FC2%set_FC2dat(fc2file, fc2indfc_file)

    Nbasis = sys%natm
    Natm_Max =  maxval(FC2%atmNum)

    allocate( FC2_Reconst(3,3,Natm_Max,Nbasis) )
    FC2_Reconst = 0.0_dp

    write(*, 100)

    call IterOver_N2(FC2, Nbasis, FC2_Reconst)

    write(*, '(/)')

    outfile = 'FC2nd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
    call wFC2_Reconst(outfile, FC2_Reconst, FC2%atmNum, FC2%cb_Indx)

    100 FORMAT('Reconstructing 2-nd order IFCs from Independent IFCs ... ')

end program main

