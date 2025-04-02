
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
    use FC4_data_mod,               only : FC4_dat, RedInFCsp
    use hdf5_wrap,                  only : wFC4_Reconst
    use ReconstFC4,                 only : IterOver_N2_N3_N4

    implicit none

    type(cell)                                                          :: sys
    character(len=256)                                                  :: filename
    real(dp)                                                            :: T

    type(FC4_dat)                                                       :: FC4
    character(len=128)                                                  :: fc4file, &
                                                                         & fc4indfc_file
    character(len=24)                                                   :: Temp_char
    character(len=8)                                                    :: frmt

    real(dp), allocatable, dimension(:,:,:,:,:,:,:,:)                   :: FC4_Reconst

    integer                                                             :: Nbasis, Natm_Max

    character(len=64)                                                   :: outfile


    call ShowWelcomeBanner()

    call get_command_line(filename, T)
    call sys%init_data(filename)

    frmt = '(F6.1)'
    write(Temp_char, frmt) T

    fc4file = 'FC_4th_common_'//trim(sys%prefix)//'_F.h5'
    fc4indfc_file = 'FD_4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

    call FC4%set_FC4dat(fc4file, fc4indfc_file)

    Nbasis = sys%natm
    Natm_Max =  maxval( FC4%atmNum )

    allocate( FC4_Reconst(3,3,3,3,Natm_Max,Natm_Max,Natm_Max,Nbasis) )
    FC4_Reconst = 0.0_dp

    write(*, 100)

    call IterOver_N2_N3_N4(FC4, Nbasis, FC4_Reconst)

    write(*, '(/)')

    outfile = 'FC4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
    call wFC4_Reconst(outfile, FC4_Reconst, FC4%atmNum, FC4%cb_Indx)

    100 FORMAT('Reconstructing 4-th order IFCs from Independent IFCs ... ')

end program main

