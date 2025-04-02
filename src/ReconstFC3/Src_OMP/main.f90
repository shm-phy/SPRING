
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
    use FC3_data_mod,               only : FC3_dat, RedInFCsp
    use hdf5_wrap,                  only : wFC3_Reconst
    use ReconstFC3,                 only : IterOver_N2_N3

    implicit none

    type(cell)                                                          :: sys
    character(len=256)                                                  :: filename
    real(dp)                                                            :: T

    type(FC3_dat)                                                       :: FC3
    character(len=128)                                                  :: fc3file, &
                                                                         & fc3indfc_file
    character(len=24)                                                   :: Temp_char
    character(len=8)                                                    :: frmt

    real(dp), allocatable, dimension(:,:,:,:,:,:)                       :: FC3_Reconst

    integer                                                             :: Nbasis, Natm_Max

    character(len=64)                                                   :: outfile


    call ShowWelcomeBanner()

    call get_command_line(filename, T)
    call sys%init_data(filename)

    frmt = '(F6.1)'
    write(Temp_char, frmt) T

    fc3file = 'FC_3rd_common_'//trim(sys%prefix)//'_F.h5'
    fc3indfc_file = 'FD_3rd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

    call FC3%set_FC3dat(fc3file, fc3indfc_file)

    Nbasis = sys%natm
    Natm_Max =  maxval( FC3%atmNum )

    allocate( FC3_Reconst(3,3,3,Natm_Max,Natm_Max,Nbasis) )
    FC3_Reconst = 0.0_dp

    write(*, 100)

    call IterOver_N2_N3(FC3, Nbasis, FC3_Reconst)

    write(*, '(/)')

    outfile = 'FC3rd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
    call wFC3_Reconst(outfile, FC3_Reconst, FC3%atmNum, FC3%cb_Indx)

    100 FORMAT('Reconstructing 3-rd order IFCs from Independent IFCs ... ')

end program main

