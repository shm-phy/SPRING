
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
    use parse_cmd_line,             only : get_command_line
    use unit_cell,                  only : cell
    use FC4_data_mod,               only : FC4_dat, RedInFCsp
    use hdf5_wrap,                  only : r_force_disp, wFD4_part, wFD_JoinAll
    use DispMatFC4,                 only : IterOver_N2_N3_N4_mp2, find_force

    implicit none

    type(cell)                                              :: sys
    character(len=256)                                      :: filename
    real(dp)                                                :: T

    type(FC4_dat)                                           :: FC4
    character(len=128)                                      :: fc4file

    character(len=24)                                       :: Temp_char
    character(len=8)                                        :: frmt
    character(len=128)                                      :: fdfile
    real(dp), dimension(:,:,:,:,:,:), allocatable           :: disp
    real(dp), dimension(:,:,:,:,:,:), allocatable           :: force

    integer                                                 :: mu, alpha
    real(dp), allocatable, dimension(:,:), codimension[:]   :: Mat
    real(dp), allocatable, dimension(:)                     :: force_mua

    integer                                                 :: Natm_mu, NFC_tot, &
                                                             & NFC_len, my_len, nstep, num_row
    integer, dimension(2)                                   :: my_start_end
    integer, dimension(3)                                   :: sup_dim
    integer                                                 :: numfc, N2, N3, N4, &
                                                             & beta, gama, delta
    integer, allocatable, dimension(:,:)                    :: all_FC_indx, my_FC_indx

    character(len=64)                                       :: outfile, file_full



    if ( this_image() == 1) write(*, '(/)')

    call get_command_line(filename, T, mu, alpha)
    call sys%init_data(filename)

    fc4file = 'FC_4th_common_'//trim(sys%prefix)//'_F.h5'
    call FC4%set_FC4dat(fc4file)


    frmt = '(F6.1)'
    write(Temp_char, frmt) T
    fdfile = 'disp_forc_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
    call r_force_disp(fdfile, disp, force)

    Natm_mu = FC4%atmNum(mu)
    NFC_tot = (Natm_mu**3) * (3**3)
    NFC_len = NFC_tot / num_images()
    nstep = size(disp, 1)
    sup_dim = sys%sup_cell
    num_row = nstep*product(sup_dim)

#ifdef _USE_DAXPY
    if ( this_image() == 1 ) write(*, 75)
    75 FORMAT( 4X, 4('-'), '(dev) Computations to be done with daxpy', 4('-') )
#else
    if ( this_image() == 1 ) write(*, 76)
    76 FORMAT( 4X, 4('-'), '(dev) Computations to be done with out daxpy', 4('-') )
#endif

    my_start_end(1) = ((this_image() - 1) * NFC_len) + 1
    my_start_end(2) = my_start_end(1) + NFC_len - 1

    last_img: if ( (this_image() == num_images()) .and. &
       & (my_start_end(2) /= NFC_tot) ) then

       my_start_end(2) = NFC_tot

    end if last_img

    my_len = my_start_end(2)-my_start_end(1)+1

    allocate( all_FC_indx(6, NFC_tot) )
    numfc = 1
    N2_loop: do N2 = 1, Natm_mu
        N3_loop: do N3 = 1, Natm_mu
            N4_loop: do N4 = 1, Natm_mu
                beta_loop: do beta = 1, 3
                    gama_loop: do gama = 1, 3
                        delta_loop: do delta = 1, 3

                            all_FC_indx(:, numfc) = (/N2, N3, N4, beta, gama, delta/)
                            numfc = numfc + 1

                        end do delta_loop
                    end do gama_loop
                end do beta_loop
            end do N4_loop
        end do N3_loop
    end do N2_loop

    allocate( my_FC_indx(6, my_len) )
    my_FC_indx = all_FC_indx( :, my_start_end(1) : my_start_end(2) )
    deallocate( all_FC_indx )

    img1_chk: if ( this_image() == 1 ) then
        write(*, *)
        call execute_command_line(' ')

        write(*, 10) mu, alpha
        call execute_command_line(' ')

        write(*, 20)
        call execute_command_line(' ')

    end if img1_chk

    allocate( Mat(num_row, FC4%ind_fc)[*] )
    SYNC ALL 

    write(*, 30) this_image(), my_start_end, my_FC_indx(:, 1), my_FC_indx(:, my_len)
    call execute_command_line(' ')

    Mat = 0.0_dp
    call IterOver_N2_N3_N4_mp2(mu, alpha, num_row, nstep, &
                             & my_FC_indx, sup_dim, FC4, disp, Mat)

    write(*, 40) this_image()
    call execute_command_line(' ')

    SYNC ALL 

        call co_sum(Mat, result_image=1)

    SYNC ALL 

    img1_chk2: if ( this_image() == 1 ) then

        write(*, 20)
        call execute_command_line(' ')

        write(*, 50) mu, alpha
        call execute_command_line(' ')

        write(*, '(/)')

        call find_force(sup_dim, nstep, num_row, mu, alpha, force, force_mua)

        outfile = 'FDpart_4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
        call wFD4_part(Mat, force_mua, outfile, mu, alpha)

        last_chk: if ( (mu == sys%natm) .and. (alpha == 3) ) then

            file_full = 'FD_4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
            call wFD_JoinAll(sys%natm, num_row, FC4%ind_fc, outfile, file_full)

        end if last_chk


    end if img1_chk2

    10 FORMAT(6X,'|', 33('='), ' Process started for mu = ',&
            & I4, ' and alpha = ',I3,' ', 33('='), '|')

    20 FORMAT(6X,'|',113X,'|')
    30 FORMAT(6X,'| ** Image No. = ', I4, ', FC covered ', I6, ' => ', I6, ' [(',6I4,') --> (',6I4,')] ** |')
    40 FORMAT(6X,'|', 41X, 'Image ', I4, ' completed execution', 42X, '|')

    50 FORMAT(6X,'|', 34('='), ' Process ended for mu = ',&
            & I4, ' and alpha = ',I3,' ', 34('='), '|')

end program main

