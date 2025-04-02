
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
    use hdf5_wrap,                  only : wFC3_Reconst, ShengBTE_FC3rd
    use ReconstFC3,                 only : IterOver_N2_N3_mp2

    implicit none

    type(cell)                                                          :: sys
    character(len=256)                                                  :: filename
    real(dp)                                                            :: T

    type(FC3_dat)                                                       :: FC3
    character(len=128)                                                  :: fc3file, &
                                                                         & fc3indfc_file, msg
    character(len=24)                                                   :: Temp_char
    character(len=8)                                                    :: frmt

    real(dp), allocatable, dimension(:,:,:,:,:,:), codimension[:]       :: FC3_Reconst

    integer                                                             :: Nbasis, Natm_mu, NFC_tot, &
                                                                         & NFC_len, my_len, istat
    integer, dimension(2)                                               :: my_start_end
    integer                                                             :: numfc, mu, N2, N3, &
                                                                         & alpha, beta, gama
    integer, allocatable, dimension(:,:)                                :: all_FC_indx, my_FC_indx

    logical                                                             :: write_sBTE
    character(len=64)                                                   :: outfile



    if ( this_image() == 1) call ShowWelcomeBanner()

    call get_command_line(filename, T, write_sBTE)
    call sys%init_data(filename)

    frmt = '(F6.1)'
    write(Temp_char, frmt) T

    fc3file = 'FC_3rd_common_'//trim(sys%prefix)//'_F.h5'
    fc3indfc_file = 'FD_3rd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

    call FC3%set_FC3dat(fc3file, fc3indfc_file)

    Nbasis = sys%natm
    Natm_mu =  maxval(FC3%atmNum)
    NFC_tot = Nbasis * (Natm_mu**2) * (3**3)
    NFC_len = NFC_tot / num_images()

    my_start_end(1) = ((this_image() - 1) * NFC_len) + 1
    my_start_end(2) = my_start_end(1) + NFC_len - 1

    last_img: if ( (this_image() == num_images()) .and. &
                 & (my_start_end(2) /= NFC_tot) ) then

       my_start_end(2) = NFC_tot

    end if last_img

    my_len = my_start_end(2)-my_start_end(1)+1

    allocate( all_FC_indx(6, NFC_tot) )
    numfc = 1
    mu_loop: do mu = 1, Nbasis
        N2_loop: do N2 = 1, Natm_mu
            N3_loop: do N3 = 1, Natm_mu
                alpha_loop: do alpha = 1, 3
                    beta_loop: do beta = 1, 3
                        gama_loop: do gama = 1, 3

                            all_FC_indx(:, numfc) = (/mu, N2, N3, alpha, beta, gama/)
                            numfc = numfc + 1

                        end do gama_loop
                    end do beta_loop
                end do alpha_loop
            end do N3_loop
        end do N2_loop
    end do mu_loop

    allocate( my_FC_indx(6, my_len) )
    my_FC_indx = all_FC_indx( :, my_start_end(1) : my_start_end(2) )
    deallocate( all_FC_indx )

    img1_chk: if ( this_image() == 1 ) then
        write(*, *)
        call execute_command_line(' ')

        write(*, 10) 
        call execute_command_line(' ')

        write(*, 20)
        call execute_command_line(' ')

    end if img1_chk

    allocate( FC3_Reconst(3,3,3,Natm_mu,Natm_mu,Nbasis)[*] )
    SYNC ALL 

    write(*, 30) this_image(), my_start_end, my_FC_indx(:, 1), my_FC_indx(:, my_len)
    call execute_command_line(' ')

    FC3_Reconst = 0.0_dp
    call IterOver_N2_N3_mp2(my_FC_indx, FC3, FC3_Reconst)

    write(*, 40) this_image()
    call execute_command_line(' ')

    SYNC ALL 

        call co_sum(FC3_Reconst, result_image=1, STAT=istat, ERRMSG=msg)
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

    SYNC ALL 

    img1_chk2: if ( this_image() == 1 ) then

        write(*, 20)
        call execute_command_line(' ')

        write(*, 50) 
        call execute_command_line(' ')

        write(*, '(/)')

        outfile = 'FC3rd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
        call wFC3_Reconst(outfile, FC3_Reconst, FC3%atmNum, FC3%cb_Indx)

        wirte_ShengBTE: if ( write_sBTE ) then

            call ShengBTE_FC3rd(sys, T, FC3_Reconst, FC3%atmNum, FC3%cb_Indx)

        end if wirte_ShengBTE

    end if img1_chk2

    10 FORMAT(6x,'|', 56('='),&
                &' Process started ', 58('='), '|')

    20 FORMAT(6X,'|',131X,'|')
    30 FORMAT(6X,'| **', 8X,' Image No. = ', I4, ', FC covered ', I7, ' => ', I7, ' [(',6I4,') --> (',6I4,')]', 8X,' ** |')
    40 FORMAT(6X,'|', 50X, 'Image ', I4, ' completed execution', 51X, '|')

    50 FORMAT(6x,'|', 57('='), &
                &' Process ended ', 59('='), '|')

end program main

