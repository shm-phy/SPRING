
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
    use Helper,                     only : ReadPointDefVar
    use hdf5_wrap,                  only : r_force_disp, wFD4_part, wFD_JoinAll
    use DispMatFC4,                 only : IterOver_N2_N3_N4_mp2, find_force

    implicit none

    type(cell)                                              :: sys
    character(len=256)                                      :: filename
    real(dp)                                                :: T

    type(FC4_dat)                                           :: FC4
    character(len=128)                                      :: fc4file, msg

    character(len=24)                                       :: Temp_char
    character(len=8)                                        :: frmt
    character(len=128)                                      :: fdfile
    real(dp), dimension(:,:,:,:,:,:), allocatable           :: disp
    real(dp), dimension(:,:,:,:,:,:), allocatable           :: force

    integer                                                 :: mu, alpha
    real(dp), allocatable, dimension(:,:), codimension[:]   :: Mat
    real(dp), allocatable, dimension(:)                     :: force_mua, &
                                                             & Mass_def !** Variable for point-defect calculation **!

    integer, allocatable, dimension(:), codimension[:]      :: flags !** Variable for point-defect calculation **!

    integer                                                 :: Natm_mu, NFC_tot, &
                                                             & nstep, num_row, istat, &
                                                             & Nbasis, NumDef_sup !** Variable for point-defect calculation **!
    integer                                                 :: numfc, N2, N3, N4, &
                                                             & beta, gama, delta
    integer, dimension(3)                                   :: sup_dim, &
                                                             & qMesh_pd !** Variable for point-defect calculation **!
    integer, allocatable, dimension(:,:)                    :: all_FC_indx, my_FC_indx, &
                                                             & cell_atm_def !** Variable for point-defect calculation **!

    logical                                                 :: PointDef !** Variable for point-defect calculation **!

    ! -------------------------------------- Static Schedule --------------------------------------- !
    integer                                                 :: NFC_len, my_len
    integer, dimension(2)                                   :: my_start_end
    ! -------------------------------------- Static Schedule --------------------------------------- !

    ! -------------------------------------- Dynamic Schedule -------------------------------------- !
    integer, allocatable, dimension(:,:), codimension[:]    :: ChnkArr
    integer                                                 :: MyChnkStrt[*]

    logical                                                 :: DynSchd
    integer                                                 :: Nchunk, StrtImg

    integer                                                 :: ch, My_strt, My_stop, &
                                                             & img, Chnklen, ChnkStride, ChnkStrt, &
                                                             & ChnkStop, Indxloop, PickFlag
    ! -------------------------------------- Dynamic Schedule -------------------------------------- !

    character(len=64)                                       :: outfile, file_full

    if ( this_image() == 1) call ShowWelcomeBanner()

    ! =================== Initialize (All images does file reading in parallel) ==================== !

    call get_command_line(filename, T, mu, alpha, DynSchd, Nchunk, StrtImg)
    call sys%init_data(filename)
    Nbasis = sys%natm !** Point defect related **!

    fc4file = 'FC_4th_common_'//trim(sys%prefix)//'_F.h5'
    call FC4%set_FC4dat(fc4file)

    ! =================== Read Point-Defect realted input from namelist file ==================== !
    
    allocate( Mass_def(Nbasis) )
    call ReadPointDefVar( filename, Nbasis, sys, qMesh_pd, Mass_def, PointDef, NumDef_sup, cell_atm_def )
    
    ! =================== Read Point-Defect realted input from namelist file ==================== !

    frmt = '(F6.1)'
    write(Temp_char, frmt) T
    fdfile = 'disp_forc_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
    call r_force_disp(fdfile, disp, force)

    Natm_mu = FC4%atmNum(mu)
    NFC_tot = (Natm_mu**3) * (3**3)

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

    allocate( Mat(num_row, FC4%ind_fc)[*], STAT=istat, ERRMSG=msg )
    if(istat /= 0) write(*, "('Memory allocation failed in image no.: ', I6, '. ERROR MESSAGE: ', A128 )") this_image(), msg
    SYNC ALL 
    Mat = 0.0_dp

    allocate( flags(num_row)[*], STAT=istat, ERRMSG=msg )
    if(istat /= 0) write(*, "('Memory allocation failed in image no.: ', I6, '. ERROR MESSAGE: ', A128 )") this_image(), msg
    SYNC ALL 
    flags = 0
    ! =================== Initialize (All images does file reading in parallel) ==================== !

    img1_chk: if ( this_image() == 1 ) then
        write(*, *)
        call execute_command_line(' ')

        write(*, 10) mu, alpha
        call execute_command_line(' ')

        write(*, 20)
        call execute_command_line(' ')

    end if img1_chk

    DynSchdChk: if ( DynSchd ) then

        Chnklen = NFC_tot / Nchunk
        ErrorStop1: if ( Chnklen == 0 ) then

            write(*, 150) NFC_tot, NChunk
            ERROR STOP

        end if ErrorStop1

        allocate( ChnkArr(4, Nchunk)[*] )
        SYNC ALL

        OnlyImg1: if ( this_image() == 1 ) then

            chnkLoop: do ch = 1, Nchunk

                My_strt = ((ch - 1) * Chnklen) + 1
                My_stop = My_strt + Chnklen - 1

                if ( (ch == Nchunk) .and. (My_stop /= NFC_tot) ) My_stop = NFC_tot

                ChnkArr(:, ch) = (/0, My_strt, My_stop, 0/)

            end do chnkLoop

            NotSingleImg1: if ( num_images() > 1 ) then

                ! ------====== Allocate ChynkStride (=Nchunk / num_images()) Chnunks to each image ======------ !

                ChnkStride = Nchunk / ( num_images() - 1 )

                ErrorStop: if ( ChnkStride == 0 ) then 

                    write(*, 120) Nchunk, ( num_images() - 1 )
                    ERROR STOP

                end if ErrorStop

                distrImage: do img = 1, (num_images() - 1)

                    ChnkStrt = ((img - 1) * ChnkStride) + 1
                    ChnkStop = ChnkStrt + ChnkStride - 1

                    if ( (img == (num_images()-1) ) .and. (ChnkStop /= Nchunk) ) ChnkStop = Nchunk

                    ! ** img+1 ** !
                    MyChnkStrt[img+1] = ChnkStrt
                    ! ** img+1 ** !

                    ChnkArr(4, ChnkStrt:ChnkStop) = img+1

                end do distrImage

                ! ------====== Allocate ChynkStride (=Nchunk / num_images()) Chnunks to each image ======------ !

                ! image 1 points to :
                MyChnkStrt[1] = ChnkStride / 2
                if ( (ChnkStride / 2) == 0 ) MyChnkStrt[1] = 2

            else NotSingleImg1

                MyChnkStrt = 1
                ChnkArr(4, :) = 1

            end if NotSingleImg1

        end if OnlyImg1

        SYNC ALL

        ! * Broadcast to all images * !
        call co_broadcast(ChnkArr, source_image=1, STAT=istat, ERRMSG=msg)
        if ( istat /= 0 ) write(*, "( 'ERROR in co_broadcast : ', A128 )") msg

        SYNC ALL

        Indxloop = MyChnkStrt

        NotSingleImag2: if ( num_images() > 1 ) then

            NotImage1: if ( this_image() >= StrtImg ) then

                DynSchdl: do
#ifdef GNU
                    AccessImg1: CRITICAL

                        PickFlag = ChnkArr(1, Indxloop) [1]
                        if ( PickFlag == 0 ) ChnkArr(1, Indxloop) [1] = 1

                    END CRITICAL AccessImg1
#else
                    ! *OpenCoarrays have not implemented atomic operations yet* !
                    call atomic_fetch_or( ChnkArr(1, Indxloop)[1], 1, PickFlag, STAT=istat )
                    if (istat /= 0) write(*, "( 'ERROR in atomic_fetch_or' )")
#endif
                    StartWork: if ( PickFlag == 0 ) then

                        my_start_end = ChnkArr(2:3, Indxloop)

                        my_len = my_start_end(2)-my_start_end(1)+1

                        allocate( my_FC_indx(6, my_len) )
                        my_FC_indx = all_FC_indx( :, my_start_end(1) : my_start_end(2) )

                        write(*, 30) this_image(), my_start_end, my_FC_indx(:, 1), my_FC_indx(:, my_len)
                        call execute_command_line(' ')

                        call IterOver_N2_N3_N4_mp2( mu, alpha, num_row, nstep, NumDef_sup, PointDef, &
                                                  & cell_atm_def, my_FC_indx, sup_dim, FC4, disp, &
                                                  & flags, Mat )

                        write(*, 40) this_image()
                        call execute_command_line(' ')

                        deallocate( my_FC_indx )

                    end if StartWork

                    Indxloop = mod( Indxloop, Nchunk ) + 1

                    if ( Indxloop == MyChnkStrt ) EXIT DynSchdl

                end do DynSchdl

            end if NotImage1

        else NotSingleImag2

            OnlyImg1Work: do 

                PickFlag = ChnkArr(1, Indxloop)

                if ( PickFlag == 0 ) ChnkArr(1, Indxloop) = 1

                StartWork2: if ( PickFlag == 0 ) then

                    my_start_end = ChnkArr(2:3, Indxloop)

                    my_len = my_start_end(2)-my_start_end(1)+1

                    allocate( my_FC_indx(6, my_len) )
                    my_FC_indx = all_FC_indx( :, my_start_end(1) : my_start_end(2) )

                    write(*, 30) this_image(), my_start_end, my_FC_indx(:, 1), my_FC_indx(:, my_len)
                    call execute_command_line(' ')

                    call IterOver_N2_N3_N4_mp2( mu, alpha, num_row, nstep, NumDef_sup, PointDef, &
                                              & cell_atm_def, my_FC_indx, sup_dim, FC4, disp, &
                                              & flags, Mat )

                    write(*, 40) this_image()
                    call execute_command_line(' ')

                    deallocate( my_FC_indx )

                end if StartWork2

                Indxloop = mod( Indxloop, Nchunk ) + 1

                if ( Indxloop == MyChnkStrt ) EXIT OnlyImg1Work

            end do OnlyImg1Work

        end if NotSingleImag2

    else DynSchdChk

        SYNC ALL 
        ! ================================= Static Scheduling ================================= !

        NFC_len = NFC_tot / num_images()
        my_start_end(1) = ((this_image() - 1) * NFC_len) + 1
        my_start_end(2) = my_start_end(1) + NFC_len - 1

        last_img: if ( (this_image() == num_images()) .and. &
                     & (my_start_end(2) /= NFC_tot) ) then

           my_start_end(2) = NFC_tot

        end if last_img

        my_len = my_start_end(2)-my_start_end(1)+1

        allocate( my_FC_indx(6, my_len) )
        my_FC_indx = all_FC_indx( :, my_start_end(1) : my_start_end(2) )
        deallocate( all_FC_indx )

        write(*, 30) this_image(), my_start_end, my_FC_indx(:, 1), my_FC_indx(:, my_len)
        call execute_command_line(' ')

        call IterOver_N2_N3_N4_mp2( mu, alpha, num_row, nstep, NumDef_sup, PointDef, &
                                  & cell_atm_def, my_FC_indx, sup_dim, FC4, disp, &
                                  & flags, Mat )

        deallocate( my_FC_indx )

        write(*, 40) this_image()
        call execute_command_line(' ')

        ! ================================= Static Scheduling ================================= !

    end if DynSchdChk

    SYNC ALL 

        call co_sum(Mat, result_image=1, STAT=istat, ERRMSG=msg)
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

    SYNC ALL 

        call co_sum(flags, result_image=1, STAT=istat, ERRMSG=msg)
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

    SYNC ALL 

    img1_chk2: if ( this_image() == 1 ) then

        write(*, 20)
        call execute_command_line(' ')

        write(*, 50) mu, alpha
        call execute_command_line(' ')

        write(*, '(/)')

        call find_force(sup_dim, nstep, num_row, mu, alpha, force, force_mua)

        outfile = 'FDpart_4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
        call wFD4_part(Mat, force_mua, flags, outfile, mu, alpha)

        last_chk: if ( (mu == sys%natm) .and. (alpha == 3) ) then

            write(*, *)
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

    120 FORMAT('Nchunk ( ', I6, ' ) must be greater than or equal to (Nprocs-1) ( ', I5, ' )')
    150 FORMAT('Number of FCs ( ', I7, ' ) must be greater than or equal to Nchunk ( ', I5, ' )')

end program main

