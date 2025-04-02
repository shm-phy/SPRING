
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

module FourthOrder

   use kinds,                   only : dp
   use constants,               only : Sctr4THz2, &
                                     & print_step
   use unit_cell,               only : cell
   use timer_class,             only : timer
   use hdf5_wrap,               only : w_Restart, r_Restart, &
                                     & w_ScatterMatEl4th, r_ScatterMatEl4th

   use FC4_mod,                 only : FC4type
   use phonon_m,                only : Phon
   use TetrahedronQ,            only : q_tetrahedron

   implicit none
   private

   public                   :: ScatterProbq1q2_per

contains

    subroutine ScatterProbq1q2_per( sys, FC4, ph_q0, ph_q1, q_tetra, restrt_timer, RestartDir, SctrEl4thWDir, &
                                  & SctrEl4thRDir, timelimit, i0, Nq1points, Nq1q2Perm, Nbasis, Ndof, my_Qsize, &
                                  & my_offset, QPerProcsMax, q1q2q3_permuted, WriteRestart, ReadRestart, &
                                  & StartFromBegin, WriteSctrEl4th, CalSctrEl4th, TimeUp, Sctrq1q2 )

        implicit none

        type(cell), intent(in)                                                  :: sys
        type(FC4type), intent(in)                                               :: FC4
        type(Phon), intent(in)                                                  :: ph_q0
        type(Phon), intent(in)                                                  :: ph_q1
        type(q_tetrahedron), intent(in)                                         :: q_tetra
        type(timer)                                                             :: restrt_timer

        character(len=*), intent(in)                                            :: RestartDir, SctrEl4thWDir, &
                                                                                 & SctrEl4thRDir
        real(dp), intent(in)                                                    :: timelimit

        integer, intent(in)                                                     :: i0, Nq1points, Nq1q2Perm, &
                                                                                 & Nbasis, Ndof, &
                                                                                 & my_Qsize, my_offset, &
                                                                                 & QPerProcsMax

        integer, dimension(3, Nq1q2Perm), intent(in)                            :: q1q2q3_permuted

        logical, intent(in)                                                     :: WriteRestart, ReadRestart, &
                                                                                 & StartFromBegin, WriteSctrEl4th, & 
                                                                                 & CalSctrEl4th

        logical, intent(out)                                                    :: TimeUp
        real(dp), dimension(Ndof,Ndof,Ndof,Ndof, QPerProcsMax), intent(inout)   :: Sctrq1q2

        ! ===================================== Local variables ======================================= !

        character(len=24)                                           :: frmt, img_num_char, q0_char, &
                                                                     & q1_char, q2_char, q3_char
        character(len=256)                                          :: Restrt_file, SctrFile4th, msg

        complex(dp), dimension(Ndof, Ndof)                          :: Evecq0, Evecq1, Evecq2, &
                                                                     & Evecq3
        complex(dp), dimension(:, :, :, :), allocatable             :: Phi_q0q1q2q3

        real(dp), dimension(3)                                      :: q1_float, q2_float, q3_float

        real(dp), dimension(Ndof)                                   :: OnebyOmega_q0, OnebyOmega_q1, &
                                                                     & OnebyOmega_q2, OnebyOmega_q3
        real(dp), dimension(:, :, :, :), allocatable                :: W_s0s1s2s3
        real(dp)                                                    :: multiply_factor

        integer                                                     :: ii, i1, ii_restrt, i0_restrt, &
                                                                     & complete, start, istat

        integer                                                     :: q1_grid_indx, q2_grid_indx, &
                                                                     & q3_grid_indx, img_num

        integer, dimension(3)                                       :: q0_int, q1_int, &
                                                                     & q2_int, q3_int

        logical, dimension(4)                                       :: q0chk

        ! ===================================== Local variables ======================================= !

        ! ***** Allocation ***** !
        allocate( Phi_q0q1q2q3(Ndof, Ndof, Ndof, Ndof), STAT=istat, ERRMSG=msg )
        if ( istat /= 0 ) then
            write(*, 30) this_image(), msg
            ERROR STOP
        end if

        allocate( W_s0s1s2s3(Ndof, Ndof, Ndof, Ndof), STAT=istat, ERRMSG=msg )
        if ( istat /= 0 ) then
            write(*, 30) this_image(), msg
            ERROR STOP
        end if
        30 FORMAT( " Memory allocation ERROR (ScatterProbq1q2_per) in image : ", I5, ". Error message: ", A128 )
        ! ***** Allocation ***** !

        TimeUp = .false. 

        multiply_factor = ( Sctr4THz2 / dble(Nq1points ** 2) )

        ! ------------------------------------------ Restart ------------------------------------------ !
        frmt = '(I6)'
        img_num = this_image()
        write(img_num_char, frmt) img_num
        Restrt_file = trim(adjustl(adjustr(RestartDir)))//'/Img'//trim(adjustl(adjustr(img_num_char)))//'_Restrt.h5'

        write(q0_char, frmt) i0

        RestartRead: if ( (.not. StartFromBegin) .and. ReadRestart ) then

            call r_Restart( Restrt_file, Ndof, QPerProcsMax, ii_restrt, i0_restrt, complete, Sctrq1q2 )

            start = ii_restrt + 1

            if ( i0_restrt /= i0 ) then
                write(*, 42) i0_restrt, i0
                42 FORMAT( "WARNING ( RESTART ) : Something wrong with Restart, i0_restart /= i0", I5, I5 )
            end if

            if ( this_image() == 1 ) then
                write(*, 328) (start-1), my_Qsize
                328 FORMAT(10X, 'RESTART: Restarting Calculating of Fourth-Order Sacttering Matrix elements, progress = ', &
                          & I6, '/ ', I6, ' ...')
            end if

        else RestartRead

            start = 1
            complete = 0

        end if RestartRead
        ! ------------------------------------------ Restart ------------------------------------------ !

        ! -------------------------------------------- q0 --------------------------------------------- !
        q0_int = q_tetra%irr_q_pnt_int(:, i0)

        Evecq0(:, :) = ph_q0%Evec(:, :, i0)
        OnebyOmega_q0(:) = ph_q0%OnebyOmega(:, i0)
        q0chk(1) = all( q0_int == 0 )
        ! -------------------------------------------- q0 --------------------------------------------- !

        q1_loop: do ii = start, my_Qsize

            if ( this_image() == 1 ) then

                if ( (ii == 1) .or. (mod(ii, print_step) == 0) ) then
                    if ( CalSctrEl4th ) then
                        write(*, 140) ii, my_Qsize
                    else
                        write(*, 142) ii, my_Qsize
                    end if
                    140 FORMAT(10X, 'Calculating Fourth-Order Sacttering Matrix elements, progress = ', I6, '/ ', I6, ' ...')
                    142 FORMAT(10X, 'Reading Fourth-Order Sacttering Matrix elements, progress = ', I6, '/ ', I6, ' ...')
                end if

            end if

            i1 = my_offset + ii     ! ** !

            ! ** Index of the q1, q2, q3 point. Index of q0 = i0
            q1_grid_indx = q1q2q3_permuted(1, i1)
            q2_grid_indx = q1q2q3_permuted(2, i1)
            q3_grid_indx = q1q2q3_permuted(3, i1)

            ReadOrCalculate: if ( CalSctrEl4th ) then

                ! ------------------------------------------ q1 ------------------------------------------ !
                q1_int = q_tetra%q_pnt_int(1:3, q1_grid_indx)
                q1_float = q_tetra%q_pnt(1:3, q1_grid_indx)

                Evecq1(:, :) = ph_q1%Evec(:, :, q1_grid_indx)
                OnebyOmega_q1(:) = ph_q1%OnebyOmega(:, q1_grid_indx)
                q0chk(2) = all( q1_int == 0 )
                ! ------------------------------------------ q1 ------------------------------------------ !

                ! ------------------------------------------ q2 ------------------------------------------ !
                q2_int = q_tetra%q_pnt_int(1:3, q2_grid_indx)
                q2_float = q_tetra%q_pnt(1:3, q2_grid_indx)

                Evecq2(:, :) = ph_q1%Evec(:, :, q2_grid_indx)
                OnebyOmega_q2(:) = ph_q1%OnebyOmega(:, q2_grid_indx)
                q0chk(3) = all( q2_int == 0 )
                ! ------------------------------------------ q2 ------------------------------------------ !

                ! ------------------------------------------ q3 ------------------------------------------ !
                q3_int = q_tetra%q_pnt_int(1:3, q3_grid_indx)
                q3_float = q_tetra%q_pnt(1:3, q3_grid_indx)

                Evecq3(:, :) = ph_q1%Evec(:, :, q3_grid_indx)
                OnebyOmega_q3(:) = ph_q1%OnebyOmega(:, q3_grid_indx)
                q0chk(4) = all( q3_int == 0 )
                ! ------------------------------------------ q3 ------------------------------------------ !

                ! ---------------------------- Scattering Matrix Calculation ----------------------------- !

                call ScatterMatEl4( q1_float, q2_float, q3_float, q0chk, Nbasis, Ndof, FC4, &
                                  & sys, Evecq0, Evecq1, Evecq2, Evecq3, Phi_q0q1q2q3 )

                call ScatterProbability4( Ndof, q0chk, multiply_factor, OnebyOmega_q0, OnebyOmega_q1, &
                                        & OnebyOmega_q2, OnebyOmega_q3, Phi_q0q1q2q3, W_s0s1s2s3 )

                ! ---------------------------- Scattering Matrix Calculation ----------------------------- !

                Sctrq1q2(:, :, :, :, ii) = W_s0s1s2s3(:, :, :, :)

                WriteScatterMat: if ( WriteSctrEl4th ) then

                    write(q1_char, frmt) q1_grid_indx
                    write(q2_char, frmt) q2_grid_indx
                    write(q3_char, frmt) q3_grid_indx

                    SctrFile4th = trim(adjustl(adjustr(SctrEl4thWDir)))//'/Sctr4thEl_'//&
                                 &trim(adjustl(adjustr(q0_char)))//'_'//&
                                 &trim(adjustl(adjustr(q1_char)))//'_'//&
                                 &trim(adjustl(adjustr(q2_char)))//'_'//&
                                 &trim(adjustl(adjustr(q3_char)))//'.h5'

                    call w_ScatterMatEl4th( SctrFile4th, Ndof, W_s0s1s2s3 )

                end if WriteScatterMat

                RestartWrite: if ( WriteRestart .and. &
                                 & ( (restrt_timer%elapsed_time() / 60.0_dp) > timelimit ) ) then

                    call w_Restart( Restrt_file, Ndof, QPerProcsMax, ii, i0, 0, Sctrq1q2 )

                    TimeUp = .true.

                    RETURN

                end if RestartWrite

            else ReadOrCalculate

                write(q1_char, frmt) q1_grid_indx
                write(q2_char, frmt) q2_grid_indx
                write(q3_char, frmt) q3_grid_indx

                SctrFile4th = trim(adjustl(adjustr(SctrEl4thRDir)))//'/Sctr4thEl_'//&
                             &trim(adjustl(adjustr(q0_char)))//'_'//&
                             &trim(adjustl(adjustr(q1_char)))//'_'//&
                             &trim(adjustl(adjustr(q2_char)))//'_'//&
                             &trim(adjustl(adjustr(q3_char)))//'.h5'
                
                call r_ScatterMatEl4th( SctrFile4th, Ndof, W_s0s1s2s3 )
                Sctrq1q2(:, :, :, :, ii) = W_s0s1s2s3(:, :, :, :)

            end if ReadOrCalculate

        end do q1_loop

        FinalWrite: if ( (complete /= 1) .and. WriteRestart ) then

            call w_Restart( Restrt_file, Ndof, QPerProcsMax, my_Qsize, i0, 1, Sctrq1q2 )

        end if FinalWrite

        deallocate( W_s0s1s2s3 )
        deallocate( Phi_q0q1q2q3 )

    end subroutine ScatterProbq1q2_per


    subroutine ScatterMatEl4( q1, q2, q3, q0chk, Nbasis, Ndof, FC4, sys, &
                            & Evecq, Evecq1, Evecq2, Evecq3, Phi_ss1s2s3 )

        implicit none

        real(dp), dimension(3), intent(in)                          :: q1, q2, q3!, q
        logical, dimension(4), intent(in)                           :: q0chk
        integer, intent(in)                                         :: Nbasis, Ndof
        type(FC4type), intent(in)                                   :: FC4
        type(cell), intent(in)                                      :: sys
        complex(dp), dimension(Ndof, Ndof), intent(in)              :: Evecq, Evecq1, Evecq2, &
                                                                     & Evecq3

        complex(dp), dimension(Ndof,Ndof,Ndof,Ndof), intent(out)    :: Phi_ss1s2s3

        !================================= Local variables ===================================!

        complex(dp)                                 :: expnq1R1, expn_mass, phi4_mass_expn, &
                                                     & wqs_mu_alpha, wq1s1_nu_beta, wq2s2_eta_gama, wq3s3_ro_delta, &
                                                     & expnq2R2, expnq3R3

        real(dp), dimension(3)                      :: N1R, N2R, N3R
        real(dp)                                    :: q1R1, mass_fac, q2R2, q3R3

        integer                                     :: Natm4, &
                                                     & rowno_mu_alpha, rowno_nu_beta, rowno_eta_gama, rowno_ro_delta
        integer                                     :: mu, N1, nu, N2, eta, N3, ro, &
                                                     & alpha, beta, gama, delta, &
                                                     & s, s1, s2, s3, ii
        integer, dimension(4)                       :: strt

        !================================= Local variables ===================================!

        Phi_ss1s2s3 = dcmplx(0.0_dp, 0.0_dp)

        strt(:) = 1
        do ii = 1, 4
            if ( q0chk(ii) ) strt(ii) = 4
        end do

        mu_loop: do mu = 1, Nbasis

            Natm4 = FC4%atmNum(mu)

            N1_loop: do N1 = 1, Natm4

                nu = FC4%basisTyp(N1, mu)
                N1R = FC4%atmPosR(:, N1, mu)

                q1R1 = dot_product(q1, N1R)
                expnq1R1 = dcmplx( dcos(q1R1), dsin(q1R1) )

                !expnq1R1 = cdexp( iu * q1R1 )

                N2_loop: do N2 = 1, Natm4

                    eta = FC4%basisTyp(N2, mu)
                    N2R = FC4%atmPosR(:, N2, mu)

                    q2R2 = dot_product(q2, N2R)
                    expnq2R2 = dcmplx( dcos(q2R2), dsin(q2R2) )

                    !expnq2R2 = cdexp( iu * q2R2 )

                    N3_loop: do N3 = 1, Natm4

                        ro = FC4%basisTyp(N3, mu)
                        N3R = FC4%atmPosR(:, N3, mu)

                        q3R3 = dot_product(q3, N3R)
                        expnq3R3 = dcmplx( dcos(q3R3), dsin(q3R3) )

                        !expnq3R3 = cdexp( iu * q3R3 )

                        !mass_fac = 1.0_dp / dsqrt( sys%mass(mu) * sys%mass(nu) * sys%mass(eta) * sys%mass(ro) )
                        mass_fac = ( sys%OnebySqrtMass(mu) * sys%OnebySqrtMass(nu) * &
                                   & sys%OnebySqrtMass(eta) * sys%OnebySqrtMass(ro) )

                        expn_mass = mass_fac * expnq1R1 * expnq2R2 * expnq3R3

                        alpha_loop: do alpha = 1, 3
                            rowno_mu_alpha = 3*(mu-1) + (alpha-1) + 1

                            s_loop: do s = strt(1), Ndof ! ToDo: This loop can run only up to 3 as optic modes do not !
                                                         ! contribute to thermal conduction
                                wqs_mu_alpha = Evecq(rowno_mu_alpha, s) !ToDo: make this memory access contiguous

                                beta_loop: do beta = 1, 3
                                    rowno_nu_beta = 3*(nu-1) + (beta-1) + 1

                                    s1_loop: do s1 = strt(2), Ndof
                                        wq1s1_nu_beta = Evecq1(rowno_nu_beta, s1) !ToDo: make this memory access contiguous

                                        gama_loop: do gama = 1, 3
                                            rowno_eta_gama = 3*(eta-1) + (gama-1) + 1

                                            s2_loop: do s2 = strt(3), Ndof
                                                wq2s2_eta_gama = Evecq2(rowno_eta_gama, s2) !ToDo: make this memory access &
                                                                                            ! contiguous
                                                delta_loop: do delta = 1, 3
                                                    rowno_ro_delta = 3*(ro - 1) + (delta-1) + 1

                                                    phi4_mass_expn = FC4%fC(delta, gama, beta, alpha, N3, N2, N1, mu) * &
                                                                   & expn_mass

                                                    s3_loop: do s3 = strt(4), Ndof
                                                        wq3s3_ro_delta = Evecq3(rowno_ro_delta, s3) !ToDo: make this memory access &
                                                                                                    ! contiguous
                                                        Phi_ss1s2s3(s3, s2, s1, s) = Phi_ss1s2s3(s3, s2, s1, s) + &
                                                      & phi4_mass_expn * wqs_mu_alpha * wq1s1_nu_beta * & 
                                                      & wq2s2_eta_gama * wq3s3_ro_delta

                                                    end do s3_loop

                                                end do delta_loop

                                            end do s2_loop

                                        end do gama_loop

                                    end do s1_loop

                                end do beta_loop

                            end do s_loop

                        end do alpha_loop

                    end do N3_loop

                end do N2_loop

            end do N1_loop

        end do mu_loop

    end subroutine ScatterMatEl4
    

    subroutine ScatterProbability4( Ndof, q0chk, multiply_factor, OnebyOmega_q, OnebyOmega_q1, &
                                  & OnebyOmega_q2, OnebyOmega_q3, Phi_ss1s2s3, W_ss1s2s3 )

        implicit none

        integer, intent(in)                                         :: Ndof
        logical, dimension(4), intent(in)                           :: q0chk
        real(dp), intent(in)                                        :: multiply_factor
        real(dp), dimension(Ndof), intent(in)                       :: OnebyOmega_q, OnebyOmega_q1, &
                                                                     & OnebyOmega_q2, OnebyOmega_q3
        complex(dp), dimension(Ndof,Ndof,Ndof,Ndof), intent(in)     :: Phi_ss1s2s3

        real(dp), dimension(Ndof,Ndof,Ndof,Ndof), intent(out)       :: W_ss1s2s3

        !! ===================================== Local variables =================================== !!

        real(dp), dimension(2)                      :: ReIm
        real(dp)                                    :: abs2, elmnt

        integer                                     :: ii, s, s1, s2, s3
        integer, dimension(4)                       :: strt

        !! ===================================== Local variables =================================== !!

        strt(:) = 1
        do ii = 1, 4
            if ( q0chk(ii) ) strt(ii) = 4
        end do

        W_ss1s2s3 = 0.0_dp

        s_loop2: do s = strt(1), Ndof
            s1_loop2: do s1 = strt(2), Ndof
                s2_loop2: do s2 = strt(3), Ndof
                    s3_loop2: do s3 = strt(4), Ndof

                        ReIm = transfer( Phi_ss1s2s3(s3, s2, s1, s), ReIm )
                        abs2 = dot_product( ReIm, ReIm )

                        elmnt = abs2 * ( OnebyOmega_q(s) * OnebyOmega_q1(s1) * OnebyOmega_q2(s2) * OnebyOmega_q3(s3) )

                        W_ss1s2s3(s3, s2, s1, s) = elmnt * multiply_factor

                    end do s3_loop2
                end do s2_loop2
            end do s1_loop2
        end do s_loop2

    end subroutine ScatterProbability4

end module FourthOrder

