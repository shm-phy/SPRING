
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


module FreeEnergy3rdTerm2

   use kinds,          only : dp
   use constants,      only : eVConv3_t2, print_step2, &
                            & print_step
   use unit_cell,      only : cell

   use Irr_q_point,    only : q_points_data
   use FC3_mod,        only : FC3type
   use phonon_m,       only : Phon

   implicit none
   private

   public                   :: ScatterMatAllqq2, FindFreeEngergy3rdT2

contains

    subroutine ScatterMatAllqq2( Qpoints, mesh, Nbasis, my_Qsize, my_offset, &
                               & sys, FC3, phonon_dat, Phiqq2 )

        implicit none


        type(q_points_data), intent(in)                             :: Qpoints
        integer, dimension(3), intent(in)                           :: mesh
        integer, intent(in)                                         :: Nbasis, my_Qsize, my_offset
        type(cell), intent(in)                                      :: sys
        type(FC3type), intent(in)                                   :: FC3
        type(Phon), intent(in)                                      :: phonon_dat

        complex(dp), dimension(:,:,:,:), intent(inout)              :: Phiqq2

        !================================= Local variable ===================================!

        complex(dp), dimension(:, :), allocatable                   :: Evecq, Evecmq, Evecq2
        complex(dp), dimension(:, :, :), allocatable                :: Phi_qmqq2

        real(dp), dimension(3)                                      :: q_float, mq_float, q2_float

        integer                                                     :: num_grid_pnt
        integer                                                     :: ii, i0, mi0, i2, PRINT_IMG
        integer, dimension(3)                                       :: q, mq, q2, mesh_mul

        logical, dimension(3)                                       :: q0chk

        !================================= Local variable ===================================!

        PRINT_IMG = ( 1 + num_images() ) / 2

        allocate( Evecq(3*Nbasis, 3*Nbasis) )
        allocate( Evecmq(3*Nbasis, 3*Nbasis) )
        allocate( Evecq2(3*Nbasis, 3*Nbasis) )

        allocate( Phi_qmqq2(3*Nbasis, 3*Nbasis, 3*Nbasis) )

        num_grid_pnt = product( mesh )
        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)

        q2_float(:) = 0.0_dp
        q2(:) = 0
        i2 = Qpoints%indx_map( dot_product( modulo( q2, mesh ), mesh_mul ) + 1 )

        Evecq2 = phonon_dat%Evec(:, :, i2)
        q0chk(3) = all( q2 == 0 )

        q_loop: do ii = 1, my_Qsize

            if ( (this_image() == PRINT_IMG) .and. (mod(ii, print_step2) == 0) ) then
                write(*, 125) ii, my_Qsize
                call execute_command_line(' ')
                125 FORMAT(10X, 'Calculating Scattering Matrix-elements(3rd), q-index = ', I5, '/ ', I5, ' ...')
            end if

            i0 = my_offset + ii

            ! ---------------------------------------------------------------------------- !
            q = Qpoints%q_pnt_int(1:3, i0)
            q_float = Qpoints%q_pnt(:, i0)
            Evecq = phonon_dat%Evec(:, :, i0)

            mq = -1 * q
            mq_float = -1.0_dp * q_float
            mi0 = Qpoints%indx_map( dot_product( modulo( mq, mesh ), mesh_mul ) + 1 )
            Evecmq = phonon_dat%Evec(:, :, mi0)

            q0chk(1) = all( q == 0 )
            q0chk(2) = all( mq == 0 )

            call ScatterMatEl( mq_float, q2_float, q0chk, Nbasis, FC3, sys, &
                            & Evecq, Evecmq, Evecq2, Phi_qmqq2)

            ! ---------------------------------------------------------------------------- !

            Phiqq2(:,:,:, i0) = Phi_qmqq2

        end do q_loop

        deallocate( Evecq, Evecmq, Evecq2 )
        deallocate( Phi_qmqq2 )

    end subroutine ScatterMatAllqq2


    subroutine FindFreeEngergy3rdT2( Qpoints, mesh, Nbasis, my_Qsize, my_qq1, &
                                   & phonon_dat, Phiqq2, FreeEng3rd_t2 )

        implicit none

        type(q_points_data), intent(in)                             :: Qpoints
        integer, dimension(3), intent(in)                           :: mesh
        integer, intent(in)                                         :: Nbasis, my_Qsize
        integer, dimension(:, :), intent(in)                        :: my_qq1

        type(Phon), intent(in)                                      :: phonon_dat

        complex(dp), dimension(:,:,:,:), intent(in)                 :: Phiqq2

        complex(dp), intent(out)                                    :: FreeEng3rd_t2

        !================================= Local variable ===================================!

        complex(dp)                                                 :: FreeEng
        complex(dp), dimension(:, :, :), allocatable                :: Phi_qmqq2, Phi_q1mq1q2

        real(dp)                                                    :: omega_cutoff
        real(dp), dimension(:), allocatable                         :: omega_q, omega_q1, omega_q2
        real(dp), dimension(:), allocatable                         :: nBE_q, nBE_q1, nBE_q2

        integer                                                     :: num_grid_pnt, PRINT_IMG
        integer                                                     :: ii, i0, i1, i2, &
                                                                     & mi0, mi1
        integer, dimension(3)                                       :: q, q1, q2, mesh_mul, &
                                                                     & mq, mq1
        logical, dimension(3)                                       :: q0chk

        !================================= Local variable ===================================!

        PRINT_IMG = ( 1 + num_images() ) / 2

        FreeEng3rd_t2 = dcmplx(0.0_dp, 0.0_dp)

        allocate( omega_q(3*Nbasis) )
        allocate( omega_q1(3*Nbasis) )
        allocate( omega_q2(3*Nbasis) )

        allocate( nBE_q(3*Nbasis) )
        allocate( nBE_q1(3*Nbasis) )
        allocate( nBE_q2(3*Nbasis) )

        allocate( Phi_qmqq2(3*Nbasis, 3*Nbasis, 3*Nbasis) )
        allocate( Phi_q1mq1q2(3*Nbasis, 3*Nbasis, 3*Nbasis) )

        omega_cutoff = phonon_dat%omega_min

        num_grid_pnt = product( mesh )
        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)

        q2(:) = 0
        i2 = Qpoints%indx_map( dot_product( modulo( q2, mesh ), mesh_mul ) + 1 )

        omega_q2 = phonon_dat%omega(:, i2)
        nBE_q2 = phonon_dat%nBE(:, i2)

        q0chk(3) = all( q2 == 0 )

        qq1_loop: do ii = 1, my_Qsize

            if ( (this_image() == PRINT_IMG) .and. (mod(ii, print_step) == 0) ) then
                write(*, 130) ii, my_Qsize
                call execute_command_line(' ')
                130 FORMAT(10X, 'Calculating Free-Energy-3rd (Term-2), progress = ', I5, '/ ', I5, ' ...')
            end if

            i0 = my_qq1(1, ii)
            i1 = my_qq1(2, ii)

            ! ---------------------------------------------------------------------------- !
            q = Qpoints%q_pnt_int(1:3, i0)
            omega_q = phonon_dat%omega(:, i0)
            nBE_q = phonon_dat%nBE(:, i0)

            q0chk(1) = all( q == 0 )

            mq = -1 * q
            mi0 = Qpoints%indx_map( dot_product( modulo( mq, mesh ), mesh_mul ) + 1 )

            Phi_qmqq2 = Phiqq2(:,:,:, i0)
            ! ---------------------------------------------------------------------------- !


            ! ---------------------------------------------------------------------------- !
            q1 = Qpoints%q_pnt_int(1:3, i1)
            omega_q1 = phonon_dat%omega(:, i1)
            nBE_q1 = phonon_dat%nBE(:, i1)

            q0chk(2) = all( q1 == 0 )

            mq1 = -1 * q1
            mi1 = Qpoints%indx_map( dot_product( modulo( mq1, mesh ), mesh_mul ) + 1 )

            Phi_q1mq1q2 = Phiqq2(:,:,:, i1)
            ! ---------------------------------------------------------------------------- !
            call MultiplyNqsOmegas( Nbasis, num_grid_pnt, q0chk, omega_cutoff, omega_q, omega_q1, omega_q2, &
                                  & nBE_q, nBE_q1, Phi_qmqq2, Phi_q1mq1q2, FreeEng )

            FreeEng3rd_t2 = FreeEng3rd_t2 + (-1.0_dp * eVConv3_t2) * FreeEng

        end do qq1_loop

        deallocate( omega_q, omega_q1, omega_q2 )
        deallocate( nBE_q, nBE_q1, nBE_q2 )
        deallocate( Phi_qmqq2, Phi_q1mq1q2 )

    end subroutine FindFreeEngergy3rdT2


    subroutine MultiplyNqsOmegas( Nbasis, Nq, q0chk, omega_cutoff, omega_q, omega_q1, omega_q2, &
                                & nBE_q, nBE_q1, Phi_qmqq2, Phi_q1mq1q2, FreeEng )

        implicit none

        integer, intent(in)                                         :: Nbasis, Nq
        logical, dimension(3), intent(in)                           :: q0chk
        real(dp), intent(in)                                        :: omega_cutoff
        real(dp), dimension(:), intent(in)                          :: omega_q, omega_q1, omega_q2
        real(dp), dimension(:), intent(in)                          :: nBE_q, nBE_q1

        complex(dp), dimension(:,:,:), intent(in)                   :: Phi_qmqq2, Phi_q1mq1q2

        complex(dp), intent(out)                                    :: FreeEng

        !! ===================================== Local variables =================================== !!

        real(dp)                                    :: omega_s, omega_s1, omega_s2, &
                                                     & n_s, n_s1
        real(dp)                                    :: denom, mul_term, denom_sctr

        integer                                     :: Ndof, s, s1, s2, ii
        integer, dimension(3)                       :: strt

        logical                                     :: not_zero

        !! ===================================== Local variables =================================== !!

        FreeEng = dcmplx(0.0_dp, 0.0_dp)

        Ndof = 3 * Nbasis

        strt(:) = 1
        do ii = 1, 3
            if ( q0chk(ii) ) strt(ii) = 4
        end do

        s_loop: do s = strt(1), Ndof
            omega_s = omega_q(s)
            n_s = nBE_q(s)

            s1_loop: do s1 = strt(2), Ndof
                omega_s1 = omega_q1(s1)
                n_s1 = nBE_q1(s1)

                s2_loop: do s2 = 4, Ndof
                    omega_s2 = omega_q2(s2)

                    denom_sctr = ( omega_s * omega_s1 * omega_s2 * dble(Nq) )

                    mul_term = 0.0_dp
                    
                    denom = omega_s2

                    not_zero = ( dabs( denom ) > omega_cutoff )

                    AvoidDivg1: if ( not_zero ) then
                        mul_term = (n_s * n_s1 + n_s + 0.25_dp) / denom

                    end if AvoidDivg1

                    AvoidDivg2: if ( denom_sctr > omega_cutoff ) then

                        mul_term = mul_term / denom_sctr

                    else AvoidDivg2

                        mul_term = 0.0_dp

                        write(*, 23)
                        23 FORMAT( 'WARNING: denom_sctr is zero !' )

                    end if AvoidDivg2

                    FreeEng = FreeEng + ( Phi_qmqq2(s2, s, s) * Phi_q1mq1q2(s2, s1, s1) * mul_term )

                end do s2_loop
            end do s1_loop
        end do s_loop

    end subroutine MultiplyNqsOmegas


    subroutine ScatterMatEl(q1, q2, q0chk, Nbasis, FC3, sys, &
                          & Evecq, Evecq1, Evecq2, Phi_ss1s2)

        implicit none

        real(dp), dimension(3), intent(in)                          :: q1, q2!, q
        logical, dimension(3), intent(in)                           :: q0chk
        integer, intent(in)                                         :: Nbasis
        type(FC3type), intent(in)                                   :: FC3
        type(cell), intent(in)                                      :: sys
        complex(dp), dimension(:, :), intent(in)                    :: Evecq, Evecq1, Evecq2

        complex(dp), dimension(:,:,:), intent(out)                  :: Phi_ss1s2

        !================================= Local variables ===================================!

        complex(dp)                                 :: expnq1R1, expnq2R2, expn_mass, phi3_mass_expn, &
                                                     & wqs_mu_alpha, wq1s1_nu_beta, wq2s2_eta_gama

        real(dp), dimension(3)                      :: N1R, N2R
        real(dp)                                    :: q1R1, q2R2, mass_fac

        integer                                     :: Natm3, Ndof, &
                                                     & rowno_mu_alpha, rowno_nu_beta, rowno_eta_gama
        integer                                     :: mu, N1, nu, N2, eta, &
                                                     & alpha, beta, gama, &
                                                     & s, s1, s2, ii
        integer, dimension(3)                       :: strt

        !================================= Local variables ===================================!

        Ndof = 3*Nbasis

        Phi_ss1s2 = dcmplx(0.0_dp, 0.0_dp)

        strt(:) = 1
        do ii = 1, 3
            if ( q0chk(ii) ) strt(ii) = 4
        end do

        mu_loop: do mu = 1, Nbasis

            Natm3 = FC3%atmNum(mu)

            N1_loop: do N1 = 1, Natm3

                nu = FC3%basisTyp(N1, mu)
                N1R = FC3%atmPosR(:, N1, mu)

                q1R1 = dot_product(q1, N1R)
                !expnq1R1 = cdexp( iu * q1R1 )

                expnq1R1 = dcmplx( dcos(q1R1), dsin(q1R1) )

                N2_loop: do N2 = 1, Natm3

                    eta = FC3%basisTyp(N2, mu)
                    N2R = FC3%atmPosR(:, N2, mu)

                    q2R2 = dot_product(q2, N2R)
                    !expnq2R2 = cdexp( iu * q2R2 )

                    expnq2R2 = dcmplx( dcos(q2R2), dsin(q2R2) )

                    mass_fac = 1.0_dp / dsqrt(  sys%mass(mu) *  sys%mass(nu) * sys%mass(eta) )

                    expn_mass = expnq1R1 * expnq2R2 * mass_fac

                    alpha_loop: do alpha = 1, 3
                        rowno_mu_alpha = 3*(mu-1) + (alpha-1) + 1

                        s_loop: do s = strt(1), Ndof
                            wqs_mu_alpha = Evecq(rowno_mu_alpha, s)

                            beta_loop: do beta = 1, 3
                                rowno_nu_beta = 3*(nu-1) + (beta-1) + 1

                                s1_loop: do s1 = strt(2), Ndof
                                    wq1s1_nu_beta = Evecq1(rowno_nu_beta, s1)

                                    gama_loop: do gama = 1, 3
                                        rowno_eta_gama = 3*(eta-1) + (gama-1) + 1

                                        phi3_mass_expn = FC3%fC(gama, beta, alpha, N2, N1, mu) * expn_mass

                                        s2_loop: do s2 = 4, Ndof
                                            wq2s2_eta_gama = Evecq2(rowno_eta_gama, s2) 
                                            
                                            Phi_ss1s2(s2, s1, s) = Phi_ss1s2(s2, s1, s) + &
                                          & phi3_mass_expn * wqs_mu_alpha * wq1s1_nu_beta * wq2s2_eta_gama

                                        end do s2_loop

                                    end do gama_loop

                                end do s1_loop

                            end do beta_loop

                        end do s_loop

                    end do alpha_loop

                end do N2_loop

            end do N1_loop

        end do mu_loop

    end subroutine ScatterMatEl


end module FreeEnergy3rdTerm2

