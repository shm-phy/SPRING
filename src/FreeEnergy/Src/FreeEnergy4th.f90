
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

module FreeEnergy4th

    use kinds,          only : dp
    use unit_cell,      only : cell
    use constants,      only : eVConv4, print_step

    use Irr_q_point,    only : q_points_data
    use FC4_mod,        only : FC4type
    use phonon_m,       only : Phon

    implicit none
    private

    public  ::  FindFreeEnergy4th

contains

    subroutine FindFreeEnergy4th( Qpoints, mesh, Nbasis, my_Qsize, my_qq1, &
                                & sys, FC4, phonon_dat, FreeEng4th )

        implicit none

        type(q_points_data), intent(in)                             :: Qpoints
        integer, dimension(3), intent(in)                           :: mesh
        integer, intent(in)                                         :: Nbasis, my_Qsize
        integer, dimension(:, :), intent(in)                        :: my_qq1

        type(cell), intent(in)                                      :: sys
        type(FC4type), intent(in)                                   :: FC4
        type(Phon), intent(in)                                      :: phonon_dat

        complex(dp), intent(out)                                    :: FreeEng4th

        !================================= Local variable ===================================!

        complex(dp)                                                 :: FreeEng
        complex(dp), dimension(:,:), allocatable                    :: Evecq, Evecmq, Evecq1, Evecmq1
        complex(dp), dimension(:,:,:,:), allocatable                :: Phiqq1

        real(dp)                                                    :: omega_cutoff
        real(dp), dimension(3)                                      :: q_float, mq_float, q1_float, mq1_float
        real(dp), dimension(:), allocatable                         :: omega_q, omega_q1
        real(dp), dimension(:), allocatable                         :: nBE_q, nBE_q1

        integer                                                     :: num_grid_pnt, PRINT_IMG
        integer                                                     :: ii, i0, mi0, i1, mi1
        integer, dimension(3)                                       :: q, mq, q1, mq1, mesh_mul

        logical, dimension(4)                                       :: q0chk

        !================================= Local variable ===================================!

        PRINT_IMG = ( 1 + num_images() ) / 2

        FreeEng4th = dcmplx(0.0_dp, 0.0_dp)

        allocate( Evecq(3*Nbasis, 3*Nbasis) )
        allocate( Evecmq(3*Nbasis, 3*Nbasis) )
        allocate( Evecq1(3*Nbasis, 3*Nbasis) )
        allocate( Evecmq1(3*Nbasis, 3*Nbasis) )

        allocate( Phiqq1(3*Nbasis, 3*Nbasis, 3*Nbasis, 3*Nbasis) )

        allocate( omega_q(3*Nbasis) )
        allocate( omega_q1(3*Nbasis) )
        allocate( nBE_q(3*Nbasis) )
        allocate( nBE_q1(3*Nbasis) )

        num_grid_pnt = product( mesh )
        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)

        omega_cutoff = phonon_dat%omega_min

        qq1_loop: do ii = 1, my_Qsize
            
            if ( (this_image() == PRINT_IMG) .and. (mod(ii, print_step) == 0) ) then
                write(*, 140) ii, my_Qsize
                140 FORMAT(10X, 'Calculating Free-Energy-4th, progress = ', I5, '/ ', I5, ' ...')
            end if

            !! -------------------------------------------------------------------------- !!
            i0 = my_qq1(1, ii)
            q = Qpoints%q_pnt_int(1:3, i0)
            q_float = Qpoints%q_pnt(:, i0)

            Evecq = phonon_dat%Evec(:,:, i0)
            omega_q = phonon_dat%omega(:, i0)
            nBE_q = phonon_dat%nBE(:, i0)

            q0chk(1) = all( q == 0 )

            mq = -1 * q
            mi0 = Qpoints%indx_map( dot_product( modulo( mq, mesh ), mesh_mul ) + 1 )
            mq_float = -1.0_dp * q_float

            Evecmq = phonon_dat%Evec(:,:, mi0)

            q0chk(2) = all( mq == 0 )
            !! -------------------------------------------------------------------------- !!

            !! -------------------------------------------------------------------------- !!
            i1 = my_qq1(2, ii)
            q1 = Qpoints%q_pnt_int(1:3, i1)
            q1_float = Qpoints%q_pnt(:, i1)

            Evecq1 = phonon_dat%Evec(:,:, i1)
            omega_q1 = phonon_dat%omega(:, i1)
            nBE_q1 = phonon_dat%nBE(:, i1)

            q0chk(3) = all( q1 == 0 )

            mq1 = -1 * q1
            mi1 = Qpoints%indx_map( dot_product( modulo( mq1, mesh ), mesh_mul ) + 1 )
            mq1_float = -1.0_dp * q1_float

            Evecmq1 = phonon_dat%Evec(:,:, mi1)

            q0chk(4) = all( mq1 == 0 )
            !! -------------------------------------------------------------------------- !!

            call ScatterMatEl4( mq_float, q1_float, mq1_float, q0chk, Nbasis, FC4, sys, &
                              & Evecq, Evecmq, Evecq1, Evecmq1, Phiqq1 )

            call MultiplynBE_omega( Nbasis, num_grid_pnt, q0chk, omega_cutoff, omega_q, omega_q1, &
                                  & nBE_q, nBE_q1, Phiqq1, FreeEng )

            FreeEng4th = FreeEng4th + eVConv4 * FreeEng

        end do qq1_loop

        deallocate( Evecq, Evecmq, Evecq1, Evecmq1 )
        deallocate( omega_q, omega_q1, nBE_q, nBE_q1 )
        deallocate( Phiqq1 )

    end subroutine FindFreeEnergy4th


    subroutine MultiplynBE_omega( Nbasis, Nq, q0chk, omega_cutoff, omega_q, omega_q1, &
                                & nBE_q, nBE_q1, Phiqq1, FreeEng )

        implicit none

        integer, intent(in)                                         :: Nbasis, Nq
        logical, dimension(4), intent(in)                           :: q0chk
        real(dp), intent(in)                                        :: omega_cutoff
        real(dp), dimension(:), intent(in)                          :: omega_q, omega_q1
        real(dp), dimension(:), intent(in)                          :: nBE_q, nBE_q1

        complex(dp), dimension(:,:,:,:), intent(in)                 :: Phiqq1

        complex(dp), intent(out)                                    :: FreeEng

        !! ===================================== Local variables =================================== !!

        real(dp)                                    :: omega_s, omega_s1, &
                                                     & n_s, n_s1
        real(dp)                                    :: mul_term, denom_sctr

        integer                                     :: Ndof, s, s1, ii
        integer, dimension(4)                       :: strt

        !! ===================================== Local variables =================================== !!

        FreeEng = dcmplx(0.0_dp, 0.0_dp)

        Ndof = 3 * Nbasis

        strt(:) = 1
        do ii = 1, 4
            if ( q0chk(ii) ) strt(ii) = 4
        end do

        s_loop: do s = strt(1), Ndof
            omega_s = omega_q(s)
            n_s = nBE_q(s)

            s1_loop: do s1 = strt(3), Ndof
                omega_s1 = omega_q1(s1)
                n_s1 = nBE_q1(s1)

                denom_sctr = ( omega_s * omega_s1 * dble(Nq) )

                mul_term = (n_s + 0.5_dp) * (n_s1 + 0.5_dp)

                AvoidDivg: if ( denom_sctr > omega_cutoff ) then

                    mul_term = mul_term / denom_sctr

                else AvoidDivg

                    mul_term = 0.0_dp

                    write(*, 24)
                    24 FORMAT( 'WARNING: denom_sctr is zero! ' )

                end if AvoidDivg

                FreeEng = FreeEng + ( 0.5_dp * Phiqq1(s1, s1, s, s) * mul_term )

            end do s1_loop

        end do s_loop

    end subroutine MultiplynBE_omega


    subroutine ScatterMatEl4( q1, q2, q3, q0chk, Nbasis, FC4, sys, &
                            & Evecq, Evecq1, Evecq2, Evecq3, Phi_ss1s2s3 )

        implicit none

        real(dp), dimension(3), intent(in)                          :: q1, q2, q3!, q
        logical, dimension(4), intent(in)                           :: q0chk
        integer, intent(in)                                         :: Nbasis
        type(FC4type), intent(in)                                   :: FC4
        type(cell), intent(in)                                      :: sys
        complex(dp), dimension(:, :), intent(in)                    :: Evecq, Evecq1, Evecq2, &
                                                                     & Evecq3

        complex(dp), dimension(:,:,:,:), intent(out)                :: Phi_ss1s2s3

        !================================= Local variables ===================================!

        complex(dp)                                 :: expnq1R1, expn_mass, phi4_mass_expn, &
                                                     & wqs_mu_alpha, wq1s1_nu_beta, wq2s2_eta_gama, wq3s3_ro_delta, &
                                                     & expnqR2_R3!, expnq2R2, expnq3R3

        real(dp), dimension(3)                      :: N1R, N2R, N3R
        real(dp)                                    :: q1R1, mass_fac, qR2_R3!, q2R2, q3R3

        integer                                     :: Natm4, Ndof, &
                                                     & rowno_mu_alpha, rowno_nu_beta, rowno_eta_gama, rowno_ro_delta
        integer                                     :: mu, N1, nu, N2, eta, N3, ro, &
                                                     & alpha, beta, gama, delta, &
                                                     & s, s1, s2, s3, ii
        integer, dimension(4)                       :: strt

        !================================= Local variables ===================================!

        Ndof = 3*Nbasis

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
                !expnq1R1 = cdexp( iu * q1R1 )

                expnq1R1 = dcmplx( dcos(q1R1), dsin(q1R1) )

                N2_loop: do N2 = 1, Natm4

                    eta = FC4%basisTyp(N2, mu)
                    N2R = FC4%atmPosR(:, N2, mu)

                    !-Off-! q2R2 = dot_product(q2, N2R)
                    !-Off-! expnq2R2 = dcmplx( dcos(q2R2), dsin(q2R2) )

                    !expnq2R2 = cdexp( iu * q2R2 )

                    N3_loop: do N3 = 1, Natm4

                        ro = FC4%basisTyp(N3, mu)
                        N3R = FC4%atmPosR(:, N3, mu)

                        !-Off-! q3R3 = dot_product(q3, N3R)
                        !-Off-! expnq3R3 = dcmplx( dcos(q3R3), dsin(q3R3) )

                        qR2_R3 = dot_product( q2, (N2R-N3R) )
                        expnqR2_R3 = dcmplx( dcos(qR2_R3), dsin(qR2_R3) )

                        mass_fac = 1.0_dp / dsqrt(  sys%mass(mu) *  sys%mass(nu) * sys%mass(eta) * sys%mass(ro) )

                        expn_mass = expnq1R1 * mass_fac * expnqR2_R3 !expnq2R2 * expnq3R3

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

                                            s2_loop: do s2 = strt(3), Ndof
                                                wq2s2_eta_gama = Evecq2(rowno_eta_gama, s2) 

                                                delta_loop: do delta = 1, 3
                                                    rowno_ro_delta = 3*(ro - 1) + (delta-1) + 1

                                                    phi4_mass_expn = FC4%fC(delta, gama, beta, alpha, N3, N2, N1, mu) * expn_mass

                                                    s3_loop: do s3 = strt(4), Ndof
                                                        wq3s3_ro_delta = Evecq3(rowno_ro_delta, s3)

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
    
end module FreeEnergy4th

