
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

module ThirdOrder

   use kinds,                   only : dp
   use constants,               only : Sctr3THz2, Onebysqrt2PI, &
                                     & print_step, zero_prec
   use unit_cell,               only : cell

   use FC3_mod,                 only : FC3type
   use phonon_m,                only : Phon
   use TetrahedronQ,            only : q_tetrahedron
   use TetrahedronIntegration,  only : FindFandOmegaAtVertices, FindgAndIver, &
                                     & FindIverxFver
   implicit none
   private

   public                   :: ScatterProbAllq1, CalculateLinewidth, &
                             & CalculateLwTetrahedron

contains

    subroutine ScatterProbAllq1( sys, FC3, ph_q0, ph_q1, q_tetra, q1mesh, &
                               & i0, Nq1points, Nq1Perm, Nbasis, Ndof, my_Qsize, &
                               & my_offset, q1q2_permuted, Sctrq1 )

        implicit none

        type(cell), intent(in)                                          :: sys
        type(FC3type), intent(in)                                       :: FC3
        type(Phon), intent(in)                                          :: ph_q0
        type(Phon), intent(in)                                          :: ph_q1
        type(q_tetrahedron), intent(in)                                 :: q_tetra

        integer, dimension(3), intent(in)                               :: q1mesh
        integer, intent(in)                                             :: i0, Nq1points, Nq1Perm, &
                                                                         & Nbasis, Ndof, &
                                                                         & my_Qsize, my_offset

        !integer, dimension(4, Nq1points), intent(in)                    :: Allq1q2
        integer, dimension(2, Nq1Perm), intent(in)                      :: q1q2_permuted

        real(dp), dimension(Ndof,Ndof,Ndof, Nq1points), intent(inout)   :: Sctrq1

        !================================= Local variables ===================================!

        complex(dp), dimension(Ndof, Ndof)                          :: Evecq0, Evecq1, Evecq2
        complex(dp), dimension(Ndof, Ndof, Ndof)                    :: Phi_q0q1q2

        real(dp), dimension(3)                                      :: q1_float, q2_float

        real(dp), dimension(Ndof)                                   :: OnebyOmega_q0, OnebyOmega_q1, &
                                                                     & OnebyOmega_q2
        real(dp), dimension(Ndof,Ndof,Ndof)                         :: W_s0s1s2, W_s0s1s2_per
        real(dp)                                                    :: multiply_factor

        integer                                                     :: ii, i1

        integer                                                     :: q1_grid_indx, q2_grid_indx

        integer, dimension(3)                                       :: mesh_mul, q0_int, q1_int, &
                                                                     & q2_int

        logical, dimension(3)                                       :: q0chk

        !================================= Local variables ===================================!

        mesh_mul = (/1, q1mesh(1), q1mesh(1)*q1mesh(2)/)

        multiply_factor = ( Sctr3THz2 / dble(Nq1points) )

        ! ------------------------------------------ q0 ------------------------------------------ !
        q0_int = q_tetra%irr_q_pnt_int(:, i0)

        Evecq0(:, :) = ph_q0%Evec(:, :, i0)
        OnebyOmega_q0(:) = ph_q0%OnebyOmega(:, i0)
        q0chk(1) = all( q0_int == 0 )
        ! ------------------------------------------ q0 ------------------------------------------ !

        q1_loop: do ii = 1, my_Qsize

            i1 = my_offset + ii     ! ** !

            ! ------------------------------------------ q1 ------------------------------------------ !
            q1_grid_indx = q1q2_permuted(1, i1)

            q1_int = q_tetra%q_pnt_int(1:3, q1_grid_indx)
            q1_float = q_tetra%q_pnt(1:3, q1_grid_indx)

            Evecq1(:, :) = ph_q1%Evec(:, :, q1_grid_indx)
            OnebyOmega_q1(:) = ph_q1%OnebyOmega(:, q1_grid_indx)
            q0chk(2) = all( q1_int == 0 )
            ! ------------------------------------------ q1 ------------------------------------------ !

            ! ------------------------------------------ q2 ------------------------------------------ !
            q2_grid_indx = q1q2_permuted(2, i1)

            q2_int = q_tetra%q_pnt_int(1:3, q2_grid_indx)
            q2_float = q_tetra%q_pnt(1:3, q2_grid_indx)

            Evecq2(:, :) = ph_q1%Evec(:, :, q2_grid_indx)
            OnebyOmega_q2(:) = ph_q1%OnebyOmega(:, q2_grid_indx)
            q0chk(3) = all( q2_int == 0 )
            ! ------------------------------------------ q2 ------------------------------------------ !

            ! ---------------------------- Scattering Matrix Calculation ----------------------------- !

            call ScatterMatEl( q1_float, q2_float, q0chk, Nbasis, Ndof, FC3, &
                             & sys, Evecq0, Evecq1, Evecq2, Phi_q0q1q2 )

            call ScatterProbability( Ndof, q0chk, OnebyOmega_q0, OnebyOmega_q1, OnebyOmega_q2, &
                                   & Phi_q0q1q2, multiply_factor, W_s0s1s2, W_s0s1s2_per )

            ! ---------------------------- Scattering Matrix Calculation ----------------------------- !

            Sctrq1(:, :, :, q1_grid_indx) = W_s0s1s2(:, :, :)

#ifdef __DEBUG
            debug_chk: if ( (q1_grid_indx == q2_grid_indx) .and. &
                          & any( (q2_int-q1_int) /= 0 ) ) then
                            
                write(*, 230) q1_grid_indx, q2_grid_indx, q1_int, q2_int
                ERROR STOP

                230 FORMAT( "ERROR: in ScatterProbAllq1. q1_grid_indx = ", I5, &
                          & "q2_grid_indx = ", I5, ". q1_grid_address = ( ", 2(I4, ' '), I4, &
                          & " ) and q1_grid_address = ( ", 2(I4, ' '), I4, " )" )

            end if debug_chk
#endif

            write_permutation: if ( q1_grid_indx /= q2_grid_indx ) then

#ifdef __DEBUG
                debug2: if ( .not. all(dabs(Sctrq1(:, :, :, q2_grid_indx)) < zero_prec) ) then
                    write(*, 232)
                    232 FORMAT( "ERROR in ScatterProbAllq1!" )
                    ERROR STOP
                end if debug2
#endif

                Sctrq1(:, :, :, q2_grid_indx) = W_s0s1s2_per(:, :, :)

            end if write_permutation

        end do q1_loop

    end subroutine ScatterProbAllq1


    subroutine CalculateLwTetrahedron( Ndof, Nq1, my_QsizeT, my_offsetT, &
                                     & q1_tetra, ph_q0, ph_q1, W_s0s1s2All, &
                                     & Polarization3rd, Allq1q2, i0, P3, lw, P3_q0 )

        implicit none

        integer, intent(in)                                     :: Ndof, Nq1, my_QsizeT, my_offsetT

        type(q_tetrahedron), intent(in)                         :: q1_tetra
        type(Phon), intent(in)                                  :: ph_q0
        type(Phon), intent(in)                                  :: ph_q1

        real(dp), dimension(Ndof, Ndof, Ndof, Nq1), intent(in)  :: W_s0s1s2All

        integer, dimension(2, 3), intent(in)                    :: Polarization3rd
        integer, dimension(4, Nq1), intent(in)                  :: Allq1q2
        integer, intent(in)                                     :: i0

        real(dp), intent(out)                                   :: P3
        real(dp), dimension(Ndof), intent(out)                  :: lw, P3_q0

        ! ==================================== Local Variables ==================================== !

        real(dp), dimension(Ndof)                           :: omega_q0
        real(dp), dimension(4)                              :: F_ver1, F_ver2, F_ver3, &
                                                             & M_ver1, M_ver2, M_ver3, &
                                                             & Iver

        real(dp)                                            :: Mq0s, g, IxF

        integer, dimension(4)                               :: tetraVert_q1, tetraVert_q2,&
                                                             & vertQ1_srt, vertQ2_srt

        integer                                             :: ii, nT, vv, s0, s1, s2

        logical                                             :: enterCase

        ! ==================================== Local Variables ==================================== !

        P3 = 0.0_dp
        lw = 0.0_dp
        P3_q0 = 0.0_dp

        omega_q0 = ph_q0%omega(:, i0)

        s0_loop: do s0 = Polarization3rd(1, 1), Polarization3rd(2, 1)

            Mq0s = omega_q0(s0)

            TetrahedronLoop: do ii = 1, my_QsizeT

                nT = my_offsetT + ii

                tetraVert_q1 = q1_tetra%tetrahedrons(:, nT)

                vertices: do vv = 1, 4
                    tetraVert_q2(vv) = Allq1q2( 4, tetraVert_q1(vv) )
                end do vertices

                s1_loop: do s1 = Polarization3rd(1, 2), Polarization3rd(2, 2)

                    s2_loop: do s2 = Polarization3rd(1, 3), Polarization3rd(2, 3)

                        call FindFandOmegaAtVertices( Ndof, Nq1, s0, s1, s2, &
                                                    & tetraVert_q1, tetraVert_q2, &
                                                    & ph_q1, W_s0s1s2All, &
                                                    & F_ver1, F_ver2, M_ver1, M_ver2, M_ver3)

                        F_ver3 = F_ver2

                        ! ====-------------------------------------------------------------------------==== !
                        vertQ1_srt = tetraVert_q1
                        vertQ2_srt = tetraVert_q2
                        call FindgAndIver( Mq0s, g, Iver, M_ver1, F_ver1, vertQ1_srt, vertQ2_srt, enterCase )
                        IxF = FindIverxFver( Iver, F_ver1, enterCase )
                        lw(s0) = lw(s0) + (g * IxF)
                        P3_q0(s0) = P3_q0(s0) + g
                        P3 = P3 + g
                        ! ====-------------------------------------------------------------------------==== !

                        ! ====-------------------------------------------------------------------------==== !
                        vertQ1_srt = tetraVert_q1
                        vertQ2_srt = tetraVert_q2
                        call FindgAndIver( Mq0s, g, Iver, M_ver2, F_ver2, vertQ1_srt, vertQ2_srt, enterCase )
                        IxF = FindIverxFver( Iver, F_ver2, enterCase )
                        lw(s0) = lw(s0) + (g * IxF)
                        P3_q0(s0) = P3_q0(s0) + g
                        P3 = P3 + g
                        ! ====-------------------------------------------------------------------------==== !

                        ! ====-------------------------------------------------------------------------==== !
                        vertQ1_srt = tetraVert_q1
                        vertQ2_srt = tetraVert_q2
                        call FindgAndIver( Mq0s, g, Iver, M_ver3, F_ver3, vertQ1_srt, vertQ2_srt, enterCase )
                        IxF = FindIverxFver( Iver, F_ver3, enterCase )
                        lw(s0) = lw(s0) - (g * IxF)
                        P3_q0(s0) = P3_q0(s0) + g ! Important => it is + ??
                        P3 = P3 + g ! Important => it is + ??
                        ! ====-------------------------------------------------------------------------==== !

                    end do s2_loop
                end do s1_loop

            end do TetrahedronLoop
        end do s0_loop

    end subroutine CalculateLwTetrahedron


    subroutine ScatterMatEl(q1, q2, q0chk, Nbasis, Ndof, FC3, sys, &
                          & Evecq, Evecq1, Evecq2, Phi_ss1s2)

        implicit none

        real(dp), dimension(3), intent(in)                          :: q1, q2!, q
        logical, dimension(3), intent(in)                           :: q0chk
        integer, intent(in)                                         :: Nbasis, Ndof
        type(FC3type), intent(in)                                   :: FC3
        type(cell), intent(in)                                      :: sys
        complex(dp), dimension(Ndof, Ndof), intent(in)              :: Evecq, Evecq1, Evecq2

        complex(dp), dimension(Ndof,Ndof,Ndof), intent(out)         :: Phi_ss1s2

        !================================= Local variables ===================================!

        complex(dp)                                 :: expnq1R1, expnq2R2, expn_mass, phi3_mass_expn, &
                                                     & wqs_mu_alpha, wq1s1_nu_beta, wq2s2_eta_gama

        real(dp), dimension(3)                      :: N1R, N2R
        real(dp)                                    :: q1R1, q2R2, mass_fac

        integer                                     :: rowno_mu_alpha, rowno_nu_beta, rowno_eta_gama, &
                                                     & Natm3
        integer                                     :: mu, N1, nu, N2, eta, &
                                                     & alpha, beta, gama, &
                                                     & s, s1, s2, ii
        integer, dimension(3)                       :: strt

        !================================= Local variables ===================================!

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

                    !mass_fac = 1.0_dp / dsqrt(  sys%mass(mu) *  sys%mass(nu) * sys%mass(eta) )
                    mass_fac = ( sys%OnebySqrtMass(mu) * sys%OnebySqrtMass(nu) * sys%OnebySqrtMass(eta) )

                    expn_mass = expnq1R1 * expnq2R2 * mass_fac

                    alpha_loop: do alpha = 1, 3
                        rowno_mu_alpha = 3*(mu-1) + (alpha-1) + 1

                        s_loop: do s = strt(1), Ndof
                            wqs_mu_alpha = Evecq(rowno_mu_alpha, s) !ToDo: make this memory access contiguous

                            beta_loop: do beta = 1, 3
                                rowno_nu_beta = 3*(nu-1) + (beta-1) + 1

                                s1_loop: do s1 = strt(2), Ndof
                                    wq1s1_nu_beta = Evecq1(rowno_nu_beta, s1) !ToDo: make this memory access contiguous

                                    gama_loop: do gama = 1, 3
                                        rowno_eta_gama = 3*(eta-1) + (gama-1) + 1

                                        phi3_mass_expn = FC3%fC(gama, beta, alpha, N2, N1, mu) * expn_mass

                                        s2_loop: do s2 = strt(3), Ndof
                                            wq2s2_eta_gama = Evecq2(rowno_eta_gama, s2) !ToDo: make this memory access contiguous
                                            
                                            Phi_ss1s2(s2, s1, s) = Phi_ss1s2(s2, s1, s) + &
                                          & ( phi3_mass_expn * wqs_mu_alpha * wq1s1_nu_beta * wq2s2_eta_gama )

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


    subroutine ScatterProbability( Ndof, q0chk, OnebyOmega_q, OnebyOmega_q1, OnebyOmega_q2, &
                                 & Phi_ss1s2, multiply_factor, W_ss1s2, W_ss1s2_per )

        implicit none

        integer, intent(in)                                         :: Ndof
        logical, dimension(3), intent(in)                           :: q0chk
        real(dp), dimension(Ndof), intent(in)                       :: OnebyOmega_q, OnebyOmega_q1, &
                                                                     & OnebyOmega_q2
        complex(dp), dimension(Ndof,Ndof,Ndof), intent(in)          :: Phi_ss1s2

        real(dp), intent(in)                                        :: multiply_factor
        real(dp), dimension(Ndof,Ndof,Ndof), intent(out)            :: W_ss1s2, W_ss1s2_per

        !! ===================================== Local variables =================================== !!

        real(dp), dimension(2)                      :: ReIm
        real(dp)                                    :: abs2, elmnt

        integer                                     :: ii, s, s1, s2
        integer, dimension(3)                       :: strt

        !! ===================================== Local variables =================================== !!

        strt(:) = 1
        do ii = 1, 3
            if ( q0chk(ii) ) strt(ii) = 4
        end do

        W_ss1s2 = 0.0_dp
        W_ss1s2_per = 0.0_dp

        s_loop2: do s = strt(1), Ndof
            s1_loop2: do s1 = strt(2), Ndof
                s2_loop2: do s2 = strt(3), Ndof

                    ReIm = transfer( Phi_ss1s2(s2, s1, s), ReIm )
                    abs2 = dot_product( ReIm, ReIm )

                    !elmnt = abs2 / ( omega_q(s) * omega_q1(s1) * omega_q2(s2) * dble(Nq) )
                    elmnt = abs2 * ( OnebyOmega_q(s) * OnebyOmega_q1(s1) * OnebyOmega_q2(s2) )

                    W_ss1s2(s2, s1, s) = elmnt * multiply_factor !Sctr3THz2

                    W_ss1s2_per(s1, s2, s) = elmnt * multiply_factor !Sctr3THz2

                end do s2_loop2
            end do s1_loop2
        end do s_loop2

    end subroutine ScatterProbability

    ! ** Obsolete **!
    subroutine CalculateLinewidth( Ndof, Nq1, my_QsizeAllq, my_offsetAllq, &
                                 & sys, q1_tetra, ph_q0, ph_q1, W_s0s1s2, &
                                 & OnebySigma, OnebySigma2, i0, q1mesh, Allq1q2, lw )

        implicit none

        integer, intent(in)                                     :: Ndof, Nq1, my_QsizeAllq, my_offsetAllq

        type(cell), intent(in)                                  :: sys
        type(q_tetrahedron), intent(in)                         :: q1_tetra
        type(Phon), intent(in)                                  :: ph_q0
        type(Phon), intent(in)                                  :: ph_q1

        real(dp), dimension(Ndof, Ndof, Ndof, Nq1), intent(in)  :: W_s0s1s2
        real(dp), intent(in)                                    :: OnebySigma, OnebySigma2

        integer, intent(in)                                     :: i0
        integer, dimension(3), intent(in)                       :: q1mesh
        integer, dimension(4, Nq1), intent(in)                  :: Allq1q2

        real(dp), dimension(Ndof), intent(out)                  :: lw

        ! ==================================== Local Variables ==================================== !

        real(dp)                                            :: DelOmega1, DelOmega2, DelOmega3, &
                                                             & delta1, delta2, delta3
        real(dp), dimension(Ndof)                           :: lw_tmp
        real(dp), dimension(Ndof)                           :: omega_q0, omega_q1, omega_q2, &
                                                             & nBE_q1, nBE_q2
        ! ** For adaptive Sigma ** !
        real(dp)                                            :: adaptiveSigma
        real(dp), dimension(3)                              :: grp_velq1s1, grp_velq2s2
        ! ** For adaptive Sigma ** !
        integer, dimension(3)                               :: strt, q1_int, q2_int, mesh_mul
        integer                                             :: ii, q1i, s0, s1, s2, q1_unq_pos, &
                                                             & q2_unq_pos

        ! ==================================== Local Variables ==================================== !

        lw = 0.0_dp

        strt(:) = 1
        mesh_mul = (/1, q1mesh(1), q1mesh(1)*q1mesh(2)/)

        ! ============---------------------------------------------------============ !
        omega_q0 = ph_q0%omega(:, i0)
        ! ============---------------------------------------------------============ !

        q1_loop: do ii = 1, my_QsizeAllq

            q1i = my_offsetAllq + ii   ! ** !

            ! ============---------------------------------------------------============ !
            q1_int = q1_tetra%q_pnt_int(1:3, q1i)
            q1_unq_pos = dot_product( modulo( q1_int, q1mesh ), mesh_mul ) + 1

            omega_q1 = ph_q1%omega(:, q1i)
            nBE_q1 = ph_q1%nBE(:, q1i)
            if (all( q1_int == 0 )) strt(2) = 4
            ! ============---------------------------------------------------============ !

            ! ============---------------------------------------------------============ !
            q2_unq_pos = Allq1q2( 4, q1_unq_pos )
            q2_int = q1_tetra%q_pnt_int(1:3, q2_unq_pos)

            omega_q2 = ph_q1%omega(:, q2_unq_pos)
            nBE_q2 = ph_q1%nBE(:, q2_unq_pos)
            if ( all(q2_int == 0) ) strt(3) = 4
            ! ============---------------------------------------------------============ !

            lw_tmp(:) = 0.0_dp

            s0_loop: do s0 = 1, Ndof
                s1_loop: do s1 = strt(2), Ndof
                    grp_velq1s1 = ph_q1%grp_vel(1:3, s1, q1_unq_pos) !* Check *!

                    s2_loop: do s2 = strt(3), Ndof
                        grp_velq2s2 = ph_q1%grp_vel(1:3, s2, q2_unq_pos) !* Check *!

                        DelOmega1 = omega_q0(s0) - omega_q1(s1) - omega_q2(s2)
                        DelOmega2 = omega_q0(s0) + omega_q1(s1) - omega_q2(s2)
                        DelOmega3 = omega_q0(s0) - omega_q1(s1) + omega_q2(s2)

                        adaptiveSigma = FindSigma( q1mesh, sys%G, grp_velq1s1, grp_velq2s2 )

                        delta1 = OnebySigma * Onebysqrt2PI * dexp( -0.5_dp * OnebySigma2 * (DelOmega1**2) )
                        delta2 = OnebySigma * Onebysqrt2PI * dexp( -0.5_dp * OnebySigma2 * (DelOmega2**2) )
                        delta3 = OnebySigma * Onebysqrt2PI * dexp( -0.5_dp * OnebySigma2 * (DelOmega3**2) )

                        lw_tmp(s0) = lw_tmp(s0) + W_s0s1s2(s2, s1, s0, q1_unq_pos) * &
                                   & ( ( (nBE_q1(s1) + nBE_q2(s2) + 1.0_dp) * delta1 ) + &
                                   &   ( (nBE_q1(s1) - nBE_q2(s2)) * (delta2 - delta3) ) )

                    end do s2_loop
                end do s1_loop
            end do s0_loop

            lw = lw + lw_tmp

        end do q1_loop

    end subroutine CalculateLinewidth
    ! ** Obsolete **!

    ! ** Obsolete **!
    Function FindSigma( mesh, G, grp_velq1s1, grp_velq2s2 ) Result ( Sigma )

        implicit none

        integer, dimension(3), intent(in)               :: mesh
        real(dp), dimension(3, 3), intent(in)           :: G
        real(dp), dimension(3), intent(in)              :: grp_velq1s1, grp_velq2s2

        real(dp)                                        :: Sigma

        ! ==================================== Local variables ==================================== !

        real(dp)                                        :: tmp1, tmp2

        integer                                         :: mu, alpha

        ! ==================================== Local variables ==================================== !

        tmp2 = 0.0_dp

        mu_loop: do mu = 1, 3

            tmp1 = 0.0_dp
            alpha_loop: do alpha = 1, 3

                tmp1 = tmp1 + (grp_velq1s1(alpha) - grp_velq2s2(alpha)) * G(alpha, mu)

            end do alpha_loop

            tmp2 = tmp2 + ( tmp1 / dble(mesh(mu)) ) **2

        end do mu_loop

        Sigma = 0.288675135_dp * dsqrt( tmp2 ) * 0.01_dp

    end Function FindSigma
    ! ** Obsolete **!

end module ThirdOrder

