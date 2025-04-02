
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

module Tetrahedron4th

    use kinds,                  only : dp
    use TetrahedronQ,           only : q_tetrahedron
    use phonon_m,               only : Phon
    use constants,              only : OnebySix

    use TetrahedronIntegration, only : FindgAndIver, &
                                     & FindIverxFver

    implicit none
    private
    public              :: TetrahedronIntegration4th

contains

    subroutine TetrahedronIntegration4th( Nq1, Ndof, QPerProcsMax, Nq1q2Perm, my_Qsize, my_offset, i0, &
                                        & q1_tetra, ph_q0, ph_q1, Sctrq1q2, Polarization4th, Allq1q2q3, &
                                        & QinWhichProcs, lw )

        implicit none
        
        integer, intent(in)                                                                 :: Nq1, Ndof, &
                                                                                             & QPerProcsMax, &
                                                                                             & Nq1q2Perm, my_Qsize, &
                                                                                             & my_offset, i0

        type(q_tetrahedron), intent(in)                                                     :: q1_tetra
        type(Phon), intent(in)                                                              :: ph_q0
        type(Phon), intent(in)                                                              :: ph_q1

        real(dp), dimension(Ndof,Ndof,Ndof,Ndof, QPerProcsMax), codimension[*], intent(in)  :: Sctrq1q2

        integer, dimension(2, 4), intent(in)                                                :: Polarization4th
        integer, dimension(3, Nq1, Nq1), intent(in)                                         :: Allq1q2q3
        !integer, dimension(3, Nq1q2Perm), intent(in)                                        :: q1q2q3_permuted
        integer, dimension(2, Nq1q2Perm), intent(in)                                        :: QinWhichProcs

        real(dp), dimension(Ndof), intent(out)                                              :: lw


        ! =========================================== Local Variables =========================================== !

        real(dp), dimension(Ndof)                       :: omega_q0, nBE_q0, omega_q1, nBE_q1

        real(dp), dimension(Ndof,Ndof,Ndof,Ndof)        :: SctrProb_ver1, SctrProb_ver2, &
                                                         & SctrProb_ver3, SctrProb_ver4

        integer, dimension(3)                           :: q0_int
        integer, dimension(4)                           :: q2_grid_indx_ver, q3_grid_indx_ver, &
                                                         & type_of_perm_ver, pos_per_sctr_ver

        integer                                         :: ii, q1_grid_indx, nT, NTetrahedrons, &
                                                         & vv

        integer                                         :: ImageNo_ver1, ImageNo_ver2, ImageNo_ver3, ImageNo_ver4, &
                                                         & pos_ver1, pos_ver2, pos_ver3, pos_ver4

        logical                                         :: q0IsZero

        ! =========================================== Local Variables =========================================== !

        lw(:) = 0.0_dp

        NTetrahedrons = 6 * Nq1

        q0_int(:) = q1_tetra%irr_q_pnt_int(:, i0)
        q0IsZero = all( q0_int == 0 )
        omega_q0(:) = ph_q0%omega(:, i0)
        nBE_q0(:) = ph_q0%nBE(:, i0)

        q1_loop: do ii = 1, my_Qsize
            
            q1_grid_indx = my_offset + ii

            omega_q1(:) = ph_q1%omega(:, q1_grid_indx)
            nBE_q1(:) = ph_q1%nBE(:, q1_grid_indx)

            Tetrahedron: do nT = 1, NTetrahedrons

                q2_grid_indx_ver(:) = q1_tetra%tetrahedrons(:, nT)

                vert_loop: do vv = 1, 4

                    type_of_perm_ver(vv) = Allq1q2q3(1, q2_grid_indx_ver(vv), q1_grid_indx)
                    pos_per_sctr_ver(vv) = Allq1q2q3(2, q2_grid_indx_ver(vv), q1_grid_indx)
                    q3_grid_indx_ver(vv) = Allq1q2q3(3, q2_grid_indx_ver(vv), q1_grid_indx)

                end do vert_loop

                ! ----------------------------- ====================== ----------------------------- !
                ImageNo_ver1 = QinWhichProcs( 1, pos_per_sctr_ver(1) )
                pos_ver1 = QinWhichProcs( 2, pos_per_sctr_ver(1) )

                if ( ImageNo_ver1 /= this_image() ) then
                    SctrProb_ver1(:, :, :, :) = Sctrq1q2(:, :, :, :, pos_ver1) [ImageNo_ver1]

                else
                    SctrProb_ver1(:, :, :, :) = Sctrq1q2(:, :, :, :, pos_ver1)

                end if
                ! ----------------------------- ====================== ----------------------------- !

                ! ----------------------------- ====================== ----------------------------- !
                ImageNo_ver2 = QinWhichProcs( 1, pos_per_sctr_ver(2) )
                pos_ver2 = QinWhichProcs( 2, pos_per_sctr_ver(2) )

                if ( ImageNo_ver2 /= this_image() ) then
                    SctrProb_ver2(:, :, :, :) = Sctrq1q2(:, :, :, :, pos_ver2) [ImageNo_ver2]

                else
                    SctrProb_ver2(:, :, :, :) = Sctrq1q2(:, :, :, :, pos_ver2)

                end if
                ! ----------------------------- ====================== ----------------------------- !

                ! ----------------------------- ====================== ----------------------------- !
                ImageNo_ver3 = QinWhichProcs( 1, pos_per_sctr_ver(3) )
                pos_ver3 = QinWhichProcs( 2, pos_per_sctr_ver(3) )

                if ( ImageNo_ver3 /= this_image() ) then
                    SctrProb_ver3(:, :, :, :) = Sctrq1q2(:, :, :, :, pos_ver3) [ImageNo_ver3]

                else
                    SctrProb_ver3(:, :, :, :) = Sctrq1q2(:, :, :, :, pos_ver3)

                end if
                ! ----------------------------- ====================== ----------------------------- !

                ! ----------------------------- ====================== ----------------------------- !
                ImageNo_ver4 = QinWhichProcs( 1, pos_per_sctr_ver(4) )
                pos_ver4 = QinWhichProcs( 2, pos_per_sctr_ver(4) )

                if ( ImageNo_ver4 /= this_image() ) then
                    SctrProb_ver4(:, :, :, :) = Sctrq1q2(:, :, :, :, pos_ver4) [ImageNo_ver4]

                else
                    SctrProb_ver4(:, :, :, :) = Sctrq1q2(:, :, :, :, pos_ver4)

                end if
                ! ----------------------------- ====================== ----------------------------- !

                call SumOverPolarization( Ndof, ph_q1, SctrProb_ver1, SctrProb_ver2, SctrProb_ver3, SctrProb_ver4, &
                                        & omega_q0, nBE_q0, omega_q1, nBE_q1, Polarization4th, type_of_perm_ver, &
                                        & q2_grid_indx_ver, q3_grid_indx_ver, q0IsZero, lw)

            end do Tetrahedron

        end do q1_loop

    end subroutine TetrahedronIntegration4th


    subroutine SumOverPolarization( Ndof, ph_q1, SctrProb_ver1, SctrProb_ver2, SctrProb_ver3, SctrProb_ver4, &
                                  & omega_q0, nBE_q0, omega_q1, nBE_q1, Polarization4th, type_of_perm_ver, &
                                  & q2_grid_indx_ver, q3_grid_indx_ver, q0IsZero, lw )

        implicit none

        integer, intent(in)                                     :: Ndof

        type(Phon), intent(in)                                  :: ph_q1

        real(dp), dimension(Ndof,Ndof,Ndof,Ndof), intent(in)    :: SctrProb_ver1, SctrProb_ver2, &
                                                                 & SctrProb_ver3, SctrProb_ver4
        real(dp), dimension(Ndof), intent(in)                   :: omega_q0, nBE_q0, omega_q1, nBE_q1

        !integer, intent(in)                                     :: i0, q1_grid_indx
        integer, dimension(2, 4), intent(in)                    :: Polarization4th
        integer, dimension(4), intent(in)                       :: type_of_perm_ver, q2_grid_indx_ver, &
                                                                 & q3_grid_indx_ver

        logical, intent(in)                                     :: q0IsZero
        real(dp), dimension(Ndof), intent(inout)                :: lw

        ! ===================================== Local Variables ===================================== !

        real(dp), dimension(4)                  :: F1_ver, F2_ver, F3_ver, &
                                                 & M1_ver, M2_ver, M3_ver, &
                                                 & I1_ver, I2_ver, I3_ver

        real(dp)                                :: Mq0s0, nBEq0s0, Mq1s1, nBEq1s1, &
                                                 & g1, g2, g3, IxF1, IxF2, IxF3

        integer, dimension(4)                   :: vertQ2_srt, vertQ3_srt
        integer                                 :: s0, s1, s2, s3, s0Strt

        logical                                 :: enterCase1, enterCase2, enterCase3

        ! ===================================== Local Variables ===================================== !

        s0Strt = Polarization4th(1, 1)
        if ( q0IsZero ) s0Strt = 4

        s0_loop: do s0 = s0Strt, Polarization4th(2, 1)

            Mq0s0 = omega_q0(s0)
            nBEq0s0 = nBE_q0(s0)

            s1_loop: do s1 = Polarization4th(1, 2), Polarization4th(2, 2)

                Mq1s1 = omega_q1(s1)
                nBEq1s1 = nBE_q1(s1)

                s2_loop: do s2 = Polarization4th(1, 3), Polarization4th(2, 3)
                    s3_loop: do s3 = Polarization4th(1, 4), Polarization4th(2, 4)

                        call FindFandOmegaVer_4th( Ndof, ph_q1, nBEq0s0, Mq1s1, nBEq1s1, SctrProb_ver1, SctrProb_ver2, &
                                                 & SctrProb_ver3, SctrProb_ver4, type_of_perm_ver, q2_grid_indx_ver, &
                                                 & q3_grid_indx_ver, s0, s1, s2, s3, F1_ver, F2_ver, F3_ver, &
                                                 & M1_ver, M2_ver, M3_ver )

                        ! ---------------------------- ============================ ---------------------------- !
                        vertQ2_srt(:) = q2_grid_indx_ver(:)
                        vertQ3_srt(:) = q3_grid_indx_ver(:)
                        call FindgAndIver( Mq0s0, g1, I1_ver, M1_ver, F1_ver, vertQ2_srt, vertQ3_srt, enterCase1 )
                        !^! Avoid divergence due to Tetrahedron scheme (due to <= and >= conditions) !^!
                        WHERE( ISNAN(I1_ver) ) I1_ver = 0.0_dp
                        if ( ISNAN(g1) ) g1 = 0.0_dp
                        !^! Avoid divergence due to Tetrahedron scheme (due to <= and >= conditions) !^!
                        IxF1 = FindIverxFver( I1_ver, F1_ver, enterCase1 )
                        lw(s0) = lw(s0) + ( OnebySix * g1 * IxF1 )
                        ! ---------------------------- ============================ ---------------------------- !
                        
                        ! ---------------------------- ============================ ---------------------------- !
                        vertQ2_srt(:) = q2_grid_indx_ver(:)
                        vertQ3_srt(:) = q3_grid_indx_ver(:)
                        call FindgAndIver( Mq0s0, g2, I2_ver, M2_ver, F2_ver, vertQ2_srt, vertQ3_srt, enterCase2 )
                        !^! Avoid divergence due to Tetrahedron scheme (due to <= and >= conditions) !^!
                        WHERE( ISNAN(I2_ver) ) I2_ver = 0.0_dp
                        if ( ISNAN(g2) ) g2 = 0.0_dp
                        !^! Avoid divergence due to Tetrahedron scheme (due to <= and >= conditions) !^!
                        IxF2 = FindIverxFver( I2_ver, F2_ver, enterCase2 )
                        lw(s0) = lw(s0) + ( 0.5_dp * g2 * IxF2 )
                        ! ---------------------------- ============================ ---------------------------- !

                        ! ---------------------------- ============================ ---------------------------- !
                        vertQ2_srt(:) = q2_grid_indx_ver(:)
                        vertQ3_srt(:) = q3_grid_indx_ver(:)
                        call FindgAndIver( Mq0s0, g3, I3_ver, M3_ver, F3_ver, vertQ2_srt, vertQ3_srt, enterCase3 )
                        !^! Avoid divergence due to Tetrahedron scheme (due to <= and >= conditions) !^!
                        WHERE( ISNAN(I3_ver) ) I3_ver = 0.0_dp
                        if ( ISNAN(g3) ) g3 = 0.0_dp
                        !^! Avoid divergence due to Tetrahedron scheme (due to <= and >= conditions) !^!
                        IxF3 = FindIverxFver( I3_ver, F3_ver, enterCase3 )
                        lw(s0) = lw(s0) + ( 0.5_dp * g3 * IxF3 )
                        ! ---------------------------- ============================ ---------------------------- !

                    end do s3_loop
                end do s2_loop
            end do s1_loop
        end do s0_loop

    end subroutine SumOverPolarization


    subroutine FindFandOmegaVer_4th( Ndof, ph_q1, nBEq0s0, Mq1s1, nBEq1s1, SctrProb_ver1, SctrProb_ver2, &
                                   & SctrProb_ver3, SctrProb_ver4, type_of_perm_ver, q2_grid_indx_ver, &
                                   & q3_grid_indx_ver, s0, s1, s2, s3, F1_ver, F2_ver, F3_ver, &
                                   & M1_ver, M2_ver, M3_ver )

        implicit none

        integer, intent(in)                                     :: Ndof

        type(Phon), intent(in)                                  :: ph_q1

        real(dp), intent(in)                                    :: nBEq0s0, Mq1s1, nBEq1s1

        real(dp), dimension(Ndof,Ndof,Ndof,Ndof), intent(in)    :: SctrProb_ver1, SctrProb_ver2, &
                                                                 & SctrProb_ver3, SctrProb_ver4
        
        integer, dimension(4), intent(in)                       :: type_of_perm_ver, q2_grid_indx_ver, &
                                                                 & q3_grid_indx_ver
        integer, intent(in)                                     :: s0, s1, s2, s3

        real(dp), dimension(4), intent(out)                     :: F1_ver, F2_ver, F3_ver, &
                                                                 & M1_ver, M2_ver, M3_ver

        ! ================================ Local Variables ================================ !

        integer, dimension(3)                       :: before_perm

        integer, dimension(3)                       :: perm_ver1, perm_ver2, perm_ver3, perm_ver4, &
                                                     & after_perm_ver1, after_perm_ver2, &
                                                     & after_perm_ver3, after_perm_ver4

        integer, dimension(3, 5)                    :: per

        ! ================================ Local Variables ================================ !

        before_perm(:) = (/s1, s2, s3/)

        per(:, 1) = (/1, 3, 2/)
        per(:, 2) = (/2, 1, 3/)
        per(:, 3) = (/2, 3, 1/)
        per(:, 4) = (/3, 1, 2/)
        per(:, 5) = (/3, 2, 1/)

        ! ---------------------------------- ================================= ---------------------------------- !
        if ( type_of_perm_ver(1) == 0 ) then
            perm_ver1(:) = (/1, 2, 3/)
        else
            perm_ver1(:) = per(:, type_of_perm_ver(1))
        end if

        after_perm_ver1(:) = (/before_perm(perm_ver1(1)), before_perm(perm_ver1(2)), before_perm(perm_ver1(3))/)

        F1_ver(1) = SctrProb_ver1( after_perm_ver1(3), after_perm_ver1(2), after_perm_ver1(1), s0 ) * &
                  & ( nBEq1s1 * ph_q1%nBE(s2, q2_grid_indx_ver(1)) * ph_q1%nBE(s3, q3_grid_indx_ver(1)) / nBEq0s0 )

        M1_ver(1) = ( Mq1s1 + ph_q1%omega(s2, q2_grid_indx_ver(1)) + ph_q1%omega(s3, q3_grid_indx_ver(1)) )

        F2_ver(1) = SctrProb_ver1( after_perm_ver1(3), after_perm_ver1(2), after_perm_ver1(1), s0 ) * &
                  & ( (1.0_dp + nBEq1s1) * ph_q1%nBE(s2, q2_grid_indx_ver(1)) * ph_q1%nBE(s3, q3_grid_indx_ver(1)) / nBEq0s0 )

        M2_ver(1) = ( ph_q1%omega(s2, q2_grid_indx_ver(1)) + ph_q1%omega(s3, q3_grid_indx_ver(1)) - Mq1s1 )

        F3_ver(1) = SctrProb_ver1( after_perm_ver1(3), after_perm_ver1(2), after_perm_ver1(1), s0 ) * &
                  & ( (1.0_dp + nBEq1s1) * (1.0_dp + ph_q1%nBE(s2, q2_grid_indx_ver(1))) * ph_q1%nBE(s3, q3_grid_indx_ver(1)) &
                  & / nBEq0s0 )

        M3_ver(1) = ( ph_q1%omega(s3, q3_grid_indx_ver(1)) - Mq1s1 - ph_q1%omega(s2, q2_grid_indx_ver(1)) )
        ! ---------------------------------- ================================= ---------------------------------- !

        ! ---------------------------------- ================================= ---------------------------------- !
        if ( type_of_perm_ver(2) == 0 ) then
            perm_ver2(:) = (/1, 2, 3/)
        else
            perm_ver2(:) = per(:, type_of_perm_ver(2))
        end if

        after_perm_ver2(:) = (/before_perm(perm_ver2(1)), before_perm(perm_ver2(2)), before_perm(perm_ver2(3))/)

        F1_ver(2) = SctrProb_ver2( after_perm_ver2(3), after_perm_ver2(2), after_perm_ver2(1), s0 ) * &
                  & ( nBEq1s1 * ph_q1%nBE(s2, q2_grid_indx_ver(2)) * ph_q1%nBE(s3, q3_grid_indx_ver(2)) / nBEq0s0 )

        M1_ver(2) = ( Mq1s1 + ph_q1%omega(s2, q2_grid_indx_ver(2)) + ph_q1%omega(s3, q3_grid_indx_ver(2)) )

        F2_ver(2) = SctrProb_ver2( after_perm_ver2(3), after_perm_ver2(2), after_perm_ver2(1), s0 ) * &
                  & ( (1.0_dp + nBEq1s1) * ph_q1%nBE(s2, q2_grid_indx_ver(2)) * ph_q1%nBE(s3, q3_grid_indx_ver(2)) / nBEq0s0 )

        M2_ver(2) = ( ph_q1%omega(s2, q2_grid_indx_ver(2)) + ph_q1%omega(s3, q3_grid_indx_ver(2)) - Mq1s1 )

        F3_ver(2) = SctrProb_ver2( after_perm_ver2(3), after_perm_ver2(2), after_perm_ver2(1), s0 ) * &
                  & ( (1.0_dp + nBEq1s1) * (1.0_dp + ph_q1%nBE(s2, q2_grid_indx_ver(2))) * ph_q1%nBE(s3, q3_grid_indx_ver(2)) &
                  & / nBEq0s0 )

        M3_ver(2) = ( ph_q1%omega(s3, q3_grid_indx_ver(2)) - Mq1s1 - ph_q1%omega(s2, q2_grid_indx_ver(2)) )
        ! ---------------------------------- ================================= ---------------------------------- !

        ! ---------------------------------- ================================= ---------------------------------- !
        if ( type_of_perm_ver(3) == 0 ) then
            perm_ver3(:) = (/1, 2, 3/)
        else
            perm_ver3(:) = per(:, type_of_perm_ver(3))
        end if

        after_perm_ver3(:) = (/before_perm(perm_ver3(1)), before_perm(perm_ver3(2)), before_perm(perm_ver3(3))/)

        F1_ver(3) = SctrProb_ver3( after_perm_ver3(3), after_perm_ver3(2), after_perm_ver3(1), s0 ) * &
                  & ( nBEq1s1 * ph_q1%nBE(s2, q2_grid_indx_ver(3)) * ph_q1%nBE(s3, q3_grid_indx_ver(3)) / nBEq0s0 )

        M1_ver(3) = ( Mq1s1 + ph_q1%omega(s2, q2_grid_indx_ver(3)) + ph_q1%omega(s3, q3_grid_indx_ver(3)) )

        F2_ver(3) = SctrProb_ver3( after_perm_ver3(3), after_perm_ver3(2), after_perm_ver3(1), s0 ) * &
                  & ( (1.0_dp + nBEq1s1) * ph_q1%nBE(s2, q2_grid_indx_ver(3)) * ph_q1%nBE(s3, q3_grid_indx_ver(3)) / nBEq0s0 )

        M2_ver(3) = ( ph_q1%omega(s2, q2_grid_indx_ver(3)) + ph_q1%omega(s3, q3_grid_indx_ver(3)) - Mq1s1 )

        F3_ver(3) = SctrProb_ver3( after_perm_ver3(3), after_perm_ver3(2), after_perm_ver3(1), s0 ) * &
                  & ( (1.0_dp + nBEq1s1) * (1.0_dp + ph_q1%nBE(s2, q2_grid_indx_ver(3))) * ph_q1%nBE(s3, q3_grid_indx_ver(3)) &
                  & / nBEq0s0 )

        M3_ver(3) = ( ph_q1%omega(s3, q3_grid_indx_ver(3)) - Mq1s1 - ph_q1%omega(s2, q2_grid_indx_ver(3)) )
        ! ---------------------------------- ================================= ---------------------------------- !

        ! ---------------------------------- ================================= ---------------------------------- !
        if ( type_of_perm_ver(4) == 0 ) then
            perm_ver4(:) = (/1, 2, 3/)
        else
            perm_ver4(:) = per(:, type_of_perm_ver(4))
        end if

        after_perm_ver4(:) = (/before_perm(perm_ver4(1)), before_perm(perm_ver4(2)), before_perm(perm_ver4(3))/)

        F1_ver(4) = SctrProb_ver4( after_perm_ver4(3), after_perm_ver4(2), after_perm_ver4(1), s0 ) * &
                  & ( nBEq1s1 * ph_q1%nBE(s2, q2_grid_indx_ver(4)) * ph_q1%nBE(s3, q3_grid_indx_ver(4)) / nBEq0s0 )

        M1_ver(4) = ( Mq1s1 + ph_q1%omega(s2, q2_grid_indx_ver(4)) + ph_q1%omega(s3, q3_grid_indx_ver(4)) )

        F2_ver(4) = SctrProb_ver4( after_perm_ver4(3), after_perm_ver4(2), after_perm_ver4(1), s0 ) * &
                  & ( (1.0_dp + nBEq1s1) * ph_q1%nBE(s2, q2_grid_indx_ver(4)) * ph_q1%nBE(s3, q3_grid_indx_ver(4)) / nBEq0s0 )

        M2_ver(4) = ( ph_q1%omega(s2, q2_grid_indx_ver(4)) + ph_q1%omega(s3, q3_grid_indx_ver(4)) - Mq1s1 )

        F3_ver(4) = SctrProb_ver4( after_perm_ver4(3), after_perm_ver4(2), after_perm_ver4(1), s0 ) * &
                  & ( (1.0_dp + nBEq1s1) * (1.0_dp + ph_q1%nBE(s2, q2_grid_indx_ver(4))) * ph_q1%nBE(s3, q3_grid_indx_ver(4)) &
                  & / nBEq0s0 )

        M3_ver(4) = ( ph_q1%omega(s3, q3_grid_indx_ver(4)) - Mq1s1 - ph_q1%omega(s2, q2_grid_indx_ver(4)) )
        ! ---------------------------------- ================================= ---------------------------------- !


    end subroutine FindFandOmegaVer_4th

end module Tetrahedron4th

