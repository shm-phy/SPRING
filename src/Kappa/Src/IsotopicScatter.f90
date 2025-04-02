
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

module IsoScatter

    use kinds,                  only : dp
    use constants,              only : PI
    use TetrahedronQ,           only : q_tetrahedron
    use phonon_m,               only : Phon
    use TetrahedronIntegration, only : FindFandOmegaAtVertices_Iso, &
                                     & FindgAndIver, FindIverxFver
    implicit none
    private

    public                  :: CalculateIsoW, CalculateIsoLwTetra

contains

    subroutine CalculateIsoW( ph_q0, ph_q1, i0, Nq1, &
                            & my_Qsize, my_offset, Nbasis, Ndof, g2iso, W_iso)

        implicit none

        type(Phon), intent(in)                                      :: ph_q0, ph_q1
        integer, intent(in)                                         :: i0, Nq1, my_Qsize, my_offset, &
                                                                     & Nbasis, Ndof
        real(dp), dimension(Nbasis), intent(in)                     :: g2iso !Explicit-shape dummy Array

        real(dp), dimension(Ndof, Ndof, Nq1), intent(inout)         :: W_iso !Explicit-shape dummy Array 

        ! ================================= Local Variables ================================= !

        complex(dp), dimension(3)           :: w_b_qs, w_b_q1s1
        complex(dp)                         :: wq0xwq1
        real(dp), dimension(2)              :: ReIm
        real(dp)                            :: Mq0s0, abs_wq0xwq1_2

        integer                             :: ii, q1_grid_indx, s0, &
                                             & s1, b, strt_cart, end_cart

        ! ================================= Local Variables ================================= !

        q1_loop: do ii = 1, my_Qsize

            q1_grid_indx = my_offset + ii

            s0_loop: do s0 = 1, Ndof

                Mq0s0 = ph_q0%omega(s0, i0)

                s1_loop: do s1 = 1, Ndof

                    b_loop: do b = 1, Nbasis

                        strt_cart = 3 * (b - 1) + 1
                        end_cart = strt_cart + 2

                        w_b_qs = ph_q0%Evec(strt_cart:end_cart, s0, i0) 
                        w_b_q1s1 = ph_q1%Evec(strt_cart:end_cart, s1, q1_grid_indx)

                        wq0xwq1 = dot_product( w_b_q1s1, w_b_qs )
                        ReIm = transfer( wq0xwq1, ReIm )
                        abs_wq0xwq1_2 = dot_product( ReIm, ReIm )

                        W_iso(s1, s0, q1_grid_indx) = W_iso(s1, s0, q1_grid_indx) + ( g2iso(b) * abs_wq0xwq1_2 )

                    end do b_loop

                    W_iso(s1, s0, q1_grid_indx) =  W_iso(s1, s0, q1_grid_indx) * PI * ( Mq0s0 ** 2 ) * 0.5_dp / dble( Nq1 )

                end do s1_loop

            end do s0_loop

        end do q1_loop

    end subroutine CalculateIsoW


    subroutine CalculateIsoLwTetra( Ndof, Nq1, my_QsizeT, my_offsetT, &
                                  & q1_tetra, ph_q0, ph_q1, &
                                  & W_s0s1All, i0, lw )

        implicit none

        integer, intent(in)                                     :: Ndof, Nq1, my_QsizeT, my_offsetT

        type(q_tetrahedron), intent(in)                         :: q1_tetra
        type(Phon), intent(in)                                  :: ph_q0
        type(Phon), intent(in)                                  :: ph_q1

        real(dp), dimension(Ndof, Ndof, Nq1), intent(in)        :: W_s0s1All

        integer, intent(in)                                     :: i0

        real(dp), dimension(Ndof), intent(out)                  :: lw

        ! ==================================== Local Variables ==================================== !

        real(dp), dimension(Ndof)                           :: omega_q0
        real(dp), dimension(4)                              :: F_ver, M_ver, Iver

        real(dp)                                            :: Mq0s, g, IxF

        integer, dimension(4)                               :: tetraVert_q1, &
                                                             & vertQ1_srt, vertQ2_srt

        integer                                             :: ii, nT, s0, s1

        logical                                             :: enterCase

        ! ==================================== Local Variables ==================================== !

        lw = 0.0_dp

        omega_q0 = ph_q0%omega(:, i0)

        s0_loop: do s0 = 1, Ndof

            Mq0s = omega_q0(s0)

            TetrahedronLoop: do ii = 1, my_QsizeT

                nT = my_offsetT + ii

                tetraVert_q1 = q1_tetra%tetrahedrons(:, nT)

                s1_loop: do s1 = 1, Ndof

                    call FindFandOmegaAtVertices_Iso( Ndof, Nq1, s0, s1, tetraVert_q1, &
                                                    & ph_q1, W_s0s1All, F_ver, M_ver )

                    ! ====-------------------------------------------------------------------------==== !
                    vertQ1_srt = tetraVert_q1
                    vertQ2_srt = tetraVert_q1
                    call FindgAndIver( Mq0s, g, Iver, M_ver, F_ver, vertQ1_srt, vertQ2_srt, enterCase )
                    !^! Avoid divergence due to Tetrahedron scheme !^!
                    WHERE( ISNAN(Iver) ) Iver = 0.0_dp
                    if ( ISNAN(g) ) g = 0.0_dp
                    !^! Avoid divergence due to Tetrahedron scheme !^!
                    IxF = FindIverxFver( Iver, F_ver, enterCase )
                    lw(s0) = lw(s0) + (g * IxF)
                    ! ====-------------------------------------------------------------------------==== !

                end do s1_loop

            end do TetrahedronLoop
        end do s0_loop

    end subroutine CalculateIsoLwTetra

end module IsoScatter

