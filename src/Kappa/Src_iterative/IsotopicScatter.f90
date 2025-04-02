
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

    use kinds,          only : dp
    use constants,      only : PI, cmp_prec, zero_prec, OneSixth
    use Irr_q_point,    only : q_points_data
    use phonon_m,       only : Phon

    implicit none

    EXTERNAL                :: dlasrt

    private

    public                  :: CalculateIsoW, TauIsoTetra, TetrahedronIterIso

contains

    subroutine CalculateIsoW( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                            & Nbasis, Ndof, Nqpoints, g2iso, W_iso )

        implicit none

        type(q_points_data), intent(in)                                     :: Qpoints
        type(Phon), intent(in)                                              :: phonon_dat
        integer, dimension(3), intent(in)                                   :: mesh
        integer, intent(in)                                                 :: my_Qsize, my_offset, Nbasis, &
                                                                             & Ndof, Nqpoints
        real(dp), dimension(Nbasis), intent(in)                             :: g2iso !Explicit-shape dummy Array

        real(dp), dimension(Ndof, Nqpoints, Ndof, my_Qsize), intent(inout)  :: W_iso !Explicit-shape dummy Array 

        ! ================================= Local Variables ================================= !

        complex(dp), dimension(3)           :: w_b_qs, w_b_q1s1
        complex(dp)                         :: wq0xwq1

        real(dp), dimension(2)              :: ReIm
        real(dp)                            :: Mqs, abs_wq0xwq1_2

        integer, dimension(3)               :: mesh_mul
        integer                             :: i, qindx, s
        integer, dimension(3)               :: q
        integer                             :: qPosIndx, i1, s1, b, strt_cart, end_cart

        ! ================================= Local Variables ================================= !

        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)

        irr_Q_loop: do i = 1, my_Qsize

            qindx = my_offset + i
            q = Qpoints%irr_q_int(1:3, qindx)

            qPosIndx = Qpoints%indx_map( dot_product( modulo( q, mesh ), mesh_mul ) + 1 )

            s_loop: do s = 1, Ndof

                Mqs = phonon_dat%omega(s, qPosIndx)

                q1_loop: do i1 = 1, Nqpoints

                    s1_loop: do s1 = 1, Ndof

                        b_loop: do b = 1, Nbasis

                            strt_cart = 3 * (b - 1) + 1
                            end_cart = strt_cart + 2

                            w_b_qs = phonon_dat%Evec(strt_cart:end_cart, s, qPosIndx) 

                            w_b_q1s1 = phonon_dat%Evec(strt_cart:end_cart, s1, i1)

                            wq0xwq1 = dot_product( w_b_q1s1, w_b_qs )
                            ReIm = transfer( wq0xwq1, ReIm )
                            abs_wq0xwq1_2 = dot_product( ReIm, ReIm )

                            W_iso(s1, i1, s, i) = W_iso(s1, i1, s, i) + ( g2iso(b) * abs_wq0xwq1_2 )

                        end do b_loop

                        W_iso(s1, i1, s, i) = W_iso(s1, i1, s, i) * ( Mqs ** 2 ) * PI * 0.5_dp / dble( Nqpoints )

                    end do s1_loop

                end do q1_loop

            end do s_loop

        end do irr_Q_loop

    end subroutine CalculateIsoW


    subroutine TauIsoTetra( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, Ndof, &
                          & Nqpoints, num_irr_q, W_iso, Inv_tauqsIso )

        implicit none

        type(q_points_data), intent(in)                                     :: Qpoints
        type(Phon), intent(in)                                              :: phonon_dat

        integer, dimension(3), intent(in)                                   :: mesh
        integer, intent(in)                                                 :: my_Qsize, my_offset, Ndof, Nqpoints, &
                                                                             & num_irr_q

        real(dp), dimension(Ndof, Nqpoints, Ndof, my_Qsize), intent(in)     :: W_iso !Explicit-shape dummy Array 

        real(dp), dimension(Ndof, num_irr_q), intent(inout)                 :: Inv_tauqsIso !Explicit-shape dummy Array

        ! =================================== Local Variables =================================== !

        real(dp), dimension(:), allocatable                 :: Mq
        real(dp), dimension(:,:), allocatable               :: Mvrtq1
        real(dp)                                            :: Mqs, g, tetra_val
        real(dp), dimension(4)                              :: Mq1s1_ver
        real(dp), dimension(4)                              :: Iver

        integer                                             :: NumTetra
        integer, dimension(3)                               :: mesh_mul, q_int
        integer                                             :: qi, qindx, qPosindx, NT

        integer, dimension(4)                               :: ver_q1, ver_q1_srt
        integer                                             :: ii, s, s1, kk, vv

        logical                                             :: enterCase

        ! =================================== Local Variables =================================== !

        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)
        NumTetra = size( Qpoints%tetrahedrons, 2 )

        allocate( Mq(Ndof) )
        allocate( Mvrtq1(Ndof, 4) )

        irr_Q: do qi = 1, my_Qsize

            qindx = my_offset + qi

            q_int = Qpoints%irr_q_int(1:3, qindx)
            qPosindx = Qpoints%indx_map( dot_product( modulo( q_int, mesh ), mesh_mul ) + 1 )

            Mq = phonon_dat%omega(:, qPosindx)

            TetraLoop: do NT = 1, NumTetra

                ver_q1 = Qpoints%tetrahedrons(:, NT)

                vert_loop1: do ii = 1, 4

                    !ver_q1(ii) = Qpoints%indx_map( ver_q1_unq(ii) )
                    ver_q1(ii) = Qpoints%indx_map( ver_q1(ii) )
                    Mvrtq1(:, ii) = phonon_dat%omega(:, ver_q1(ii))

                end do vert_loop1

                s_loop: do s = 1, Ndof
                    Mqs = Mq(s)

                    s1_loop: do s1 = 1, Ndof

                        vert_loop2: do kk = 1, 4

                            Mq1s1_ver(kk) = Mvrtq1(s1, kk)

                        end do vert_loop2

                        ver_q1_srt = ver_q1

                        call my_sort2(Mq1s1_ver, ver_q1_srt)

                        enterCase = .false.
                        tetra_val = 0.0_dp

                        tetra_coef_case: if ( (Mqs >= Mq1s1_ver(1)) .and. (Mqs < Mq1s1_ver(2)) ) then

                            ! ------------------- value of r and values of I at vertices ------------------- !
                            g = ( 3.0_dp * ( (Mqs - Mq1s1_ver(1)) ** 2 ) ) / &
                              & ( (Mq1s1_ver(2) - Mq1s1_ver(1)) * & 
                              &   (Mq1s1_ver(3) - Mq1s1_ver(1)) * &
                              &   (Mq1s1_ver(4) - Mq1s1_ver(1)) )

                            Iver(1) = ( ( (Mqs - Mq1s1_ver(2)) / (Mq1s1_ver(1) - Mq1s1_ver(2)) ) + &
                                      & ( (Mqs - Mq1s1_ver(3)) / (Mq1s1_ver(1) - Mq1s1_ver(3)) ) + &
                                      & ( (Mqs - Mq1s1_ver(4)) / (Mq1s1_ver(1) - Mq1s1_ver(4)) ) ) / 3.0_dp

                            do vv = 2, 4
                                Iver(vv) = ( (Mqs - Mq1s1_ver(1)) / (Mq1s1_ver(vv) - Mq1s1_ver(1)) ) / 3.0_dp
                            end do
                            ! ------------------- value of r and values of I at vertices ------------------- !

                            enterCase = .true.


                        else if ( (Mqs >= Mq1s1_ver(2)) .and. (Mqs < Mq1s1_ver(3)) ) then tetra_coef_case

                            ! ------------------- value of r and values of I at vertices ------------------- !
                            g = ( 3.0_dp / ( (Mq1s1_ver(2) - Mq1s1_ver(3)) * &
                              &              (Mq1s1_ver(4) - Mq1s1_ver(1)) ) ) * &
                              & ( ( ( (Mqs - Mq1s1_ver(1)) * (Mqs - Mq1s1_ver(3)) ) / &
                              &       (Mq1s1_ver(3) - Mq1s1_ver(1)) ) - &
                              &   ( ( (Mqs - Mq1s1_ver(2)) * (Mqs - Mq1s1_ver(4)) ) / &
                              &       (Mq1s1_ver(2) - Mq1s1_ver(4)) ) )

                            Iver(1) = ( ( (Mqs - Mq1s1_ver(4)) / (Mq1s1_ver(1) - Mq1s1_ver(4)) ) / 3.0_dp ) + &
                                    & ( ( ( (Mqs - Mq1s1_ver(3)) / (Mq1s1_ver(1) - Mq1s1_ver(3)) ) * &
                                    &   ( (Mqs - Mq1s1_ver(1)) / (Mq1s1_ver(3) - Mq1s1_ver(1)) ) * &
                                    &   ( (Mqs - Mq1s1_ver(3)) / (Mq1s1_ver(2) - Mq1s1_ver(3)) ) ) / &
                                    &   ( g * (Mq1s1_ver(4) - Mq1s1_ver(1)) ) )

                            Iver(2) = ( ( (Mqs - Mq1s1_ver(3)) / (Mq1s1_ver(2) - Mq1s1_ver(3)) ) / 3.0_dp ) + &
                                    & ( ( ( ((Mqs - Mq1s1_ver(4)) / (Mq1s1_ver(2) - Mq1s1_ver(4))) ** 2 ) * &
                                    &   ( (Mqs - Mq1s1_ver(2)) / (Mq1s1_ver(3) - Mq1s1_ver(2)) ) ) / &
                                    &   ( g * (Mq1s1_ver(4) - Mq1s1_ver(1)) ) )

                            Iver(3) = ( ( (Mqs - Mq1s1_ver(2)) / (Mq1s1_ver(3) - Mq1s1_ver(2)) ) / 3.0_dp ) + &
                                    & ( ( ( ((Mqs - Mq1s1_ver(1)) / (Mq1s1_ver(3) - Mq1s1_ver(1))) ** 2 ) * &
                                    &   ( (Mqs - Mq1s1_ver(3)) / (Mq1s1_ver(2) - Mq1s1_ver(3)) ) ) / &
                                    &   ( g * (Mq1s1_ver(4) - Mq1s1_ver(1)) ) )
                        
                            Iver(4) = ( ( (Mqs - Mq1s1_ver(1)) / (Mq1s1_ver(4) - Mq1s1_ver(1)) ) / 3.0_dp ) + &
                                    & ( ( ( (Mqs - Mq1s1_ver(2)) / (Mq1s1_ver(4) - Mq1s1_ver(2)) ) * &
                                    &   ( (Mqs - Mq1s1_ver(4)) / (Mq1s1_ver(2) - Mq1s1_ver(4)) ) * &
                                    &   ( (Mqs - Mq1s1_ver(2)) / (Mq1s1_ver(3) - Mq1s1_ver(2)) ) ) / &
                                    &   ( g * (Mq1s1_ver(4) - Mq1s1_ver(1)) ) )
                            ! ------------------- value of r and values of I at vertices ------------------- !

                            enterCase = .true.


                        else if ( (Mqs >= Mq1s1_ver(3)) .and. (Mqs <= Mq1s1_ver(4)) ) then tetra_coef_case

                            ! ------------------- value of r and values of I at vertices ------------------- !
                            g = ( 3.0_dp * ((Mqs - Mq1s1_ver(4)) ** 2) ) / &
                              & ( ( Mq1s1_ver(4) - Mq1s1_ver(1) ) * ( Mq1s1_ver(4) - Mq1s1_ver(2) ) * &
                              &   ( Mq1s1_ver(4) - Mq1s1_ver(3) ) )

                            do vv = 1, 3
                                Iver(vv) = ( (Mqs - Mq1s1_ver(4)) / (Mq1s1_ver(vv) - Mq1s1_ver(4)) ) / 3.0_dp
                            end do

                            Iver(4) = ( ( (Mqs - Mq1s1_ver(1)) / (Mq1s1_ver(4) - Mq1s1_ver(1)) ) + &
                                    &   ( (Mqs - Mq1s1_ver(2)) / (Mq1s1_ver(4) - Mq1s1_ver(2)) ) + &
                                    &   ( (Mqs - Mq1s1_ver(3)) / (Mq1s1_ver(4) - Mq1s1_ver(3)) ) ) / 3.0_dp
                            ! ------------------- value of r and values of I at vertices ------------------- !

                            enterCase = .true.

                        else tetra_coef_case

                            does_not_enters: if ( (Mqs >= Mq1s1_ver(1)) .and. &
                                                & (Mqs <= Mq1s1_ver(4)) .and. (.not. enterCase) ) then
                                
                                write(*, 55) Mqs
                                write(*, 65) Mq1s1_ver

                                55 FORMAT("WARNING: Mqs = ", F12.5, "THz. It does not enters any cases.")
                                65 FORMAT("         The values of Mq1s1_ver at vertices: [", 3(F12.5, ', '), F12.5, "] THz")

                            end if does_not_enters

                            enterCase = .false.

                        end if tetra_coef_case


                        chk_enterCase: if ( enterCase ) then

                            !^! Avoid divergence due to Tetrahedron scheme !^!
                            WHERE( ISNAN(Iver) ) Iver = 0.0_dp
                            if ( ISNAN(g) ) g = 0.0_dp
                            !^! Avoid divergence due to Tetrahedron scheme !^!

                            chk_sumeq1: if ( dabs(sum(Iver) - 1.0_dp) < zero_prec ) then

                                ! ----------------------------- sum over the 4 vertices ------------------------------ !
                                vert_loop3: do vv = 1, 4

                                    tetra_val = tetra_val + ( Iver(vv) * W_iso(s1, ver_q1_srt(vv), s, qi) )

                                end do vert_loop3
                                ! ----------------------------- sum over the 4 vertices ------------------------------ !
                                tetra_val = OneSixth * g * tetra_val

                                Inv_tauqsIso(s, qindx) = Inv_tauqsIso(s, qindx) + tetra_val ! ** !

                            else chk_sumeq1
                                
                                tetra_val = 0.0_dp
                                !-! write(*, 45) Iver

                            end if chk_sumeq1

                        else chk_enterCase

                            tetra_val = 0.0_dp

                        end if chk_enterCase

                    end do s1_loop

                end do s_loop

            end do TetraLoop

        end do irr_Q

        !-! 45 FORMAT("WARNING: sum( Iver(:) ) .ne. 1. Values of Iver(1:4): [", 3(F12.5, ', '), F12.5, "]")

        deallocate( Mq )
        deallocate( Mvrtq1 )

    end subroutine TauIsoTetra


    subroutine my_sort2(Omega_vert, tetra_ver_q1)

        implicit none

        real(dp), dimension(4), intent(inout)                       :: Omega_vert
        integer, dimension(4), intent(inout)                        :: tetra_ver_q1

        !! ===================================== Local variables =================================== !!

        real(dp), dimension(4)                      :: Omega_vertc
        integer, dimension(4)                       :: tetra_ver_q1c

        integer, dimension(1)                       :: indx
        integer, dimension(1)                       :: indx_back
        integer                                     :: info, ii, strt_pnt, ActualIndx

        logical, dimension(4)                       :: covered
        logical                                     :: back_chk

        !! ===================================== Local variables =================================== !!

        covered(:) = .false.
        back_chk = .false.

        Omega_vertc = Omega_vert
        tetra_ver_q1c = tetra_ver_q1

        call dlasrt( 'I', 4, Omega_vert, info )

        vert_loop: do ii = 1, 4

            strt_pnt = 0
            loop: do

                strt_pnt = mod( (strt_pnt + 1), 4)
                if (strt_pnt == 0) strt_pnt = 4

                indx = findloc( Omega_vertc(strt_pnt:4), Omega_vert(ii) )
                ActualIndx = strt_pnt + (indx(1) - 1)

                if ( .not. covered(ActualIndx) ) EXIT loop

            end do loop

            covered( ActualIndx ) = .true.

            tetra_ver_q1(ii) = tetra_ver_q1c( ActualIndx )

            debug: if ( back_chk ) then

                indx_back = findloc(Omega_vertc, Omega_vert(ii), back=back_chk)

                if ( ActualIndx /= indx_back(1) ) then

                    write(*, 80)
                    write(*, 90) Omega_vert

                    80 FORMAT("WARNING: {-/+ Omega(q1) + Omega(q2)} is same in two vertices of tetrahedron")
                    90 FORMAT("         The values of sorted Mq1q2s at vertices: [", 3(F12.5, ', '), F12.5, "] THz")

                end if

            end if debug

        end do vert_loop

    end subroutine my_sort2


    subroutine TetrahedronIterIso(Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                                & Ndof, Nqpoints, num_irr_q, W_iso, F_old, dFnew)

        implicit none

        type(q_points_data), intent(in)                                 :: Qpoints
        type(Phon), intent(in)                                          :: phonon_dat

        integer, dimension(3), intent(in)                               :: mesh
        integer, intent(in)                                             :: my_Qsize, my_offset, Ndof, Nqpoints, &
                                                                         & num_irr_q 

        real(dp), dimension(Ndof, Nqpoints, Ndof, my_Qsize), intent(in) :: W_iso !Explicit-shape dummy Array

        real(dp), dimension(3, Ndof, num_irr_q), intent(in)             :: F_old !Explicit-shape dummy Array

        real(dp), dimension(3, Ndof, num_irr_q), intent(inout)          :: dFnew !Explicit-shape dummy Array

        ! =================================== Local Variables =================================== !

        real(dp), dimension(:), allocatable                 :: Mq
        real(dp), dimension(:,:), allocatable               :: Mvrtq1
        real(dp)                                            :: Mqs, g
        real(dp), dimension(3)                              :: tetra_val
        real(dp), dimension(3)                              :: Flmb1
        real(dp), dimension(4)                              :: Mq1s1_ver
        real(dp), dimension(4)                              :: Iver

        integer                                             :: NumTetra
        integer, dimension(3)                               :: mesh_mul, q_int
        integer                                             :: qi, qindx, qPosindx, NT
        integer, dimension(4)                               :: ver_q1, ver_q1_srt
        integer                                             :: ii, s, s1, kk, vv
        integer                                             :: posInIrr_q1, existsInIrr_q1

        logical                                             :: enterCase

        ! =================================== Local Variables =================================== !

        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)
        NumTetra = size( Qpoints%tetrahedrons, 2 )

        allocate( Mq(Ndof) )
        allocate( Mvrtq1(Ndof, 4) )

        irr_Q: do qi = 1, my_Qsize

            qindx = my_offset + qi

            q_int = Qpoints%irr_q_int(1:3, qindx)
            qPosindx = Qpoints%indx_map( dot_product( modulo( q_int, mesh ), mesh_mul ) + 1 )

            Mq = phonon_dat%omega(:, qPosindx)

            TetraLoop: do NT = 1, NumTetra

                ver_q1 = Qpoints%tetrahedrons(:, NT)

                vert_loop1: do ii = 1, 4

                    !ver_q1(ii) = Qpoints%indx_map( ver_q1_unq(ii) )
                    ver_q1(ii) = Qpoints%indx_map( ver_q1(ii) )
                    Mvrtq1(:, ii) = phonon_dat%omega(:, ver_q1(ii))

                end do vert_loop1

                s_loop: do s = 1, Ndof
                    Mqs = Mq(s)

                    s1_loop: do s1 = 1, Ndof

                        vert_loop2: do kk = 1, 4

                            Mq1s1_ver(kk) = Mvrtq1(s1, kk)

                        end do vert_loop2

                        ver_q1_srt = ver_q1

                        call my_sort2(Mq1s1_ver, ver_q1_srt)

                        enterCase = .false.
                        tetra_val = 0.0_dp

                        tetra_coef_case: if ( (Mqs >= Mq1s1_ver(1)) .and. (Mqs < Mq1s1_ver(2)) ) then

                            ! ------------------- value of r and values of I at vertices ------------------- !
                            g = ( 3.0_dp * ( (Mqs - Mq1s1_ver(1)) ** 2 ) ) / &
                              & ( (Mq1s1_ver(2) - Mq1s1_ver(1)) * & 
                              &   (Mq1s1_ver(3) - Mq1s1_ver(1)) * &
                              &   (Mq1s1_ver(4) - Mq1s1_ver(1)) )

                            Iver(1) = ( ( (Mqs - Mq1s1_ver(2)) / (Mq1s1_ver(1) - Mq1s1_ver(2)) ) + &
                                      & ( (Mqs - Mq1s1_ver(3)) / (Mq1s1_ver(1) - Mq1s1_ver(3)) ) + &
                                      & ( (Mqs - Mq1s1_ver(4)) / (Mq1s1_ver(1) - Mq1s1_ver(4)) ) ) / 3.0_dp

                            do vv = 2, 4
                                Iver(vv) = ( (Mqs - Mq1s1_ver(1)) / (Mq1s1_ver(vv) - Mq1s1_ver(1)) ) / 3.0_dp
                            end do
                            ! ------------------- value of r and values of I at vertices ------------------- !

                            enterCase = .true.


                        else if ( (Mqs >= Mq1s1_ver(2)) .and. (Mqs < Mq1s1_ver(3)) ) then tetra_coef_case

                            ! ------------------- value of r and values of I at vertices ------------------- !
                            g = ( 3.0_dp / ( (Mq1s1_ver(2) - Mq1s1_ver(3)) * &
                              &              (Mq1s1_ver(4) - Mq1s1_ver(1)) ) ) * &
                              & ( ( ( (Mqs - Mq1s1_ver(1)) * (Mqs - Mq1s1_ver(3)) ) / &
                              &       (Mq1s1_ver(3) - Mq1s1_ver(1)) ) - &
                              &   ( ( (Mqs - Mq1s1_ver(2)) * (Mqs - Mq1s1_ver(4)) ) / &
                              &       (Mq1s1_ver(2) - Mq1s1_ver(4)) ) )

                            Iver(1) = ( ( (Mqs - Mq1s1_ver(4)) / (Mq1s1_ver(1) - Mq1s1_ver(4)) ) / 3.0_dp ) + &
                                    & ( ( ( (Mqs - Mq1s1_ver(3)) / (Mq1s1_ver(1) - Mq1s1_ver(3)) ) * &
                                    &   ( (Mqs - Mq1s1_ver(1)) / (Mq1s1_ver(3) - Mq1s1_ver(1)) ) * &
                                    &   ( (Mqs - Mq1s1_ver(3)) / (Mq1s1_ver(2) - Mq1s1_ver(3)) ) ) / &
                                    &   ( g * (Mq1s1_ver(4) - Mq1s1_ver(1)) ) )

                            Iver(2) = ( ( (Mqs - Mq1s1_ver(3)) / (Mq1s1_ver(2) - Mq1s1_ver(3)) ) / 3.0_dp ) + &
                                    & ( ( ( ((Mqs - Mq1s1_ver(4)) / (Mq1s1_ver(2) - Mq1s1_ver(4))) ** 2 ) * &
                                    &   ( (Mqs - Mq1s1_ver(2)) / (Mq1s1_ver(3) - Mq1s1_ver(2)) ) ) / &
                                    &   ( g * (Mq1s1_ver(4) - Mq1s1_ver(1)) ) )

                            Iver(3) = ( ( (Mqs - Mq1s1_ver(2)) / (Mq1s1_ver(3) - Mq1s1_ver(2)) ) / 3.0_dp ) + &
                                    & ( ( ( ((Mqs - Mq1s1_ver(1)) / (Mq1s1_ver(3) - Mq1s1_ver(1))) ** 2 ) * &
                                    &   ( (Mqs - Mq1s1_ver(3)) / (Mq1s1_ver(2) - Mq1s1_ver(3)) ) ) / &
                                    &   ( g * (Mq1s1_ver(4) - Mq1s1_ver(1)) ) )
                        
                            Iver(4) = ( ( (Mqs - Mq1s1_ver(1)) / (Mq1s1_ver(4) - Mq1s1_ver(1)) ) / 3.0_dp ) + &
                                    & ( ( ( (Mqs - Mq1s1_ver(2)) / (Mq1s1_ver(4) - Mq1s1_ver(2)) ) * &
                                    &   ( (Mqs - Mq1s1_ver(4)) / (Mq1s1_ver(2) - Mq1s1_ver(4)) ) * &
                                    &   ( (Mqs - Mq1s1_ver(2)) / (Mq1s1_ver(3) - Mq1s1_ver(2)) ) ) / &
                                    &   ( g * (Mq1s1_ver(4) - Mq1s1_ver(1)) ) )
                            ! ------------------- value of r and values of I at vertices ------------------- !

                            enterCase = .true.


                        else if ( (Mqs >= Mq1s1_ver(3)) .and. (Mqs <= Mq1s1_ver(4)) ) then tetra_coef_case

                            ! ------------------- value of r and values of I at vertices ------------------- !
                            g = ( 3.0_dp * ((Mqs - Mq1s1_ver(4)) ** 2) ) / &
                              & ( ( Mq1s1_ver(4) - Mq1s1_ver(1) ) * ( Mq1s1_ver(4) - Mq1s1_ver(2) ) * &
                              &   ( Mq1s1_ver(4) - Mq1s1_ver(3) ) )

                            do vv = 1, 3
                                Iver(vv) = ( (Mqs - Mq1s1_ver(4)) / (Mq1s1_ver(vv) - Mq1s1_ver(4)) ) / 3.0_dp
                            end do

                            Iver(4) = ( ( (Mqs - Mq1s1_ver(1)) / (Mq1s1_ver(4) - Mq1s1_ver(1)) ) + &
                                    &   ( (Mqs - Mq1s1_ver(2)) / (Mq1s1_ver(4) - Mq1s1_ver(2)) ) + &
                                    &   ( (Mqs - Mq1s1_ver(3)) / (Mq1s1_ver(4) - Mq1s1_ver(3)) ) ) / 3.0_dp
                            ! ------------------- value of r and values of I at vertices ------------------- !

                            enterCase = .true.

                        else tetra_coef_case

                            does_not_enters: if ( (Mqs >= Mq1s1_ver(1)) .and. &
                                                & (Mqs <= Mq1s1_ver(4)) .and. (.not. enterCase) ) then
                                
                                write(*, 55) Mqs
                                write(*, 65) Mq1s1_ver

                                55 FORMAT("WARNING: Mqs = ", F12.5, "THz. It does not enters any cases.")
                                65 FORMAT("         The values of Mq1s1_ver at vertices: [", 3(F12.5, ', '), F12.5, "] THz")

                            end if does_not_enters

                            enterCase = .false.

                        end if tetra_coef_case


                        chk_enterCase: if ( enterCase ) then

                            !^! Avoid divergence due to Tetrahedron scheme !^!
                            WHERE( ISNAN(Iver) ) Iver = 0.0_dp
                            if ( ISNAN(g) ) g = 0.0_dp
                            !^! Avoid divergence due to Tetrahedron scheme !^!

                            chk_sumeq1: if ( dabs(sum(Iver) - 1.0_dp) < zero_prec ) then

                                ! ----------------------------- sum over the 4 vertices ------------------------------ !
                                vert_loop3: do vv = 1, 4

                                    posInIrr_q1 = Qpoints%q_pnt_int( 4, ver_q1_srt(vv) )
                                    existsInIrr_q1 = Qpoints%q_pnt_int( 5, ver_q1_srt(vv) )

                                    DirectOrIndirect_q1: if ( existsInIrr_q1 == 1 ) then
                                        Flmb1 = F_old(:, s1, posInIrr_q1)

                                    else DirectOrIndirect_q1
                                        Flmb1 = -1.0_dp * F_old(:, s1, posInIrr_q1)

                                    end if DirectOrIndirect_q1

                                    tetra_val = tetra_val + ( Iver(vv) * W_iso(s1, ver_q1_srt(vv), s, qi) * Flmb1 )

                                end do vert_loop3
                                ! ----------------------------- sum over the 4 vertices ------------------------------ !
                                tetra_val = OneSixth * g * tetra_val

                                dFnew(:, s, qindx) = dFnew(:, s, qindx) + tetra_val ! ** !

                            else chk_sumeq1
                                
                                tetra_val = 0.0_dp
                                !-! write(*, 45) Iver

                            end if chk_sumeq1

                        else chk_enterCase

                            tetra_val = 0.0_dp

                        end if chk_enterCase

                    end do s1_loop

                end do s_loop

            end do TetraLoop

        end do irr_Q

        !-! 45 FORMAT("WARNING: sum( Iver(:) ) .ne. 1. Values of Iver(1:4): [", 3(F12.5, ', '), F12.5, "]")

        deallocate( Mq )
        deallocate( Mvrtq1 )

    end subroutine TetrahedronIterIso

end module IsoScatter

