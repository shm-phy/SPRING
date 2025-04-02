
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


subroutine TetrahedronIterPlus(this, Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                             & Ndof, num_irr_q, F_old, dFnew)

    implicit none

    class(AllQ), intent(in)                                     :: this

    type(q_points_data), intent(in)                             :: Qpoints
    type(Phon), intent(in)                                      :: phonon_dat

    integer, dimension(3), intent(in)                           :: mesh
    integer, intent(in)                                         :: my_Qsize, my_offset, &
                                                                 & Ndof, num_irr_q

    real(dp), dimension(3, Ndof, num_irr_q), intent(in)         :: F_old

    real(dp), dimension(3, Ndof, num_irr_q), intent(inout)      :: dFnew    

    !================================= Local variable ===================================!

    integer                                             :: NumTetra

    integer                                             :: NT, qi, qindx, qindx2, &
                                                         & s, ii, s1, s2, kk, vv

    integer, dimension(3)                               :: mesh_mul, q_int
    logical, dimension(3)                               :: q0chk

    real(dp), dimension(:), allocatable                 :: Mq

    real(dp)                                            :: Mqs, g
    real(dp), dimension(3)                              :: tetra_val
    logical                                             :: enterCase

    real(dp), dimension(4)                              :: Iver

    integer, dimension(4)                               :: ver_q1, ver_q2, sctrPos, &
                                                         & ver_q1_srt, ver_q2_srt, sctrPos_srt !,ver_q1_unq

    logical, dimension(4)                               :: drctFlag, drctFlag_srt

    real(dp), dimension(:,:), allocatable               :: Mvrtq1, Mvrtq2

    real(dp), dimension(4)                              :: Mq1q2s

    integer                                             :: posInIrr_q1, existsInIrr_q1, &
                                                         & posInIrr_q2, existsInIrr_q2

    real(dp), dimension(3)                              :: Flmb1, Flmb2, Fl2_minus_l1

    !================================= Local variable ===================================!

    mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)

    NumTetra = size( Qpoints%tetrahedrons, 2 )

    allocate( Mq(Ndof) )
    allocate( Mvrtq1(Ndof, 4) )
    allocate( Mvrtq2(Ndof, 4) )

    irr_q_loop: do qi = 1, my_Qsize

        qindx = my_offset + qi

        q_int = Qpoints%irr_q_int(1:3, qindx)
        qindx2 = Qpoints%indx_map( dot_product( modulo( q_int, mesh ), mesh_mul ) + 1 )

        Mq = phonon_dat%omega(:, qindx2)

        q0chk(1) = all( q_int == 0 )

        !-!s_loop: do s = 1, Ndof

        !-!    Mqs = Mq(s)

        loop_Tetra: do NT = 1, NumTetra

            !ver_q1_unq = Qpoints%tetrahedrons(:, NT)
            ver_q1 = Qpoints%tetrahedrons(:, NT)

            vert_loop1: do ii = 1, 4

                !ver_q1(ii) = Qpoints%indx_map( ver_q1_unq(ii) )
                ver_q1(ii) = Qpoints%indx_map( ver_q1(ii) )

                Mvrtq1(:, ii) = phonon_dat%omega(:, ver_q1(ii))

                sctrPos(ii) = this%ConsrvPlus(qi)%ConsrvQ( ver_q1(ii) )%SctrElPos

                drctFlag(ii) = this%ConsrvPlus(qi)%ConsrvQ( ver_q1(ii) )%flag

                ver_q2(ii) = this%ConsrvPlus(qi)%ConsrvQ( ver_q1(ii) )%q2indx
                Mvrtq2(:, ii) = phonon_dat%omega(:, ver_q2(ii))

            end do vert_loop1

            s_loop: do s = 1, Ndof
                Mqs = Mq(s)

                s1_loop: do s1 = 1, Ndof
                    s2_loop: do s2 = 1, Ndof

                        vert_loop2: do kk = 1, 4

                            !Mq1q2s(kk) = ( -1.0_dp * Mvrtq1(s1, kk) ) + Mvrtq2(s2, kk)

                            Mq1q2s(kk) = ( Mvrtq2(s2, kk) - Mvrtq1(s1, kk) )

                        end do vert_loop2

                        ver_q1_srt = ver_q1
                        ver_q2_srt = ver_q2
                        sctrPos_srt = sctrPos
                        drctFlag_srt = drctFlag

                        call my_sort(Mq1q2s, ver_q1_srt, ver_q2_srt, sctrPos_srt, drctFlag_srt)

                        enterCase = .false.
                        tetra_val = 0.0_dp

                        tetra_coef_case: if ( (Mqs > Mq1q2s(1)) .and. (Mqs < Mq1q2s(2)) ) then

                            ! ------------------- value of r and values of I at vertices ------------------- !
                            g = ( 3.0_dp * ( (Mqs - Mq1q2s(1)) ** 2 ) ) / &
                              & ( (Mq1q2s(2) - Mq1q2s(1)) * & 
                              &   (Mq1q2s(3) - Mq1q2s(1)) * &
                              &   (Mq1q2s(4) - Mq1q2s(1)) )

                            Iver(1) = ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(1) - Mq1q2s(2)) ) + &
                                      & ( (Mqs - Mq1q2s(3)) / (Mq1q2s(1) - Mq1q2s(3)) ) + &
                                      & ( (Mqs - Mq1q2s(4)) / (Mq1q2s(1) - Mq1q2s(4)) ) ) / 3.0_dp

                            do vv = 2, 4
                                Iver(vv) = ( (Mqs - Mq1q2s(1)) / (Mq1q2s(vv) - Mq1q2s(1)) ) / 3.0_dp
                            end do
                            ! ------------------- value of r and values of I at vertices ------------------- !

                            enterCase = .true.


                        else if ( (Mqs > Mq1q2s(2)) .and. (Mqs < Mq1q2s(3)) ) then tetra_coef_case

                            ! ------------------- value of r and values of I at vertices ------------------- !
                            g = ( 3.0_dp / ( (Mq1q2s(2) - Mq1q2s(3)) * &
                              &              (Mq1q2s(4) - Mq1q2s(1)) ) ) * &
                              & ( ( ( (Mqs - Mq1q2s(1)) * (Mqs - Mq1q2s(3)) ) / &
                              &       (Mq1q2s(3) - Mq1q2s(1)) ) - &
                              &   ( ( (Mqs - Mq1q2s(2)) * (Mqs - Mq1q2s(4)) ) / &
                              &       (Mq1q2s(2) - Mq1q2s(4)) ) )

                            Iver(1) = ( ( (Mqs - Mq1q2s(4)) / (Mq1q2s(1) - Mq1q2s(4)) ) / 3.0_dp ) + &
                                    & ( ( ( (Mqs - Mq1q2s(3)) / (Mq1q2s(1) - Mq1q2s(3)) ) * &
                                    &   ( (Mqs - Mq1q2s(1)) / (Mq1q2s(3) - Mq1q2s(1)) ) * &
                                    &   ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) ) / &
                                    &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )

                            Iver(2) = ( ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) / 3.0_dp ) + &
                                    & ( ( ( ((Mqs - Mq1q2s(4)) / (Mq1q2s(2) - Mq1q2s(4))) ** 2 ) * &
                                    &   ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) ) / &
                                    &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )

                            Iver(3) = ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) / 3.0_dp ) + &
                                    & ( ( ( ((Mqs - Mq1q2s(1)) / (Mq1q2s(3) - Mq1q2s(1))) ** 2 ) * &
                                    &   ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) ) / &
                                    &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )
                        
                            Iver(4) = ( ( (Mqs - Mq1q2s(1)) / (Mq1q2s(4) - Mq1q2s(1)) ) / 3.0_dp ) + &
                                    & ( ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(4) - Mq1q2s(2)) ) * &
                                    &   ( (Mqs - Mq1q2s(4)) / (Mq1q2s(2) - Mq1q2s(4)) ) * &
                                    &   ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) ) / &
                                    &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )
                            ! ------------------- value of r and values of I at vertices ------------------- !

                            enterCase = .true.


                        else if ( (Mqs > Mq1q2s(3)) .and. (Mqs < Mq1q2s(4)) ) then tetra_coef_case

                            ! ------------------- value of r and values of I at vertices ------------------- !
                            g = ( 3.0_dp * ((Mqs - Mq1q2s(4)) ** 2) ) / &
                              & ( ( Mq1q2s(4) - Mq1q2s(1) ) * ( Mq1q2s(4) - Mq1q2s(2) ) * &
                              &   ( Mq1q2s(4) - Mq1q2s(3) ) )

                            do vv = 1, 3
                                Iver(vv) = ( (Mqs - Mq1q2s(4)) / (Mq1q2s(vv) - Mq1q2s(4)) ) / 3.0_dp
                            end do

                            Iver(4) = ( ( (Mqs - Mq1q2s(1)) / (Mq1q2s(4) - Mq1q2s(1)) ) + &
                                    &   ( (Mqs - Mq1q2s(2)) / (Mq1q2s(4) - Mq1q2s(2)) ) + &
                                    &   ( (Mqs - Mq1q2s(3)) / (Mq1q2s(4) - Mq1q2s(3)) ) ) / 3.0_dp
                            ! ------------------- value of r and values of I at vertices ------------------- !

                            enterCase = .true.

                        else tetra_coef_case

                            does_not_enters: if ( (Mqs >= Mq1q2s(1)) .and. &
                                                & (Mqs <= Mq1q2s(4)) .and. (.not. enterCase) ) then
                                
                                write(*, 55) Mqs
                                write(*, 65) Mq1q2s

                                55 FORMAT("WARNING: Mqs = ", F12.5, "THz. It does not enters any cases.")
                                65 FORMAT("         The values of Mq1q2s at vertices: [", 3(F12.5, ', '), F12.5, "] THz")

                            end if does_not_enters

                            enterCase = .false.

                        end if tetra_coef_case


                        chk_enterCase: if ( enterCase ) then

                            chk_sumeq1: if ( dabs(sum(Iver) - 1.0_dp) < zero_prec ) then

                                ! ----------------------------- sum over the 4 vertices ------------------------------ !
                                vert_loop3: do vv = 1, 4

                                    !  What is l1 and l2 at vertices vv ( let's say vv = 1), and find F(l2) - F(l1)  !

                                    posInIrr_q1 = Qpoints%q_pnt_int( 4, ver_q1_srt(vv) )
                                    existsInIrr_q1 = Qpoints%q_pnt_int( 5, ver_q1_srt(vv) )

                                    DirectOrIndirect_q1: if ( existsInIrr_q1 == 1 ) then
                                        Flmb1 = F_old(:, s1, posInIrr_q1)

                                    else DirectOrIndirect_q1
                                        Flmb1 = -1.0_dp * F_old(:, s1, posInIrr_q1)

                                    end if DirectOrIndirect_q1

                                    posInIrr_q2 = Qpoints%q_pnt_int( 4, ver_q2_srt(vv) )
                                    existsInIrr_q2 = Qpoints%q_pnt_int( 5, ver_q2_srt(vv) )

                                    DirectOrIndirect_q2: if ( existsInIrr_q2 == 1 ) then
                                        Flmb2 = F_old(:, s2, posInIrr_q2)

                                    else DirectOrIndirect_q2
                                        Flmb2 = -1.0_dp * F_old(:, s2, posInIrr_q2)

                                    end if DirectOrIndirect_q2

                                    Fl2_minus_l1 = Flmb2 - Flmb1

                                    per_symm: if ( drctFlag_srt(vv) ) then

                                        tetra_val = tetra_val + ( (Iver(vv) * &
                                                  & this%WPlus(qi)%q1q2( sctrPos_srt(vv) )%ss1s2(s2, s1, s)) * &
                                                  & Fl2_minus_l1 )

                                    else per_symm

                                        !*! tetra_val = tetra_val + ( (Iver(vv) * &
                                        !*!           & (-1.0_dp) * this%WPlus(qi)%q1q2( sctrPos_srt(vv) )%ss1s2(s1, s2, s)) * &
                                        !*!           & Fl2_minus_l1 )

                                        tetra_val = tetra_val - ( (Iver(vv) * &
                                                  & this%WPlus(qi)%q1q2( sctrPos_srt(vv) )%ss1s2(s1, s2, s)) * &
                                                  & Fl2_minus_l1 )

                                    end if per_symm

                                end do vert_loop3
                                ! ----------------------------- sum over the 4 vertices ------------------------------ !
                                tetra_val = OneSixth * g * tetra_val

                                dFnew(:, s, qindx) = dFnew(:, s, qindx) + tetra_val ! ** !

                            else chk_sumeq1
                                
                                tetra_val = 0.0_dp
                                write(*, 45) Iver

                            end if chk_sumeq1

                        else chk_enterCase

                            tetra_val = 0.0_dp

                        end if chk_enterCase

                    end do s2_loop
                end do s1_loop
            end do s_loop

        end do loop_Tetra

    end do irr_q_loop

    45 FORMAT("WARNING: sum( Iver(:) ) .ne. 1. Values of Iver(1:4): [", 3(F12.5, ', '), F12.5, "]")

    deallocate( Mq )
    deallocate( Mvrtq1, Mvrtq2 )

end subroutine TetrahedronIterPlus


subroutine TetrahedronIterMinus(this, Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                              & Ndof, num_irr_q, F_old, dFnew)

    implicit none

    class(AllQ), intent(in)                                     :: this

    type(q_points_data), intent(in)                             :: Qpoints
    type(Phon), intent(in)                                      :: phonon_dat

    integer, dimension(3), intent(in)                           :: mesh
    integer, intent(in)                                         :: my_Qsize, my_offset, &
                                                                 & Ndof, num_irr_q

    real(dp), dimension(3, Ndof, num_irr_q), intent(in)         :: F_old

    real(dp), dimension(3, Ndof, num_irr_q), intent(inout)      :: dFnew    


    !================================= Local variable ===================================!

    integer                                             :: NumTetra

    integer                                             :: NT, qi, qindx, qindx2, &
                                                         & s, ii, s1, s2, kk, vv

    integer, dimension(3)                               :: mesh_mul, q_int
    logical, dimension(3)                               :: q0chk

    real(dp), dimension(:), allocatable                 :: Mq

    real(dp)                                            :: Mqs, g
    real(dp), dimension(3)                              :: tetra_val
    logical                                             :: enterCase

    real(dp), dimension(4)                              :: Iver

    integer, dimension(4)                               :: ver_q1, ver_q2, sctrPos, &
                                                         & ver_q1_srt, ver_q2_srt, sctrPos_srt

    logical, dimension(4)                               :: drctFlag, drctFlag_srt

    real(dp), dimension(:,:), allocatable               :: Mvrtq1, Mvrtq2

    real(dp), dimension(4)                              :: Mq1q2s

    integer                                             :: posInIrr_q1, existsInIrr_q1, &
                                                         & posInIrr_q2, existsInIrr_q2

    real(dp), dimension(3)                              :: Flmb1, Flmb2, Fl2_plus_l1

    !================================= Local variable ===================================!

    mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)

    NumTetra = size( Qpoints%tetrahedrons, 2 )

    allocate( Mq(Ndof) )
    allocate( Mvrtq1(Ndof, 4) )
    allocate( Mvrtq2(Ndof, 4) )

    irr_q_loop: do qi = 1, my_Qsize

        qindx = my_offset + qi

        q_int = Qpoints%irr_q_int(1:3, qindx)
        qindx2 = Qpoints%indx_map( dot_product( modulo( q_int, mesh ), mesh_mul ) + 1 )

        Mq = phonon_dat%omega(:, qindx2)

        q0chk(1) = all( q_int == 0 )

        !-!s_loop: do s = 1, Ndof

        !-!    Mqs = Mq(s)

        loop_Tetra: do NT = 1, NumTetra

            ver_q1 = Qpoints%tetrahedrons(:, NT)

            vert_loop1: do ii = 1, 4

                ver_q1(ii) = Qpoints%indx_map( ver_q1(ii) ) ! ** !
                Mvrtq1(:, ii) = phonon_dat%omega(:, ver_q1(ii))

                sctrPos(ii) = this%ConsrvMinus(qi)%ConsrvQ( ver_q1(ii) )%SctrElPos

                drctFlag(ii) = this%ConsrvMinus(qi)%ConsrvQ( ver_q1(ii) )%flag

                ver_q2(ii) = this%ConsrvMinus(qi)%ConsrvQ( ver_q1(ii) )%q2indx
                Mvrtq2(:, ii) = phonon_dat%omega(:, ver_q2(ii))

            end do vert_loop1

            s_loop: do s = 1, Ndof
                Mqs = Mq(s)

                s1_loop: do s1 = 1, Ndof
                    s2_loop: do s2 = 1, Ndof

                        vert_loop2: do kk = 1, 4

                            Mq1q2s(kk) = ( Mvrtq1(s1, kk) + Mvrtq2(s2, kk) )

                        end do vert_loop2

                        ver_q1_srt = ver_q1
                        ver_q2_srt = ver_q2
                        sctrPos_srt = sctrPos
                        drctFlag_srt = drctFlag

                        call my_sort(Mq1q2s, ver_q1_srt, ver_q2_srt, sctrPos_srt, drctFlag_srt)

                        enterCase = .false.
                        tetra_val = 0.0_dp

                        tetra_coef_case: if ( (Mqs > Mq1q2s(1)) .and. (Mqs < Mq1q2s(2)) ) then

                            ! ------------------- value of r and values of I at vertices ------------------- !
                            g = ( 3.0_dp * ( (Mqs - Mq1q2s(1)) ** 2 ) ) / &
                              & ( (Mq1q2s(2) - Mq1q2s(1)) * & 
                              &   (Mq1q2s(3) - Mq1q2s(1)) * &
                              &   (Mq1q2s(4) - Mq1q2s(1)) )

                            Iver(1) = ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(1) - Mq1q2s(2)) ) + &
                                      & ( (Mqs - Mq1q2s(3)) / (Mq1q2s(1) - Mq1q2s(3)) ) + &
                                      & ( (Mqs - Mq1q2s(4)) / (Mq1q2s(1) - Mq1q2s(4)) ) ) / 3.0_dp

                            do vv = 2, 4
                                Iver(vv) = ( (Mqs - Mq1q2s(1)) / (Mq1q2s(vv) - Mq1q2s(1)) ) / 3.0_dp
                            end do
                            ! ------------------- value of r and values of I at vertices ------------------- !

                            enterCase = .true.


                        else if ( (Mqs > Mq1q2s(2)) .and. (Mqs < Mq1q2s(3)) ) then tetra_coef_case

                            ! ------------------- value of r and values of I at vertices ------------------- !
                            g = ( 3.0_dp / ( (Mq1q2s(2) - Mq1q2s(3)) * &
                              &              (Mq1q2s(4) - Mq1q2s(1)) ) ) * &
                              & ( ( ( (Mqs - Mq1q2s(1)) * (Mqs - Mq1q2s(3)) ) / &
                              &       (Mq1q2s(3) - Mq1q2s(1)) ) - &
                              &   ( ( (Mqs - Mq1q2s(2)) * (Mqs - Mq1q2s(4)) ) / &
                              &       (Mq1q2s(2) - Mq1q2s(4)) ) )

                            Iver(1) = ( ( (Mqs - Mq1q2s(4)) / (Mq1q2s(1) - Mq1q2s(4)) ) / 3.0_dp ) + &
                                    & ( ( ( (Mqs - Mq1q2s(3)) / (Mq1q2s(1) - Mq1q2s(3)) ) * &
                                    &   ( (Mqs - Mq1q2s(1)) / (Mq1q2s(3) - Mq1q2s(1)) ) * &
                                    &   ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) ) / &
                                    &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )

                            Iver(2) = ( ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) / 3.0_dp ) + &
                                    & ( ( ( ((Mqs - Mq1q2s(4)) / (Mq1q2s(2) - Mq1q2s(4))) ** 2 ) * &
                                    &   ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) ) / &
                                    &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )

                            Iver(3) = ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) / 3.0_dp ) + &
                                    & ( ( ( ((Mqs - Mq1q2s(1)) / (Mq1q2s(3) - Mq1q2s(1))) ** 2 ) * &
                                    &   ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) ) / &
                                    &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )
                        
                            Iver(4) = ( ( (Mqs - Mq1q2s(1)) / (Mq1q2s(4) - Mq1q2s(1)) ) / 3.0_dp ) + &
                                    & ( ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(4) - Mq1q2s(2)) ) * &
                                    &   ( (Mqs - Mq1q2s(4)) / (Mq1q2s(2) - Mq1q2s(4)) ) * &
                                    &   ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) ) / &
                                    &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )
                            ! ------------------- value of r and values of I at vertices ------------------- !

                            enterCase = .true.


                        else if ( (Mqs > Mq1q2s(3)) .and. (Mqs < Mq1q2s(4)) ) then tetra_coef_case

                            ! ------------------- value of r and values of I at vertices ------------------- !
                            g = ( 3.0_dp * ((Mqs - Mq1q2s(4)) ** 2) ) / &
                              & ( ( Mq1q2s(4) - Mq1q2s(1) ) * ( Mq1q2s(4) - Mq1q2s(2) ) * &
                              &   ( Mq1q2s(4) - Mq1q2s(3) ) )

                            do vv = 1, 3
                                Iver(vv) = ( (Mqs - Mq1q2s(4)) / (Mq1q2s(vv) - Mq1q2s(4)) ) / 3.0_dp
                            end do

                            Iver(4) = ( ( (Mqs - Mq1q2s(1)) / (Mq1q2s(4) - Mq1q2s(1)) ) + &
                                    &   ( (Mqs - Mq1q2s(2)) / (Mq1q2s(4) - Mq1q2s(2)) ) + &
                                    &   ( (Mqs - Mq1q2s(3)) / (Mq1q2s(4) - Mq1q2s(3)) ) ) / 3.0_dp
                            ! ------------------- value of r and values of I at vertices ------------------- !

                            enterCase = .true.

                        else tetra_coef_case

                            does_not_enters: if ( (Mqs >= Mq1q2s(1)) .and. &
                                                & (Mqs <= Mq1q2s(4)) .and. (.not. enterCase) ) then
                                
                                write(*, 55) Mqs
                                write(*, 65) Mq1q2s

                                55 FORMAT("WARNING: Mqs = ", F12.5, "THz. It does not enters any cases.")
                                65 FORMAT("         The values of Mq1q2s at vertices: [", 3(F12.5, ', '), F12.5, "] THz")

                            end if does_not_enters

                            enterCase = .false.

                        end if tetra_coef_case


                        chk_enterCase: if ( enterCase ) then

                            chk_sumeq1: if ( dabs(sum(Iver) - 1.0_dp) < zero_prec ) then

                                ! ----------------------------- sum over the 4 vertices ------------------------------ !
                                vert_loop3: do vv = 1, 4

                                    !  What is l1 and l2 at vertices vv ( let's say vv = 1), and find F(l2) - F(l1)  !

                                    posInIrr_q1 = Qpoints%q_pnt_int( 4, ver_q1_srt(vv) )
                                    existsInIrr_q1 = Qpoints%q_pnt_int( 5, ver_q1_srt(vv) )

                                    DirectOrIndirect_q1: if ( existsInIrr_q1 == 1 ) then
                                        Flmb1 = F_old(:, s1, posInIrr_q1)

                                    else DirectOrIndirect_q1
                                        Flmb1 = -1.0_dp * F_old(:, s1, posInIrr_q1)

                                    end if DirectOrIndirect_q1

                                    posInIrr_q2 = Qpoints%q_pnt_int( 4, ver_q2_srt(vv) )
                                    existsInIrr_q2 = Qpoints%q_pnt_int( 5, ver_q2_srt(vv) )

                                    DirectOrIndirect_q2: if ( existsInIrr_q2 == 1 ) then
                                        Flmb2 = F_old(:, s2, posInIrr_q2)

                                    else DirectOrIndirect_q2
                                        Flmb2 = -1.0_dp * F_old(:, s2, posInIrr_q2)

                                    end if DirectOrIndirect_q2

                                    Fl2_plus_l1 = Flmb2 + Flmb1

                                    per_symm: if ( drctFlag_srt(vv) ) then

                                        tetra_val = tetra_val + ( (Iver(vv) * &
                                                  & this%WPlus(qi)%q1q2( sctrPos_srt(vv) )%ss1s2(s2, s1, s)) * &
                                                  & Fl2_plus_l1 )

                                    else per_symm

                                        tetra_val = tetra_val + ( (Iver(vv) * &
                                                  & this%WMinus(qi)%q1q2( sctrPos_srt(vv) )%ss1s2(s1, s2, s)) * &
                                                  & Fl2_plus_l1 )

                                    end if per_symm

                                end do vert_loop3
                                ! ----------------------------- sum over the 4 vertices ------------------------------ !
                                tetra_val = 0.5_dp * OneSixth * g * tetra_val

                                dFnew(:, s, qindx) = dFnew(:, s, qindx) + tetra_val ! ** !

                            else chk_sumeq1
                                
                                tetra_val = 0.0_dp
                                write(*, 45) Iver

                            end if chk_sumeq1

                        else chk_enterCase

                            tetra_val = 0.0_dp

                        end if chk_enterCase

                    end do s2_loop
                end do s1_loop
            end do s_loop

        end do loop_Tetra

    end do irr_q_loop

    45 FORMAT("WARNING: sum( Iver(:) ) .ne. 1. Values of Iver(1:4): [", 3(F12.5, ', '), F12.5, "]")

    deallocate( Mq )
    deallocate( Mvrtq1, Mvrtq2 )

end subroutine TetrahedronIterMinus

