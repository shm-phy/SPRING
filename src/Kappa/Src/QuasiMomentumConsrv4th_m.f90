
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

module MomentumConsrv4th_mod

    use TetrahedronQ,            only : q_tetrahedron

    implicit none
    private

    public              :: PermutationReducedQs, MomentumConsrv4th

contains

    subroutine MomentumConsrv4th( qTetra, Nq1, Nq1q2Perm, i0, mesh, Allq1q2q3, q1q2q3_permuted )

        implicit none

        type(q_tetrahedron), intent(in)                             :: qTetra
        integer, intent(in)                                         :: Nq1, Nq1q2Perm, i0
        integer, dimension(3), intent(in)                           :: mesh

        integer, dimension(3, Nq1, Nq1), intent(out)                :: Allq1q2q3
        integer, dimension(3, Nq1q2Perm), intent(out)               :: q1q2q3_permuted

        !============================== Local variables =============================!

        character(len=128)                                          :: msg
        integer, parameter                                          :: order = 4, num_per = 5 != factorial(order - 1) - 1
        integer, dimension(order, num_per)                          :: permutation

        integer                                                     :: i1, iG, q1_grid_indx, &
                                                                     & num_cell, q3_cnt, pp, &
                                                                     & c1, c2, c3, q2_grid_indx, &
                                                                     & q1q2q3_count, G_cnt, &
                                                                     & i2, q3_grid_indx, total_count
        
        integer, dimension(3)                                       :: bound, q0_int, G_vec, &
                                                                     & q1_int, q2_int, mesh_mul, &
                                                                     & q3_int

        integer, allocatable, dimension(:)                          :: sgn_mul
        integer, dimension(:, :), allocatable                       :: recp_latt

        integer                                                     :: istat

        logical                                                     :: WithinBZ
        !logical, dimension(Nq1, Nq1)                                :: AlreadyCovered
        logical, dimension(:, :), allocatable                       :: AlreadyCovered

        !============================== Local variables =============================!

        allocate( AlreadyCovered(Nq1, Nq1), STAT=istat, ERRMSG=msg )
        if ( istat /= 0 ) then
            write(*, 28) this_image(), msg
            28 FORMAT( " Memory allocation ERROR in image : ", I5, ". Error message: ", A128 )
            ERROR STOP
        end if

        q1q2q3_permuted(:, :) = 0

        AlreadyCovered(:, :) = .false.

        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)
        bound = mesh / 2

        !=========== Find all the possible reciprocal lattice vectors (G) ===========!

        allocate( sgn_mul(3) )
        sgn_mul = (/0, 1, -1/)

        !allocate( sgn_mul(5) )
        !sgn_mul = (/0, 1, -1, 2, -2/)

        num_cell = size( sgn_mul )

        allocate( recp_latt(3, num_cell**3) )

        G_cnt = 0
        c3_loop: do c3 = 1, num_cell
            c2_loop: do c2 = 1, num_cell
                c1_loop: do c1 = 1, num_cell

                    G_cnt = G_cnt + 1
                    recp_latt(:, G_cnt) = (/sgn_mul(c1)*mesh(1), sgn_mul(c2)*mesh(2), sgn_mul(c3)*mesh(3)/)

                end do c1_loop
            end do c2_loop
        end do c3_loop

        G_cnt = num_cell**3

        !=========== Find all the possible reciprocal lattice vectors (G) ===========!

        q0_int(:) = qTetra%irr_q_pnt_int(:, i0)

        q1q2q3_count = 0
        total_count = 0
        q1_loop: do i1 = 1, Nq1

            q1_int = qTetra%q_pnt_int(1:3, i1)

            q1_grid_indx = dot_product( modulo( q1_int, mesh ), mesh_mul ) + 1

            q2_loop: do i2 = 1, Nq1

                q2_int = qTetra%q_pnt_int(1:3, i2)
                q2_grid_indx = dot_product( modulo( q2_int, mesh ), mesh_mul ) + 1

                ChkIfNotAlradyCovered: if ( .not. (AlreadyCovered(q2_grid_indx, q1_grid_indx)) ) then

                    q3_cnt = 0
                    G_loop: do iG = 1, G_cnt

                        G_vec = recp_latt(:, iG)

                        q3_int = G_vec - (q0_int + q1_int + q2_int)
                        q3_grid_indx = dot_product( modulo( q3_int, mesh ), mesh_mul ) + 1

                        WithinBZ = CheckIfWithinBZ( bound, q3_int, q3_grid_indx, Nq1 )

                        if ( WithinBZ ) then
                            q3_cnt = q3_cnt + 1
                            EXIT G_loop
                        end if

                    end do G_loop

                    chk_q2_found: if ( q3_cnt /= 0 ) then

                        q1q2q3_count = q1q2q3_count + 1
                        !Allq1q2q3(1:5, q2_grid_indx, q1_grid_indx) = (/1, q1q2q3_count, q1_grid_indx, q2_grid_indx, q3_grid_indx/)

                        Allq1q2q3(1:3, q2_grid_indx, q1_grid_indx) = (/0, q1q2q3_count, q3_grid_indx/)
                        AlreadyCovered( q2_grid_indx, q1_grid_indx ) = .true.
                        total_count = total_count + 1

                        q1q2q3_permuted(1:3, q1q2q3_count) = (/q1_grid_indx, q2_grid_indx, q3_grid_indx/)

                        permutation = Permutation4th( q1_int, q2_int, q3_int, &
                                                    & q1_grid_indx, q2_grid_indx, q3_grid_indx ) !, G_vec, q0_int )

                        AllPossiblePermutation: do pp = 1, num_per

                            if ( permutation(1, pp) == 1 ) then

                                Allq1q2q3(1:3, permutation(3, pp), permutation(2, pp)) = &
                                 & (/pp, q1q2q3_count, permutation(4, pp)/)
                                AlreadyCovered( permutation(3, pp), permutation(2, pp) ) = .true.
                                total_count = total_count + 1

                            end if

                        end do AllPossiblePermutation

                    else chk_q2_found

                        write(*, 66) i0, q0_int, q1_int, q2_int
                        66 FORMAT( "WARNING: No q3 found conserving quasi-momentum. q0 (i0 = ", I5, &
                                 & ") = (", 2(I4, ', '),I4, ") , q1 = (", &
                                 & 2(I4, ', '), I4, ") and q2 = (", 2(I4, ', '), I4, ")" )

                    end if chk_q2_found

                end if ChkIfNotAlradyCovered

            end do q2_loop

        end do q1_loop

        !-! write(*, *) total_count, q1q2q3_count
        ! ** debug ** !
        if ( (total_count /= (Nq1**2)) .or. (q1q2q3_count /= Nq1q2Perm) ) then

            write(*, 77) total_count, (Nq1**2)
            77 FORMAT( "ERROR:  total_count /= (Nq1**2) ", I9, '  ', I9, &
                     & ", or q1q2q3_count /= Nq1q2Perm ", I9, '   ', I9 )
            write(*, 88)
            88 FORMAT( "If it is second case use prime number as N in NxNxN mesh" )

            ERROR STOP

        end if
        ! ** debug ** !

        deallocate( sgn_mul, recp_latt )
        deallocate( AlreadyCovered )

    end subroutine MomentumConsrv4th


    Function CheckIfWithinBZ( bound, q3_int, q3_grid_indx, Nq1 ) Result( WithinBZ )

        implicit none

        integer, dimension(3), intent(in)               :: bound, q3_int
        integer, intent(in)                             :: q3_grid_indx, Nq1

        logical                                         :: WithinBZ !Result

        WithinBZ = ( ( .not. (any(q3_int > bound) .or. any(q3_int < -bound)) ) .and. &
                   & ( (q3_grid_indx <= Nq1) .and. (q3_grid_indx >= 1) ) )

    end Function CheckIfWithinBZ


    Function Permutation4th( q1_int, q2_int, q3_int, &
                           & q1_grid_indx, q2_grid_indx, q3_grid_indx ) Result(permutation)

        implicit none

        integer, dimension(3), intent(in)                   :: q1_int, q2_int, q3_int!, G_vec, q0_int
        integer, intent(in)                                 :: q1_grid_indx, q2_grid_indx, &
                                                             & q3_grid_indx

        integer, dimension(4, 5)                            :: permutation !Result
        
        ! ============================== Local Variables ============================== !

        integer, parameter                                  :: order = 4, &
                                                             & num_per = 5 != factorial(order - 1) - 1

        integer, dimension(order-1, num_per)                :: per
        integer, dimension(3, order-1)                      :: initial_q1q2q3, permuted_q1q2q3
        integer, dimension(order-1)                         :: initial_grid_indx, &
                                                             & permuted_grid_indx

        integer                                             :: pp, qq

        logical                                             :: NoNewPermutation

        ! ============================== Local Variables ============================== !

        permutation(:, :) = 0

        per(:, 1) = (/1, 3, 2/)
        per(:, 2) = (/2, 1, 3/)
        per(:, 3) = (/2, 3, 1/)
        per(:, 4) = (/3, 1, 2/)
        per(:, 5) = (/3, 2, 1/)

        initial_q1q2q3(:, 1) = q1_int
        initial_q1q2q3(:, 2) = q2_int
        initial_q1q2q3(:, 3) = q3_int

        initial_grid_indx(:) = (/q1_grid_indx, q2_grid_indx, q3_grid_indx/)

        do pp = 1, num_per
            
            permuted_q1q2q3(:, 1) = initial_q1q2q3( :, per(1, pp) )
            permuted_q1q2q3(:, 2) = initial_q1q2q3( :, per(2, pp) )
            permuted_q1q2q3(:, 3) = initial_q1q2q3( :, per(3, pp) )

            permuted_grid_indx(1) = initial_grid_indx( per(1, pp) )
            permuted_grid_indx(2) = initial_grid_indx( per(2, pp) )
            permuted_grid_indx(3) = initial_grid_indx( per(3, pp) )

            permutation(2:4, pp) = permuted_grid_indx(:)

            NoNewPermutation = all( (permuted_q1q2q3 - initial_q1q2q3) == 0 ) !.or. &
                             !& all( (permuted_grid_indx(:) - initial_grid_indx(:)) == 0 )
            !NoNewPermutation = all( (permuted_grid_indx(:) - initial_grid_indx(:)) == 0 ) 

            inner_loop: do qq = 1, (pp - 1)

                NoNewPermutation = (NoNewPermutation) .or. &
                                 & all( (permuted_grid_indx(:) - permutation(2:4, qq)) == 0 )

            end do inner_loop

            if ( .not. NoNewPermutation ) then

                permutation(1, pp) = 1

            end if

            !** debug **!
            !-! write(*, 5)  q0_int, initial_q1q2(:, 1), initial_q1q2(:, 2), G_vec
            !-! 5 FORMAT( "( ", 2(I4, ' '), I4, ") + (", 2(I4, ' '), I4, ") + (", &
            !-!         & "( ", 2(I4, ' '), I4, " ) = (", 2(I4, ' '), I4, " )" )
            !-! write(*, 10)  permuted_q1q2(:, 1), permuted_q1q2(:, 2), initial_q1q2(:, 1), initial_q1q2(:, 2), &
            !-!             & permutation(1, pp), any( (permuted_q1q2 - initial_q1q2) /= 0 )

            !-! 10 FORMAT( "( ", 2(I4, ' '), I4, " | ", 2(I4, ' '), I4, ")", " - ", &
            !-!          & "( ", 2(I4, ' '), I4, " | ", 2(I4, ' '), I4, ")", "  ", I3, "   ", L)
            !-! write(*, 20) (permuted_q1q2 - initial_q1q2)
            !-! 20 FORMAT(6(I2, ' '))
            !-! write(*, 30) ((permuted_q1q2 - initial_q1q2) /= 0)
            !-! 30 FORMAT(6(L, ' '))
            !-! write(*, *) any( (permuted_q1q2 - initial_q1q2) /= 0 )
            !** debug **!

        end do

    end Function Permutation4th

    ! ** Use this when N (in NxNxN mesh ) is not a prime number ** !
    Function PermutationReducedQs( qTetra, Nq1, i0, mesh) Result( Nq1q2q3_per )

        implicit none

        type(q_tetrahedron), intent(in)                             :: qTetra
        integer, intent(in)                                         :: Nq1, i0
        integer, dimension(3), intent(in)                           :: mesh

        integer                                                     :: Nq1q2q3_per ! Result

        !============================== Local variables =============================!

        character(len=128)                                          :: msg

        integer, parameter                                          :: order = 4, num_per = 5 != factorial(order - 1) - 1
        integer, dimension(order, num_per)                          :: permutation

        integer                                                     :: i1, iG, q1_grid_indx, &
                                                                     & num_cell, q3_cnt, pp, &
                                                                     & c1, c2, c3, q2_grid_indx, &
                                                                     & G_cnt, i2, q3_grid_indx, &
                                                                     & total_count
        integer                                                     :: istat

        integer, dimension(3)                                       :: bound, q0_int, G_vec, &
                                                                     & q1_int, q2_int, mesh_mul, &
                                                                     & q3_int

        integer, allocatable, dimension(:)                          :: sgn_mul
        integer, dimension(:, :), allocatable                       :: recp_latt

        logical                                                     :: WithinBZ
        logical, dimension(:, :), allocatable                       :: AlreadyCovered

        !============================== Local variables =============================!

        allocate( AlreadyCovered(Nq1, Nq1), STAT=istat, ERRMSG=msg )
        if ( istat /= 0 ) then
            write(*, 28) this_image(), msg
            28 FORMAT( " Memory allocation ERROR in image : ", I5, ". Error message: ", A128 )
            ERROR STOP
        end if

        AlreadyCovered(:, :) = .false.

        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)
        bound = mesh / 2

        !=========== Find all the possible reciprocal lattice vectors (G) ===========!

        allocate( sgn_mul(3) )
        sgn_mul = (/0, 1, -1/)

        !allocate( sgn_mul(5) )
        !sgn_mul = (/0, 1, -1, 2, -2/)

        num_cell = size( sgn_mul )

        allocate( recp_latt(3, num_cell**3) )

        G_cnt = 0
        c3_loop: do c3 = 1, num_cell
            c2_loop: do c2 = 1, num_cell
                c1_loop: do c1 = 1, num_cell

                    G_cnt = G_cnt + 1
                    recp_latt(:, G_cnt) = (/sgn_mul(c1)*mesh(1), sgn_mul(c2)*mesh(2), sgn_mul(c3)*mesh(3)/)

                end do c1_loop
            end do c2_loop
        end do c3_loop

        G_cnt = num_cell**3

        !=========== Find all the possible reciprocal lattice vectors (G) ===========!

        q0_int(:) = qTetra%irr_q_pnt_int(:, i0)

        Nq1q2q3_per = 0
        total_count = 0
        q1_loop: do i1 = 1, Nq1

            q1_int = qTetra%q_pnt_int(1:3, i1)

            q1_grid_indx = dot_product( modulo( q1_int, mesh ), mesh_mul ) + 1

            q2_loop: do i2 = 1, Nq1

                q2_int = qTetra%q_pnt_int(1:3, i2)
                q2_grid_indx = dot_product( modulo( q2_int, mesh ), mesh_mul ) + 1

                ChkIfNotAlradyCovered: if ( .not. (AlreadyCovered(q2_grid_indx, q1_grid_indx)) ) then

                    q3_cnt = 0
                    G_loop: do iG = 1, G_cnt

                        G_vec = recp_latt(:, iG)

                        q3_int = G_vec - (q0_int + q1_int + q2_int)
                        q3_grid_indx = dot_product( modulo( q3_int, mesh ), mesh_mul ) + 1

                        WithinBZ = CheckIfWithinBZ( bound, q3_int, q3_grid_indx, Nq1 )

                        if ( WithinBZ ) then
                            q3_cnt = q3_cnt + 1
                            EXIT G_loop
                        end if

                    end do G_loop

                    chk_q2_found: if ( q3_cnt /= 0 ) then

                        Nq1q2q3_per = Nq1q2q3_per + 1

                        AlreadyCovered( q2_grid_indx, q1_grid_indx ) = .true.
                        total_count = total_count + 1

                        permutation = Permutation4th( q1_int, q2_int, q3_int, &
                                                    & q1_grid_indx, q2_grid_indx, q3_grid_indx ) !, G_vec, q0_int )

                        AllPossiblePermutation: do pp = 1, num_per

                            if ( permutation(1, pp) == 1 ) then

                                AlreadyCovered( permutation(3, pp), permutation(2, pp) ) = .true.
                                total_count = total_count + 1

                            end if

                        end do AllPossiblePermutation

                    else chk_q2_found

                        write(*, 66) i0, q0_int, q1_int, q2_int
                        66 FORMAT( "WARNING: No q3 found conserving quasi-momentum. q0 (i0 = ", I5, &
                                 & ") = (", 2(I4, ', '),I4, ") , q1 = (", &
                                 & 2(I4, ', '), I4, ") and q2 = (", 2(I4, ', '), I4, ")" )

                    end if chk_q2_found

                end if ChkIfNotAlradyCovered

            end do q2_loop

        end do q1_loop

        ! ** debug ** !
        if ( total_count /= (Nq1**2) ) then
            write(*, 77) total_count, (Nq1**2)
            77 FORMAT( "ERROR:  total_count /= (Nq1**2) ", I9, '  ', I9)
            ERROR STOP
        end if
        ! ** debug ** !

        deallocate( sgn_mul, recp_latt )
        deallocate( AlreadyCovered )

    end Function PermutationReducedQs
    ! ** Use this when N (in NxNxN mesh ) is not a prime number ** !

end module MomentumConsrv4th_mod

