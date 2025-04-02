
module MomentumConsrv3rd_mod

    use TetrahedronQ,            only : q_tetrahedron

    implicit none
    private

    public              :: MomentumConsrv3rd

contains

    subroutine MomentumConsrv3rd( qTetra, Nq1, Nq1Perm, i0, mesh, Allq1q2, q1q2_permuted )

        implicit none

        type(q_tetrahedron), intent(in)                             :: qTetra
        integer, intent(in)                                         :: Nq1, Nq1Perm, i0
        integer, dimension(3), intent(in)                           :: mesh

        integer, dimension(4, Nq1), intent(out)                     :: Allq1q2
        integer, dimension(2, Nq1Perm), intent(out)                 :: q1q2_permuted

        !============================== Local variables =============================!

        integer, parameter                                          :: order = 3, num_per = 1 != factorial(order - 1) - 1
        integer, dimension(order, num_per)                          :: permutation

        integer                                                     :: i1, iG, q1_grid_indx, &
                                                                     & num_cell, q2_cnt, pp, &
                                                                     & c1, c2, c3, q2_grid_indx, &
                                                                     & q1q2_count, G_cnt, total_count
        
        integer, dimension(3)                                       :: bound, q0_int, G_vec, &
                                                                     & q1_int, q2_int, mesh_mul

        integer, allocatable, dimension(:)                          :: sgn_mul
        integer, dimension(:, :), allocatable                       :: recp_latt

        logical                                                     :: WithinBZ
        logical, dimension(Nq1)                                     :: AlreadyCovered

        !============================== Local variables =============================!

        AlreadyCovered(:) = .false.

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

        total_count = 0
        q1q2_count = 0
        q1_loop: do i1 = 1, Nq1


            q1_int = qTetra%q_pnt_int(1:3, i1)

            q1_grid_indx = dot_product( modulo( q1_int, mesh ), mesh_mul ) + 1

            ChkIfNotAlradyCovered: if ( .not. (AlreadyCovered(q1_grid_indx)) ) then

                q2_cnt = 0
                G_loop: do iG = 1, G_cnt

                    G_vec = recp_latt(:, iG)

                    q2_int = G_vec - (q0_int + q1_int)
                    q2_grid_indx = dot_product( modulo( q2_int, mesh ), mesh_mul ) + 1

                    WithinBZ = CheckIfWithinBZ( bound, q2_int, q2_grid_indx, Nq1 )

                    if ( WithinBZ ) then
                        q2_cnt = q2_cnt + 1
                        EXIT G_loop
                    end if

                end do G_loop

                chk_q2_found: if ( q2_cnt /= 0 ) then

                    q1q2_count = q1q2_count + 1
                    Allq1q2(1:4, q1_grid_indx) = (/1, q1q2_count, q1_grid_indx, q2_grid_indx/)

                    AlreadyCovered( q1_grid_indx ) = .true.
                    total_count = total_count + 1

                    q1q2_permuted(1:2, q1q2_count) = (/q1_grid_indx, q2_grid_indx/)

                    permutation = Permutation3rd( q1_int, q2_int, q1_grid_indx, q2_grid_indx) !, G_vec, q0_int )

                    AllPossiblePermutation: do pp = 1, num_per

                        if ( permutation(1, pp) == 1 ) then

                            Allq1q2(1:4, permutation(2, pp)) = (/0, q1q2_count, permutation(2, pp), permutation(3, pp)/)

                            AlreadyCovered( permutation(2, pp) ) = .true.
                            total_count = total_count + 1

                        end if

                    end do AllPossiblePermutation

                else chk_q2_found

                    write(*, 66) i0, q0_int, q1_int
                    66 FORMAT( "WARNING: No q2 found conserving quasi-momentum. q0 (i0 = ", I5, &
                             & ") = (", 2(I4, ', '),I4, ") and q1 = (", &
                             & 2(I4, ', '), I4, ")." )

                end if chk_q2_found

            end if ChkIfNotAlradyCovered

        end do q1_loop

        !-! write(*, *) total_count, q1q2_count
        ! ** debug ** !
        if ( (q1q2_count /= Nq1Perm) .or. (total_count /= Nq1) ) then
            write(*, 77) q1q2_count, Nq1Perm, total_count, Nq1
            77 FORMAT( "ERROR: q1q2_count /= Nq1Perm ", I6, '  ', I6, &
                     & ", or total_count /= Nq1 ", I6, '  ', I6 )
            ERROR STOP
        end if
        ! ** debug ** !

        deallocate( sgn_mul, recp_latt )

    end subroutine MomentumConsrv3rd


    Function CheckIfWithinBZ( bound, q2_int, q2_grid_indx, Nq1 ) Result( WithinBZ )

        implicit none

        integer, dimension(3), intent(in)               :: bound, q2_int
        integer, intent(in)                             :: q2_grid_indx, Nq1

        logical                                         :: WithinBZ !Result

        WithinBZ = ( ( .not. (any(q2_int > bound) .or. any(q2_int < -bound)) ) .and. &
                   & ( (q2_grid_indx <= Nq1) .and. (q2_grid_indx >= 1) ) )

    end Function CheckIfWithinBZ


    Function Permutation3rd( q1_int, q2_int, q1_grid_indx, q2_grid_indx ) Result(permutation)

        implicit none

        integer, dimension(3), intent(in)                   :: q1_int, q2_int!, G_vec, q0_int
        integer, intent(in)                                 :: q1_grid_indx, q2_grid_indx

        integer, dimension(3, 1)                            :: permutation !Result
        
        ! ============================== Local Variables ============================== !

        integer, parameter                                  :: order = 3, &
                                                             & num_per = 1 != factorial(order - 1) - 1

        integer, dimension(order-1, num_per)                :: per
        integer, dimension(3, order-1)                      :: initial_q1q2, permuted_q1q2
        integer, dimension(order-1)                         :: initial_grid_indx, &
                                                             & permuted_grid_indx

        integer                                             :: pp

        ! ============================== Local Variables ============================== !

        permutation(:, :) = 0

        per(:, 1) = (/2, 1/)

        initial_q1q2(:, 1) = q1_int
        initial_q1q2(:, 2) = q2_int

        initial_grid_indx(:) = (/q1_grid_indx, q2_grid_indx/)

        do pp = 1, num_per
            
            permuted_q1q2(:, 1) = initial_q1q2( :, per(1, pp) )
            permuted_q1q2(:, 2) = initial_q1q2( :, per(2, pp) )

            permuted_grid_indx(1) = initial_grid_indx( per(1, pp) )
            permuted_grid_indx(2) = initial_grid_indx( per(2, pp) )

            permutation(2:3, pp) = permuted_grid_indx(:)

            if ( any( (permuted_q1q2 - initial_q1q2) /= 0 ) ) then
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

    end Function Permutation3rd

end module MomentumConsrv3rd_mod

