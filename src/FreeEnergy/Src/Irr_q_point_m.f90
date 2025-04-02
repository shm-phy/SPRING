
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


module Irr_q_point

    use kinds,          only : dp
    use unit_cell,      only : cell
    use IndxChng,       only : fold_indx

    implicit none
    private

    ! The datatype q_points_data contains all the q-points information. !
    type, public        :: q_points_data

        integer, allocatable, dimension(:)              :: indx_map
        integer, allocatable, dimension(:,:)            :: q_pnt_int
        integer, allocatable, dimension(:,:)            :: irr_q_int
        integer, allocatable, dimension(:,:)            :: tetrahedrons
        real(dp), allocatable, dimension(:,:)           :: q_pnt
        real(dp), allocatable, dimension(:,:)           :: irr_q

        contains
            procedure, public, pass                     :: make_tetrahedron, &
                                                         & clean_Qdata
            !FINAL                                       :: clean_Qdata

    end type q_points_data

contains

    subroutine make_tetrahedron(this, sys, mesh, shift, time_rev)

        implicit none

        class(q_points_data)                                        :: this

        type(cell), intent(in)                                      :: sys
        integer, dimension(3), intent(in)                           :: mesh, shift
        integer, intent(in)                                         :: time_rev

        !===================================== Local variable =========================================!

        integer, dimension( mesh(1)*mesh(2)*mesh(3) )               :: map

        integer, dimension(:), allocatable                          :: tmp_arr
        integer                                                     :: num_irr_q, num_grid_pnt, irr_q_count, &
                                                                     & q_indx_unq, qminus_indx_unq, &
                                                                     & indx_qfirst, indx_minusq, my_num
                                                                     
        integer, dimension(3)                                       :: mesh_mul, grid_pnt
        real(dp), dimension(3)                                      :: qvec

        !! ******************************** Tetrahedron ******************************** !!

        integer, dimension(3, 8)                                    :: box
        integer, dimension(4, 8)                                    :: box_trans
        integer, dimension(3)                                       :: indx, vec_trans, &
                                                                     & vert, vert_fold

        logical, dimension( mesh(1)*mesh(2)*mesh(3) )               :: uniq_count
        integer, dimension( mesh(1)*mesh(2)*mesh(3) )               :: grid_indx_tetra
        integer, dimension( 3, mesh(1)*mesh(2)*mesh(3) )            :: grid_pnt_tetra
        integer, dimension( 4, mesh(1)*mesh(2)*mesh(3)*6 )          :: tetrahedrons

        integer                                                     :: grd_indx_unq, my_count, &
                                                                     & my_tetra_count
        integer                                                     :: ix, iy, iz, v

        !! ******************************** Tetrahedron ******************************** !!
                                                                     
        integer                                                     :: ii

        !*! Debug !*!
        !*! real(dp), dimension(3)                                      :: q1, q2, q3, q4, c
        !*! real(dp)                                                    :: vol
        !*! integer                                                     :: kk
        !*! Debug !*!

        !===================================== Local variable =========================================!

        num_grid_pnt = product( mesh )
        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)

        !! =================================== devide the whole BZ into tetrahedrons ========================================= !!
        box(:, 1) = (/0, 0, 0/) !1
        box(:, 2) = (/1, 0, 0/) !5
        box(:, 3) = (/0, 0, 1/) !2
        box(:, 4) = (/0, 1, 1/) !4
        box(:, 5) = (/1, 0, 1/) !6
        box(:, 6) = (/1, 1, 1/) !8
        box(:, 7) = (/1, 1, 0/) !7
        box(:, 8) = (/0, 1, 0/) !3

        uniq_count(:) = .true.
        my_count = 1
        my_tetra_count = 1

        iz_loop: do iz = 0, (mesh(3) - 1)
            iy_loop: do iy = 0, (mesh(2) - 1)
                ix_loop: do ix = 0, (mesh(1) - 1)

                    indx = (/ix, iy, iz/)
                    call fold_indx(mesh, indx, vec_trans)

                    v_loop: do v = 1, 8

                        vert = box(:, v) + vec_trans
                        call fold_indx(mesh, vert, vert_fold)

                        grd_indx_unq = dot_product( modulo( vert_fold, mesh ), mesh_mul ) + 1

                        !*! degub !*!
                        !*! debg_chk: if ( grd_indx_unq > num_grid_pnt ) then

                        !*!     write(*, 40)
                        !*!     40 FORMAT("ERROR: grd_indx_unq out of bound")

                        !*! end if debg_chk
                        !*! degub !*!

                        box_trans(1:3, v) = vert_fold
                        box_trans(4, v) = grd_indx_unq

                        unq_test: if ( uniq_count(grd_indx_unq) ) then

                            grid_indx_tetra(grd_indx_unq) = my_count ! $ !
                            grid_pnt_tetra(:, my_count) = vert_fold

                            my_count = my_count + 1

                            uniq_count( grd_indx_unq ) = .false.

                        end if unq_test

                    end do v_loop

                    ! To recover the vertices(v) of the tetrahedron numbered tetra_count =>         !
                    ! qv = grid_pnt_tetra( :, grid_indx_tetra( tetrahedrons(v, tetra_count) ) )     !
                    ! q1 = this%q_pnt( :, grid_indx_tetra( tetrahedrons(1, kk) ) )                  !
                    
                    ! $ !
                    tetrahedrons(:, my_tetra_count) = (/box_trans(4, 1), box_trans(4, 2), &
                                                      & box_trans(4, 3), box_trans(4, 4)/)  !(1, 5, 2, 4)

                    tetrahedrons(:, my_tetra_count+1) = (/box_trans(4, 3), box_trans(4, 4), &
                                                      & box_trans(4, 5), box_trans(4, 2)/)  !(2, 4, 6, 5)

                    tetrahedrons(:, my_tetra_count+2) = (/box_trans(4, 6), box_trans(4, 4), &
                                                      & box_trans(4, 5), box_trans(4, 2)/)  !(8, 4, 6, 5)

                    tetrahedrons(:, my_tetra_count+3) = (/box_trans(4, 6), box_trans(4, 4), &
                                                      & box_trans(4, 7), box_trans(4, 2)/)  !(8, 4, 7, 5)

                    tetrahedrons(:, my_tetra_count+4) = (/box_trans(4, 7), box_trans(4, 4), &
                                                      & box_trans(4, 8), box_trans(4, 2)/)  !(7, 4, 3, 5)

                    tetrahedrons(:, my_tetra_count+5) = (/box_trans(4, 1), box_trans(4, 4), &
                                                      & box_trans(4, 8), box_trans(4, 2)/)  !(1, 4, 3, 5)
                    ! $ !
                    my_tetra_count = my_tetra_count + 6

                end do ix_loop
            end do iy_loop
        end do iz_loop

        ! indx_map is a 1d array provides the map from unq_indx to where the grid point is actually !
        ! located (my_pos) in q_pnt_int. Here unq_indx is given by:                                 ! 
        ! unq_indx = dot_product(modulo(grid_int{folded}, mesh), mesh_mul) + 1. So the 3-integer    !
        ! grid point of a particular q-point is given by grid_int = q_pnt_int(1:3, my_pos).         ! 
        ! This can be useful to store only a single ineger (my_pos) corresponding to a grid_int     !
        ! (3-integers). Advantages: (1) Storage. (2) Helpful during the quasi momentum conservation !
        ! (q2 specifically). (3) This map (indx_map) makes the points of tetrahedrons consecutive   !
        allocate( this%indx_map(num_grid_pnt) )

        ! tetrahedrons(:, tetra_count) = (/ Four integers/), which are unq_indx and corresponds to  !
        ! 4 vertices of tetrahedron numbered tetra_count. To recover the vertices '1' of the        !
        ! tetrahedron numbered 'tetra_count': q1 = q_pnt(:, indx_map( tetrahedrons(1, tetra_count)) )!
        allocate( this%tetrahedrons(4, num_grid_pnt*6) )

        this%indx_map = grid_indx_tetra
        this%tetrahedrons = tetrahedrons

        !*! debug !*!
        !*! write(*, *) "my_count: ", my_count
        !*! write(*, *) "my_tetra_count: ", my_tetra_count
        !*! debug !*!
        !! =================================== devide the whole BZ into tetrahedrons ========================================= !!

        !! ============================== Irreducible BZ only considering time reversal symmetry ============================= !!

        ! 'map' is a 1d array of size ( mesh(1)*mesh(2)*mesh(3) ). 'map' stores a single integer corresponding to symetrically  !
        ! equivalent set of q-points ( here q and -q). The unique integer corresponding to this set of symmetrically related    !
        ! q-points is obtained from indx_map (or grid_indx_tetra). So, basically this is the position of q-point in q_pnt_int.  !
        ! Which 'qfirst' in the symmetrically related set of qs is                                                              !
        ! encountered first, indx_map of that 'qfirst' is assigned to all the qs in the set.                                    !
        ! indx_qfirst = indx_map( dot_product( modulo( grid_pnt_tetra(:, ii), mesh ), mesh_mul ) + 1 )                          !
        ! Both q and -q in the 'map' is assigned indx_qfirst here.                                                              !
        map = 0
        num_irr_q = (num_grid_pnt / 2) + 1

        q_grid_loop: do ii = 1, num_grid_pnt

            q_indx_unq = dot_product( modulo( grid_pnt_tetra(:, ii), mesh ), mesh_mul ) + 1
            indx_qfirst = grid_indx_tetra(q_indx_unq)

            ! indx_qfirst and ii are same integer, because q_point is taken from grid_pnt_tetra !

            !*! write(*, *) 'indx_qfirst = ', indx_qfirst, 'ii = ', ii

            chk_exists: if ( map(indx_qfirst) == 0 ) then

                ! If map(indx_qfirst) = 0, the corresponding q and -q is not encountered earlier. Assign the    !
                ! map(indx_qfirst) to indx_qfirst and also assign the map(indx_minusq) to same indx_qfirst      !
                map(indx_qfirst) = indx_qfirst

                qminus_indx_unq = dot_product( modulo( -1*grid_pnt_tetra(:, ii), mesh ), mesh_mul ) + 1
                indx_minusq = grid_indx_tetra(qminus_indx_unq)

                chk_non_zero: if ( indx_minusq /= indx_qfirst ) then

                    ! If map(indx_qfirst) = 0 the map(indx_minusq) must also be 0   !
                    debug: if ( map(indx_minusq) /= 0 ) then
                        write(*, 20)
                        20 FORMAT ('Error: error in irredcible BZ')
                        ERROR STOP
                    end if debug

                    map(indx_minusq) = indx_qfirst

                end if chk_non_zero

            end if chk_exists

            !*! debug !*!
            !*! write(*, 43) ii, grid_pnt_tetra(:, ii), map(ii), grid_pnt_tetra(:, map(ii)), map(map(ii))
            !*! write(*, *)
            !*! 43 FORMAT(I5, '-->', 3I3, ' | ', I5, ' | ', 3I3, ' || ', I5) 
            !*! debug !*!

        end do q_grid_loop

        !! ============================== Irreducible BZ only considering time reversal symmetry ============================= !!

        !*! write(*, *) map

        ! 'tmp_arr' is a 1d array of size ( mesh(1)*mesh(2)*mesh(3) ).                                                          ! 
        ! All the symmetrically related qs shares the same map number in the map array, So =>                                   !
        ! For a single q among the symmetrically related set of qs                                                              !
        ! (which arises first, say qfirst) it is assigned to a nonzero number (which is increasing from 1 to irr_q_count).      !
        ! The position where this non-zero number is assigned in 'tmp_arr' is given by map(my_pos of the qfirst). This will     !
        ! cover all the qs of the set of symmetrically related qs, as map(my_pos of the qfirst) = map(my_pos of any q which are !
        ! symmetrically related to qfirst).                                                                                     !
        allocate( tmp_arr(num_grid_pnt) )
        ! Initially make all the tmp_arr to 0   !
        tmp_arr = 0

        ! this%q_pnt_int(1:3, :) = grid_pnt_tetra(:, :). The fourth value this%q_pnt_int(4, :) is a integer which corresponds   !
        ! to the position of the symmetrically equivalent q-point of this%irr_q_pnt. So, this%q_pnt_int(1:3, ii) and            !
        ! this%irr_q_int( 1:3, this%q_pnt_int(4, ii) ) are symmetrically related. All the properties of grid_pnt_tetra is       !
        ! available to this%q_pnt_int.                                                                                          !
        ! The fifth value this%q_pnt_int(5, :) indicates whether it is in irr_q_pnt or not. If this fifth value is 1 then it    !
        ! exists in irr_q_pnt and if it is 0 then it does not exist in irr_q_pnt                                                !
        allocate( this%q_pnt_int(5, num_grid_pnt) )
        this%q_pnt_int = 0

        ! this%irr_q_int(1:3, :) correspond to the coordinate of irreducible q-point. The fourth value this%irr_q_int(4, :)     !
        ! corresponds to the degeneracy, i.e., number of q-points that are symmetrically equivalent to this irr_q_pnt           !
        allocate( this%irr_q_int(4, num_irr_q) ) 
        this%irr_q_int = 0

        allocate( this%q_pnt(3, num_grid_pnt) )
        allocate( this%irr_q(3, num_irr_q) ) 

        this%q_pnt = 0.0_dp
        this%irr_q = 0.0_dp

        irr_q_count = 0        
        grid_pnt_loop: do ii = 1, num_grid_pnt

            grid_pnt = grid_pnt_tetra(1:3, ii)

            qvec = matmul( sys%G, &
                        & ((dble(grid_pnt) + dble(shift)*0.5_dp) / dble(mesh)) )

            this%q_pnt_int(1:3, ii) = grid_pnt

            this%q_pnt(:, ii) =  qvec

            ! If tmp_arr( map(my_pos or ii) ) is zero then the representative q from the set of symmetrically realted qs    !
            ! is not recorded yet in irr_q_int (and irr_q). Also make the tmp_arr( map(my_pos or ii) ) to the current       !
            ! count of irreducible q. So if later any q which is from the set of symmetrically related qs appear it will    !
            ! encounter a non-zero values. All the symmetrically related q will encounter this non-zero value as            !
            ! map(my_pos of the qfirst) = map(my_pos of any q which are symmetrically related to qfirst).                   !
            ! So here tmp_arr( map(my_pos of q) ) and tmp_arr( map(my_pos of -q) ) have the nonzero value as                !
            ! map(my_pos of q) = map(my_pos of -q)                                                                          !
            ! Here ii and grid_indx_tetra(q_indx_unq) are same integer, as grid_point is taken from grid_pnt_tetra          ! 
            my_num = tmp_arr( map(ii) )

            chk_irr_q_covered: if ( (my_num == 0) .and. (irr_q_count < num_irr_q) ) then

                irr_q_count = irr_q_count + 1
                tmp_arr( map(ii) ) = irr_q_count

                my_num = irr_q_count ! Reuse of my_num !

                debug2: if ( any(this%irr_q_int(:, my_num) > 0) ) then
                    write(*, 32)
                    ERROR STOP
                end if debug2

                this%irr_q_int(1:3, my_num) = this%q_pnt_int(1:3, ii)
                this%irr_q(:, my_num) = qvec

                ! ***** !
                this%q_pnt_int(5, ii) = 1
                ! ***** !

            else if ( (my_num == 0) .and. (irr_q_count >= num_irr_q) ) then chk_irr_q_covered

                write(*, 32)
                ERROR STOP

            end if chk_irr_q_covered

            this%q_pnt_int(4, ii) = my_num
            this%irr_q_int(4, my_num) = this%irr_q_int(4, my_num) + 1

        end do grid_pnt_loop
        32 FORMAT('Error: Something going wrong in irreducible q point generation')

        !*! debug !*!

        !*! write(*, *) 
        !*! write(*, *) 'Debug'
        !*! do ii = 1, num_grid_pnt

        !*!     write(*, 22) this%q_pnt_int(1:3, ii) , this%q_pnt_int(5, ii), this%q_pnt_int(4, ii), &
        !*!                & this%irr_q_int(1:3, this%q_pnt_int(4, ii)), &
        !*!                & this%irr_q_int(4, this%q_pnt_int(4, ii))
        !*!     
        !*!     22 FORMAT('(', 2(I3, ', '), I3, ') ', '[ ', I1, ' ]', ' ( ', I0, ' ) --> ', '(', 2(I3, ', '), I3, ') ', '[', I0, ']',/)

        !*! end do

        !*! write(*, *) maxval( this%q_pnt_int(4, :) )
        !*! write(*, *) sum( this%irr_q_int(4, :) )
        !*! write(*, *) sum( this%q_pnt_int(5, :) )

        !*! write(*, *) 'Volume'

        !*! do kk = 1, (num_grid_pnt*6)
        !*!     q1 = this%q_pnt( :, grid_indx_tetra( tetrahedrons(1, kk) ) ) !a
        !*!     q2 = this%q_pnt( :, grid_indx_tetra( tetrahedrons(2, kk) ) ) !b
        !*!     q3 = this%q_pnt( :, grid_indx_tetra( tetrahedrons(3, kk) ) ) !c
        !*!     q4 = this%q_pnt( :, grid_indx_tetra( tetrahedrons(4, kk) ) ) !d

        !*!     c = cross((q2-q4), (q3-q4))
        !*!     vol = dot_product( (q1-q4), c)  

        !*!     write(*, *) kk, dabs(vol)/6.0_dp

        !*! end do

        !*! debug !*!

        deallocate( tmp_arr )

    end subroutine make_tetrahedron

    function cross(a, b)

        implicit none
        real(dp), intent(in) :: a(:), b(:)
        real(dp), dimension(3) :: cross

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)

    end function cross

    subroutine clean_Qdata( this )

        ! ** Finalizer ** !

        implicit none

        CLASS(q_points_data)                    :: this
        !TYPE(q_points_data)                     :: this

        ! ============================= Local Variable ============================= !

        integer                                 :: istat

        ! ============================= Local Variable ============================= !

        deallocate( this%indx_map, STAT=istat )
        deallocate( this%q_pnt_int, this%irr_q_int, this%tetrahedrons, STAT=istat )
        deallocate( this%q_pnt, this%irr_q, STAT=istat )

    end subroutine clean_Qdata

end module Irr_q_point

