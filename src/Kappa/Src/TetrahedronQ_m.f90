
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


module TetrahedronQ

    use kinds,          only : dp
    use unit_cell,      only : cell
    use IndxChng,       only : fold_indx, MapAwayFromBZboundary
    use Helper,         only : unique_sort, FindPositionAndMultiplicity

    use spglib_f08,     only : spg_get_ir_reciprocal_mesh

    implicit none
    private

    ! The datatype q_points_data contains all the q-points information. !
    type, public        :: q_tetrahedron

        real(dp), allocatable, dimension(:,:)           :: q_pnt
        real(dp), allocatable, dimension(:,:)           :: irr_q_pnt

        integer, allocatable, dimension(:,:)            :: q_pnt_int
        integer, allocatable, dimension(:,:)            :: irr_q_pnt_int
        integer, allocatable, dimension(:,:)            :: tetrahedrons

        integer, allocatable, dimension(:)              :: pos
        integer, allocatable, dimension(:)              :: multiplicity

        contains
            procedure, public, pass                     :: make_tetrahedron, &
                                                         & FindIrrQ, clean_Qdata
            !FINAL                                       :: clean_Qdata

    end type q_tetrahedron

contains

    subroutine make_tetrahedron(this, sys, mesh, shift, ZeroCentered)

        implicit none

        class(q_tetrahedron)                                        :: this

        type(cell), intent(in)                                      :: sys
        integer, dimension(3), intent(in)                           :: mesh, shift
        logical, intent(in)                                         :: ZeroCentered

        !===================================== Local variable =========================================!

        real(dp), dimension(3)                                      :: qvec
        real(dp), dimension(3, 3)                                   :: Recp_cell_unit
        real(dp), dimension( 3, mesh(1)*mesh(2)*mesh(3) )           :: q_flaot

        integer, dimension(3)                                       :: mesh_mul

        !! ******************************** Tetrahedron ******************************** !!

        integer, dimension(3, 8)                                    :: box
        integer, dimension(4, 8)                                    :: box_trans
        integer, dimension(3)                                       :: indx, vec_trans, &
                                                                     & vert, vert_fold

        integer, dimension( 3, mesh(1)*mesh(2)*mesh(3) )            :: grid_pnt_tetra
        integer, dimension( 4, mesh(1)*mesh(2)*mesh(3)*6 )          :: tetrahedrons

        integer                                                     :: grd_indx_unq, my_count, &
                                                                     & my_tetra_count, num_grid_pnt
        integer                                                     :: i1, ix, iy, iz, v

        logical, dimension( mesh(1)*mesh(2)*mesh(3) )               :: uniq_count

        !! ******************************** Tetrahedron ******************************** !!
                                                                     
        !===================================== Local variable =========================================!

        do i1 = 1, 3
            Recp_cell_unit(:, i1) = sys%G(:, i1) / dble(mesh(i1))
        end do

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

                    if ( ZeroCentered ) then
                        call fold_indx(mesh, indx, vec_trans)
                    else
                        vec_trans = indx
                    end if

                    v_loop: do v = 1, 8

                        vert = box(:, v) + vec_trans

                        if ( ZeroCentered ) then
                            call fold_indx(mesh, vert, vert_fold)
                        else
                            vert_fold = MapAwayFromBZboundary(mesh, vert)
                        end if

                        grd_indx_unq = dot_product( modulo( vert_fold, mesh ), mesh_mul ) + 1

                        !*! degub !*!
                        debg_chk: if ( grd_indx_unq > num_grid_pnt ) then

                            write(*, 40) vert_fold, grd_indx_unq
                            40 FORMAT("ERROR: grd_indx_unq out of bound: (", 2(I4, '  '), I4, " ) => ", I4)
                            ERROR STOP

                        end if debg_chk
                        !*! degub !*!

                        box_trans(1:3, v) = vert_fold
                        box_trans(4, v) = grd_indx_unq

                        unq_test: if ( uniq_count(grd_indx_unq) ) then

                            grid_pnt_tetra(:, grd_indx_unq) = vert_fold

                            qvec = matmul( Recp_cell_unit, (dble(vert_fold) + dble(shift)*0.5_dp) )

                            q_flaot(:, grd_indx_unq) = qvec

                            my_count = my_count + 1
                            uniq_count( grd_indx_unq ) = .false.

                        end if unq_test

                    end do v_loop

                    ! To recover the vertices(v) of the tetrahedron numbered 'tetra_count' =>       !
                    ! qv = grid_pnt_tetra( :,  tetrahedrons(v, tetra_count) )                       !
                    ! q1 = this%q_pnt( :, tetrahedrons(1, kk) )                                     !
                    
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

        allocate( this%q_pnt(3, num_grid_pnt) )

        ! tetrahedrons(:, tetra_count) = (/ Four integers/), which are unq_indx and corresponds to  !
        ! 4 vertices of tetrahedron numbered tetra_count. To recover the vertices '1' of the        !
        ! tetrahedron numbered 'tetra_count': q1 = q_pnt(:, tetrahedrons(1, tetra_count) )          !
        allocate( this%tetrahedrons(4, num_grid_pnt*6) )

        ! this%q_pnt_int(1:3, :) = grid_pnt_tetra(:, :)   !
        allocate( this%q_pnt_int(3, num_grid_pnt) )

        this%q_pnt = q_flaot
        this%tetrahedrons = tetrahedrons
        this%q_pnt_int = grid_pnt_tetra

        !*! debug !*!
        debug2: if ( ( (my_count-1) /= num_grid_pnt ) .or. &
                   & ( (my_tetra_count-1) /= num_grid_pnt*6 ) ) then
            
            write(*, 111)
            111 FORMAT("ERROR: In make_tetrahedron")

        end if debug2
        !*! debug !*!
        !! =================================== devide the whole BZ into tetrahedrons ========================================= !!

    end subroutine make_tetrahedron


    subroutine FindIrrQ( this, sys, mesh, shift, time_rev, symprec )

        implicit none

        class(q_tetrahedron)                                        :: this

        type(cell), intent(in)                                      :: sys
        integer, dimension(3), intent(in)                           :: mesh, shift
        integer, intent(in)                                         :: time_rev
        real(dp), intent(in)                                        :: symprec

        !================================= Local variables ===================================!

        real(dp), dimension(3, 3)                                   :: lattice_t, Recp_cell_unit
        real(dp), dimension(3)                                      :: irr_q_float

        integer, dimension( 3, mesh(1)*mesh(2)*mesh(3) )            :: grid_pnt_int
        integer, dimension( mesh(1)*mesh(2)*mesh(3) )               :: map, pos

        integer, dimension(3)                                       :: grid_address
        integer                                                     :: num_irr_q, num_grid_pnt, &
                                                                     & ii, grid_indx

        integer, dimension(:), allocatable                          :: unique_irrq, multiplicity

        !================================= Local variables ===================================!

        num_grid_pnt = product( mesh )

        do ii = 1, 3
            Recp_cell_unit(:, ii) = sys%G(:, ii) / dble(mesh(ii))
        end do

        !=================================== spglib call =====================================!

        lattice_t = transpose( sys%latvec )

        num_irr_q = spg_get_ir_reciprocal_mesh( grid_pnt_int, map, mesh, shift, &
                                              & time_rev, lattice_t, sys%basis_frac, &
                                              & sys%typ_indx, sys%natm, symprec )

        !=================================== spglib call =====================================!

        call unique_sort( map, unique_irrq )
        unique_irrq = unique_irrq + 1
        map = map + 1

        !write(*, *) ( (grid_pnt_int - this%q_pnt_int) /= 0 )
        !write(*, *) .not. all( (grid_pnt_int - this%q_pnt_int) == 0 )

        debug: if ( any( (grid_pnt_int - this%q_pnt_int) /= 0 ) .or. &
                  & ( size(unique_irrq) /= num_irr_q ) ) then

            write(*, 10)
            ERROR STOP
            10 FORMAT(" ERROR: In FindIrrQ ")

        end if debug

        if ( this_image() == 1 ) then

            Write(*, 20) mesh, num_irr_q
            20 FORMAT( "Number of irreducible q-points for ( ", 2(I4, ' X'), I4, " ) mesh is: ", I4 )

        end if

        allocate( this%irr_q_pnt(3, num_irr_q) )
        allocate( this%irr_q_pnt_int(3, num_irr_q) )

        irr_q_loop: do ii = 1, num_irr_q

            grid_indx = unique_irrq( ii )
            grid_address = grid_pnt_int(:, grid_indx)
            irr_q_float = matmul( Recp_cell_unit, (dble(grid_address) + dble(shift)*0.5_dp) )

            this%irr_q_pnt(:, ii) = irr_q_float
            this%irr_q_pnt_int(:, ii) = grid_address

        end do irr_q_loop

        allocate( multiplicity(num_irr_q) )
        call FindPositionAndMultiplicity( num_grid_pnt, num_irr_q, map, unique_irrq, &
                                        & pos, multiplicity )

        allocate( this%pos(num_grid_pnt) )
        this%pos(:) = pos(:)

        !call move_alloc( multiplicity, this%multiplicity )
        allocate( this%multiplicity(num_irr_q) )
        this%multiplicity(:) = multiplicity(:)

        !** debug **!
        !-! do ii = 1, num_grid_pnt
        !-!     write(*, 50) ii, this%q_pnt_int(:, ii), pos(ii), this%irr_q_pnt_int(:, pos(ii)), multiplicity( pos(ii) )
        !-!     50 FORMAT( "[", I6, "]->( ", 2(I4, ', '), I4, " ) ==> [", I6, "]->( ", 2(I4, ', '), I4, " ) ", "{ ", I6, " }" )
        !-! end do
        !-! write(*, *) sum( multiplicity ), num_grid_pnt
        !-! write(*, *)
        !-! do ii = 1, num_irr_q
        !-!     write(*, 45) this%irr_q_pnt_int(:, ii), this%irr_q_pnt(:, ii)
        !-!     45 FORMAT(2(I3, ' '), I3, ' --> ', 2(G15.7, ' '), G15.7)
        !-! end do
        !-! write(*, *)
        !** debug **!

        deallocate( unique_irrq, multiplicity )

    end subroutine FindIrrQ


    subroutine clean_Qdata( this )

        ! ** Finalizer ** !

        implicit none

        CLASS(q_tetrahedron)                    :: this
        !TYPE(q_tetrahedron)                     :: this

        ! ============================= Local Variable ============================= !

        integer                                 :: istat

        ! ============================= Local Variable ============================= !

        deallocate( this%q_pnt_int, this%irr_q_pnt_int, this%tetrahedrons, STAT=istat )
        deallocate( this%q_pnt, this%irr_q_pnt, STAT=istat )

    end subroutine clean_Qdata

end module TetrahedronQ

