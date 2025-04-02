
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

module Symmetry

    use kinds,              only : dp
    use unit_cell,          only : cell
    use spglib_f08,         only : spg_get_symmetry
    use mklWrap,            only : dinverse

    implicit none

    private

    type, public    :: SpaceGr_typ

        integer                                             :: Nsym
        real(dp), dimension(:,:,:), allocatable             :: R
        real(dp), dimension(:,:), allocatable               :: T
        real(dp), dimension(:,:,:), allocatable             :: Rxyz

        contains
            procedure, public, pass                         :: set_SpaceGrSPGlib

    end type SpaceGr_typ

contains

    subroutine set_SpaceGrSPGlib(this, sys, symprec)

        implicit none
        integer, parameter                                  :: max_num_sym=192

        class(SpaceGr_typ)                                  :: this

        type(cell), intent(in)                              :: sys
        real(dp), intent(in)                                :: symprec

        ! ===================================== Local variable ========================================= !

        integer                                             :: num_atom
        real(dp), dimension(3, 3)                           :: lattice_t
        real(dp), dimension(:, :), allocatable              :: positions
        integer, dimension(:), allocatable                  :: atom_types

        integer, dimension(3, 3, max_num_sym)               :: rotations
        real(dp), dimension(3, max_num_sym)                 :: translations

        real(dp), dimension(3, 3)                           :: latc, latc_inv, sym_mat, &
                                                             & sym_matT, symmatT_latcinv
        integer                                             :: ns

        integer                                             :: num_symm

        ! ===================================== Local variable ========================================= !

        num_atom = sys%natm

        lattice_t = transpose( sys%latvec )

        allocate( positions(3, num_atom) )
        positions = sys%basis_frac

        allocate( atom_types(num_atom) )
        atom_types = sys%typ_indx

        num_symm = spg_get_symmetry( rotations, translations, max_num_sym, lattice_t, &
                                   & positions, atom_types, num_atom, symprec )

        this%Nsym = num_symm

        img_chk: if ( this_image() == 1 ) then
            write(*, *)
            write(*, 10) num_symm, symprec
        end if img_chk

        allocate( this%R(3, 3, num_symm) )
        this%R(:,:,:) = dble( rotations(:,:,1:num_symm) )

        allocate( this%T(3, num_symm) )
        this%T(:,:) = translations(:, 1:num_symm)

        allocate( this%Rxyz(3, 3, num_symm) )

        latc = sys%latvec 
        latc_inv = latc
        call dinverse(latc_inv)

        sym_op_loop: do ns = 1, num_symm
            
            sym_mat = this%R(:, :, ns)
            sym_matT = transpose( sym_mat )

            symmatT_latcinv = matmul(sym_matT, latc_inv)

            this%Rxyz(:, :, ns) = matmul( latc, symmatT_latcinv )

        end do sym_op_loop

        10 FORMAT("Number of Space Group Symmetry operations: ", I5, " for symmetry precision: ", E10.3)

    end subroutine set_SpaceGrSPGlib

end module Symmetry

