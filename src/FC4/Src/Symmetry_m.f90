
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

    type, public    :: PntPer_typ

        integer                                             :: Nsym
        real(dp), dimension(:,:,:), allocatable             :: R
        real(dp), dimension(:,:), allocatable               :: T
        real(dp), dimension(:,:,:), allocatable             :: Rxyz
        integer, dimension(:,:), allocatable                :: Per

        contains
            procedure, public, pass                         :: set_PntPerSPGlib

    end type PntPer_typ

contains

    subroutine set_PntPerSPGlib(this, sys, symprec)

        implicit none
        integer, parameter                                  :: max_num_sym=192

        class(PntPer_typ)                                   :: this

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

        !*! debug !*!
        !*! integer                                             :: ii, jj
        !*! debug !*!

        ! ===================================== Local variable ========================================= !

        num_atom = sys%natm

        lattice_t = transpose( sys%latvec )

        allocate( positions(3, num_atom) )
        positions = sys%basis_frac

        allocate( atom_types(num_atom) )
        atom_types = sys%typ_indx

        num_symm = spg_get_symmetry( rotations, translations, max_num_sym, lattice_t, &
                                   & positions, atom_types, num_atom, symprec)

        this%Nsym = num_symm

        write(*, *)
        write(*, 10) num_symm, symprec

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

        !*! debug !*!
        !*! sym_op_loop1: do ii = 1, num_symm

        !*!    write(*, 20) ii
        !*!    write(*, *)
        !*!    write(*, 25)

        !*!    do jj = 1, 3
        !*!        write(*, 30) this%R(:, jj, ii)
        !*!    end do

        !*!    write(*, *) 
        !*!    write(*, 35)
        !*!    write(*, 40) this%T(:, ii)

        !*!    write(*, *)

        !*! end do sym_op_loop1

        !*! 20 FORMAT("Space group number: ", I5)
        !*! 25 FORMAT("Rotation: ")
        !*! 30 FORMAT("[", 2(F5.2, ', '), F5.2, "]")
        !*! 35 FORMAT("Translation: ")
        !*! 40 FORMAT("[", 2(F6.5, ', '), F6.3, "]")
        !*! debug !*!

        allocate( this%Per(4, 24) )

        this%Per(:, 1) = (/1, 2, 3, 4/)
        this%Per(:, 2) = (/1, 2, 4, 3/)
        this%Per(:, 3) = (/1, 3, 2, 4/)
        this%Per(:, 4) = (/1, 3, 4, 2/)
        this%Per(:, 5) = (/1, 4, 2, 3/)
        this%Per(:, 6) = (/1, 4, 3, 2/)
        this%Per(:, 7) = (/2, 1, 3, 4/)
        this%Per(:, 8) = (/2, 1, 4, 3/)
        this%Per(:, 9) = (/2, 3, 1, 4/)
        this%Per(:, 10) = (/2, 3, 4, 1/)
        this%Per(:, 11) = (/2, 4, 1, 3/)
        this%Per(:, 12) = (/2, 4, 3, 1/)
        this%Per(:, 13) = (/3, 1, 2, 4/)
        this%Per(:, 14) = (/3, 1, 4, 2/)
        this%Per(:, 15) = (/3, 2, 1, 4/)
        this%Per(:, 16) = (/3, 2, 4, 1/)
        this%Per(:, 17) = (/3, 4, 1, 2/)
        this%Per(:, 18) = (/3, 4, 2, 1/)
        this%Per(:, 19) = (/4, 1, 2, 3/)
        this%Per(:, 20) = (/4, 1, 3, 2/)
        this%Per(:, 21) = (/4, 2, 1, 3/)
        this%Per(:, 22) = (/4, 2, 3, 1/)
        this%Per(:, 23) = (/4, 3, 1, 2/)
        this%Per(:, 24) = (/4, 3, 2, 1/)

        10 FORMAT("Number of Space Group Symmetry operations: ", I5, " for symmetry precision: ", E10.3)

    end subroutine set_PntPerSPGlib

end module Symmetry

