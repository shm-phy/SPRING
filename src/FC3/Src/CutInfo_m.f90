
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

module CutInfo

    use kinds,                  only : dp
    use unit_cell,              only : cell
    use IndxChng,               only : fold_indx

    implicit none

    private

    type, public        :: CellBasis

        integer, dimension(:,:), allocatable            :: cbIndx

    end type CellBasis

    type, public        :: CutOff

        type(CellBasis), dimension(:), allocatable      :: atmIndx
        integer, dimension(:), allocatable              :: numAtom
        real(dp), dimension(:), allocatable             :: CutRad

    contains
        procedure, public, pass                         :: set_CutOff
        procedure, private, nopass                      :: SortUnique

    end type CutOff

contains

    subroutine set_CutOff(this, sys, near_neigh)

        implicit none

        class(CutOff)                               :: this

        type(cell), intent(in)                      :: sys
        integer, dimension(:), intent(in)           :: near_neigh

        ! ===================================== Local variable ========================================= !

        integer, dimension(3)                       :: SupCell, N2flat, N2fold
        integer                                     :: Nbasis, cnt
        real(dp), dimension(3)                      :: atm1cord, atm2cord
        real(dp), dimension(:), allocatable         :: atmDist
        real(dp), dimension(:), allocatable         :: outSrtUnq
        real(dp)                                    :: dist
        integer, dimension(:,:), allocatable        :: cell_basis

        integer                                     :: mu, nx, ny, nz, nu !, ii

        ! ===================================== Local variable ========================================= !

        SupCell = sys%sup_cell
        Nbasis = sys%natm

        allocate( this%CutRad(Nbasis) )
        
        allocate( atmDist( Nbasis*product(SupCell) ) )

        write(*, *)
        out_mu_loop: do mu = 1, Nbasis

            cnt = 1
            atmDist(:) = 0.0_dp

            nx_loop: do nx = 0, (SupCell(1) - 1)
                ny_loop: do ny = 0, (SupCell(2) - 1)
                    nz_loop: do nz = 0, (SupCell(3) - 1)

                        N2flat = (/nx, ny, nz/)
                        N2fold = fold_indx(SupCell, N2flat)

                        atm1cord = matmul( sys%latvec, &
                                        & (sys%basis_frac(:, mu)) )

                        basis_loop: do nu = 1, Nbasis

                            atm2cord = matmul( sys%latvec, &
                                            & (dble( N2fold ) + sys%basis_frac(:, nu) ) )

                            atmDist(cnt) = norm2( (atm2cord-atm1cord) )
                            cnt = cnt + 1

                        end do basis_loop

                    end do nz_loop
                end do ny_loop
            end do nx_loop

            call SortUnique( atmDist, outSrtUnq )

            CutOffOutCheck: if ( near_neigh(mu) >= (size( outSrtUnq ) - 1) ) then
                write(*, 420)  mu, near_neigh(mu), (size( outSrtUnq ) - 1)
                420 FORMAT( "ERROR: CutOff lies outside box. [near_neigh(", I3, ") = ", I5, "] >= [SIZE(outSrtUnq) = ", I5, "]" )
                write(*, 430)
                430 FORMAT( "Cannot proceed, decrease near_neigh! (CutOff distance should not exceed L_box/2)" )
                ERROR STOP
            end if CutOffOutCheck

            this%CutRad(mu) = 0.5_dp * ( outSrtUnq(near_neigh(mu)+2) + outSrtUnq(near_neigh(mu)+1) )

            write(*, 110) mu, near_neigh(mu), this%CutRad(mu)
            110 FORMAT("[Basis Index = ", I4, "] => Cut-off radius for ", I4, "th nearest neighbor is:", F12.6, " Ang")

            !*! debug !*!
            !*! write(*, *) 
            !*! write(*, *) atmDist
            !*! write(*, *)
            !*! write(*, *) outSrtUnq
            !*! debug !*!

            deallocate( outSrtUnq )

        end do out_mu_loop

        allocate( cell_basis(4,Nbasis*product(SupCell)) )
        allocate( this%numAtom(Nbasis) )
        allocate( this%atmIndx(Nbasis) )

        write(*, *)
        out_mu_loop2: do mu = 1, Nbasis

            cnt = 0

            cell_basis(:,:) = 0
            nx_loop2: do nx = 0, (SupCell(1) - 1)
                ny_loop2: do ny = 0, (SupCell(2) - 1)
                    nz_loop2: do nz = 0, (SupCell(3) - 1)

                        N2flat = (/nx, ny, nz/)
                        N2fold = fold_indx(SupCell, N2flat)

                        atm1cord = matmul( sys%latvec, &
                                        & (sys%basis_frac(:, mu)) )

                        basis_loop2: do nu = 1, Nbasis

                            atm2cord = matmul( sys%latvec, &
                                            & (dble( N2fold ) + sys%basis_frac(:, nu) ) )

                            dist = norm2( (atm2cord-atm1cord) )

                            cut_rad_chk: if ( dist < this%CutRad(mu) ) then

                                cnt = cnt + 1
                                cell_basis(1:3, cnt) = N2fold
                                cell_basis(4, cnt) = nu

                            end if cut_rad_chk

                        end do basis_loop2

                    end do nz_loop2
                end do ny_loop2
            end do nx_loop2

            this%numAtom(mu) = cnt
            allocate( this%atmIndx(mu)%cbIndx(4, cnt) )
            this%atmIndx(mu)%cbIndx(:, :) = cell_basis(:, 1:cnt)

            write(*, 115) mu, near_neigh(mu), this%numAtom(mu)
            115 FORMAT("[Basis Index = ", I4, "] => Number of atoms with in ", I4, "th nearest neighbor is:", I4)

            !*! debug !*!
            !*! write(*, *) 
            !*! write(*, *) "Cell Basis Index: "

            !*! do ii = 1, cnt
            !*!     write(*, 125) this%atmIndx(mu)%cbIndx(:, ii)
            !*! end do

            !*! 125 FORMAT("[", 3(I3, ', '), I3, "]")
            !*! debug !*!

        end do out_mu_loop2

    end subroutine set_CutOff


    subroutine SortUnique( atmDist, outSrtUnq )

        implicit none

        real(dp), dimension(:), intent(in)                      :: atmDist
        real(dp), dimension(:), allocatable, intent(out)        :: outSrtUnq

        ! =============================== Local variables =============================== !

        real(dp), dimension(:), allocatable                     :: uniq
        real(dp)                                                :: min_val, max_val, &
                                                                 & accu

        integer                                                 :: arr_len, i

        ! =============================== Local variables =============================== !

        arr_len = size( atmDist )
        allocate( uniq(arr_len) )

        min_val = minval(atmDist) - 0.7425_dp
        max_val = maxval(atmDist)

        accu = 1.0E-4_dp
        i = 0
        do while ( min_val < max_val )

            i = i + 1
            min_val = minval( atmDist, mask= ( (atmDist-min_val)>accu ) )

            uniq(i) = min_val

        end do 

        allocate( outSrtUnq(i), source=uniq(1:i) )

    end subroutine SortUnique

end module CutInfo

