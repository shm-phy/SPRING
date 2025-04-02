
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


module Helper

    use kinds,              only : dp
    use constants,          only : EPS_nzero
    use unit_cell,          only : cell
    use CutInfo,            only : CutOff, CellBasis

    implicit none

    private

    public          :: TransformAtoms, get_atom_num, kron_prod

contains

    subroutine TransformAtoms(sys, atm1cb, atm2cb, atm3cb, atm4cb, &
                            & Rot, Trans)

        implicit none

        integer, parameter                              :: fc_ord = 4

        type(cell), intent(in)                          :: sys

        integer, dimension(4), intent(inout)            :: atm1cb, atm2cb, atm3cb, &
                                                         & atm4cb
        real(dp), dimension(3,3), intent(in)            :: Rot
        real(dp), dimension(3), intent(in)              :: Trans

        ! ==================================== Local Variables ==================================== !

        integer, dimension(4, fc_ord)                   :: atmcb
        real(dp), dimension(3, fc_ord)                  :: atom, atomc, atm_frac
        integer, dimension(3, fc_ord)                   :: atm_int

        integer, dimension(fc_ord)                      :: basis_indx
        integer, dimension(3)                           :: CentralCell

        integer                                         :: Nbasis, cnt
        integer                                         :: nn, ii, nb

        ! ==================================== Local Variables ==================================== !

        atmcb(:, 1) = atm1cb
        atmcb(:, 2) = atm2cb
        atmcb(:, 3) = atm3cb
        atmcb(:, 4) = atm4cb

        atm_loop1: do nn = 1, fc_ord

            !*! write(*, *)
            !*! write(*, *)

            atom(:, nn) = dble( atmcb(1:3, nn) ) + sys%basis_frac(:, atmcb(4, nn))
            atom(:, nn) = matmul( atom(:, nn), Rot ) + Trans

            atomc(:, nn) = atom(:, nn)

            atm_int(:, nn) = int( atom(:, nn) )
            atm_frac(:, nn) = atom(:, nn) - dble( atm_int(:, nn) )

            !*! write(*, *) 
            !*! write(*, *) atomc(:, nn)

            !*! write(*, *) "************************************"
            !*! write(*, *) atm_int(:, nn)
            !*! write(*, *) atm_frac(:, nn)

            comp_loop: do ii = 1, 3

                if (  ( (atm_int(ii, nn)) <= 0 )  .and. (atm_frac(ii, nn) < 0.0_dp) &
                    & .and. ( dabs(atm_frac(ii, nn)) >= EPS_nzero ) ) then

                    atom(ii, nn) = atom(ii, nn) - 1.0_dp

                end if

            end do comp_loop

            atm_int(:, nn) = int( atom(:, nn) )
            atm_frac(:, nn) = atom(:, nn) - dble( atm_int(:, nn) )

            !*! write(*, *) "************************************"
            !*! write(*, *) atm_int(:, nn)
            !*! write(*, *) atm_frac(:, nn)

            comp_loop2: do ii = 1, 3

                if ( atm_frac(ii, nn) < 0.0_dp ) then

                    atm_frac(ii, nn) = 1.0_dp + atm_frac(ii, nn)

                end if

            end do comp_loop2

            !*! write(*, *) "************************************"
            !*! write(*, *) atm_int(:, nn)
            !*! write(*, *) atm_frac(:, nn)

            if ( any ( dabs( dble(atm_int(:, nn)) + atm_frac(:, nn) - &
                     & atomc(:, nn) ) > EPS_nzero ) ) then  

                write(*, 10) 
                write(*, 20) atm_int(:, nn), atm_frac(:, nn), atomc(:, nn)
                STOP

            end if

        end do atm_loop1

        CentralCell(:) = atm_int(:, 1)
        Nbasis = sys%natm
        
        atm_loop2: do nn = 1, fc_ord

            atm_int(:, nn) = atm_int(:, nn) - CentralCell(:)

            cnt = 0
            basis_loop: do nb = 1, Nbasis

                basis_chk: if ( all( dabs( atm_frac(:, nn) - &
                                         & sys%basis_frac(:, nb) ) < EPS_nzero ) ) then

                    cnt = cnt + 1
                    basis_indx( nn ) = nb

                end if basis_chk

            end do basis_loop

            err_chk: if ( (cnt == 0) .or. (cnt > 1) ) then

                write(*, 10)
                write(*, 30) cnt
                STOP

            end if err_chk

        end do atm_loop2

        atm1cb(1:3) = atm_int(:, 1)
        atm1cb(4) = basis_indx(1)

        atm2cb(1:3) = atm_int(:, 2)
        atm2cb(4) = basis_indx(2)

        atm3cb(1:3) = atm_int(:, 3)
        atm3cb(4) = basis_indx(3)

        atm4cb(1:3) = atm_int(:, 4)
        atm4cb(4) = basis_indx(4)

        10 FORMAT("ERROR")
        20 FORMAT("[", 2(I3, ', '), I3, "] + [", 2(F6.3, ', '), F6.3, "] /= [", 2(F6.3, ', '), F6.3, "]")
        30 FORMAT("0 or more than 1 basis batch: ", I4)

    end subroutine TransformAtoms


    Function get_atom_num(mu_dep, dep_N_bs, Cut_Off) Result(atmnum)

        implicit none

        integer, intent(in)                     :: mu_dep
        integer, dimension(4), intent(in)       :: dep_N_bs
        type(CutOff), intent(in)                :: Cut_Off

        integer                                 :: atmnum !Result

        ! ====================================== Local Variables ====================================== !

        integer                         :: Natm, cnt
        integer                         :: na

        ! ====================================== Local Variables ====================================== !

        ! * Initialize * !
        atmnum = 0
        ! * Initialize * !

        Natm = Cut_Off%numAtom(mu_dep)
        cnt = 0
        atm_loop: do na = 1, Natm

            exist_chk: if ( all( (Cut_Off%atmIndx(mu_dep)%cbIndx(:, na) - dep_N_bs) == 0 ) ) then

                atmnum = na
                cnt = cnt + 1

            end if exist_chk

        end do atm_loop

        no_exist: if ( cnt == 0 ) then
            atmnum = 0

        else if ( cnt > 1 ) then no_exist
            write(*, 24) cnt
            STOP

        end if no_exist

        24 FORMAT("ERROR(in get_atom_num): Map to multiple atoms: ", I4)

    end Function get_atom_num


    subroutine kron_prod(invec1, invec2, outvec)

        implicit none

        real(dp), dimension(:), intent(in)                  :: invec1, invec2
        real(dp), dimension(:), allocatable, intent(out)    :: outvec

        ! ================================== Local Variables ================================== !

        integer                                             :: vec_size1, vec_size2
        integer                                             :: vi, strtp, endp

        ! ================================== Local Variables ================================== !

        vec_size1 = size( invec1 )
        vec_size2 = size( invec2 )

        allocate( outvec(vec_size1*vec_size2) )
        outvec = 0.0_dp

        strtp = 1
        endp = vec_size2

        invec1_loop: do vi = 1, vec_size1

            outvec(strtp : endp) = invec1(vi) * invec2(:)

            strtp = endp + 1
            endp = endp + vec_size2

        end do invec1_loop

    end subroutine kron_prod

end module Helper

