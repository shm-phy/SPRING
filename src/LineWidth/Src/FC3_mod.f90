
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

module FC3_mod

    use kinds,          only : dp
    use unit_cell,      only : cell
    use hdf5_wrap,      only : r3rdFC

    implicit none
    private

    type, public    ::  FC3type
        
        integer, dimension(:), allocatable                  :: atmNum
        integer, dimension(:,:,:), allocatable              :: atmIndx
        integer, dimension(:,:), allocatable                :: basisTyp
        real(dp), dimension(:,:,:), allocatable             :: atmPosR
        real(dp), dimension(:,:,:,:,:,:), allocatable       :: fC

    contains
        procedure, public, pass         :: set_FC3, set_Zero

        GENERIC     ::  ASSIGNMENT(=)   =>  FC3_assignment
        GENERIC     ::  OPERATOR(+)     =>  FC3_add
        procedure, private, pass        :: FC3_assignment, FC3_add

    end type FC3type

contains

    subroutine set_FC3(this, sys, filename)

        implicit none

        class(FC3type)                                      :: this
        type(cell), intent(in)                              :: sys
        character(len=*), intent(in)                        :: filename

        ! =========================== Local variables =========================== !
        integer, dimension(:), allocatable                  :: atmNum3
        integer, dimension(:,:,:), allocatable              :: atmIndx3
        real(dp), dimension(:,:,:,:,:,:), allocatable       :: FC3rd

        integer, dimension(3)                               :: Natm3Indx
        integer                                             :: Nbasis, Natm3, Natm3max, &
                                                             & mu, atm
        ! =========================== Local variables =========================== !

        call r3rdFC(filename, atmNum3, atmIndx3, FC3rd)

        this%atmNum = atmNum3
        this%atmIndx = atmIndx3
        this%fC = FC3rd

        Nbasis = size( atmNum3 )
        Natm3max = maxval( atmNum3 )

        allocate( this%basisTyp( Natm3max, Nbasis ) )
        allocate( this%atmPosR( 3, Natm3max, Nbasis ) )

        this%basisTyp(:, :) = 0
        this%atmPosR(:, :, :) = 0.0_dp

        mu_loop: do mu = 1, Nbasis

            Natm3 = atmNum3( mu )

            atm_loop: do atm = 1, Natm3

                Natm3Indx = atmIndx3(1:3, atm, mu)
                this%atmPosR(:, atm, mu) = matmul( sys%latvec, dble(Natm3Indx) )

                this%basisTyp(atm, mu) = atmIndx3(4, atm, mu)

            end do atm_loop

        end do mu_loop

        deallocate(atmNum3, atmIndx3, FC3rd)

    end subroutine set_FC3


    subroutine set_Zero(this)

        implicit none

        class(FC3type)              :: this

        allct_chk: if ( allocated(this%atmNum) .and. allocated(this%atmIndx) &
                      & .and. allocated(this%fC) ) then
            
            this%fC = 0.0_dp

        else allct_chk

            write(*, 355)
            355 FORMAT('ERROR in setting FC3 to zero: not allocated')

        end if allct_chk

    end subroutine set_Zero

    subroutine FC3_assignment(this, FC3_in)

        implicit none

        class(FC3type), intent(out)        :: this
        type(FC3type), intent(in)          :: FC3_in

        allct_chk: if ( allocated(FC3_in%atmNum) .and. allocated(FC3_in%atmIndx) &
                      & .and. allocated(FC3_in%fC) ) then

            this%atmNum = FC3_in%atmNum
            this%atmIndx = FC3_in%atmIndx
            this%fC = FC3_in%fC

        else allct_chk

            write(*, 455)
            455 FORMAT('ERROR: Assignment to a non-allocated FC3type')

        end if allct_chk

    end subroutine FC3_assignment

    !subroutine FC2_assignment_c(this, FC2_in)

    !    implicit none

    !    class(FC2type), intent(out)                     :: this
    !    complex(dp), dimension(:,:,:,:,:), intent(in)   :: FC2_in

    !    allct_chk: if ( allocated(this%atmNum) .and. allocated(this%atmIndx) ) then

    !        this%fC = dble(FC2_in)

    !    else allct_chk

    !        write(*, 455)
    !        455 FORMAT('ERROR: Assignment to a non-allocated FC2')

    !    end if allct_chk

    !end subroutine FC2_assignment_c


    Function FC3_add(this, FC3in2) RESULT(FC3out)

        implicit none

        class(FC3type), intent(in)      :: this
        type(FC3type), intent(in)       :: FC3in2
        type(FC3type)                   :: FC3out

        allct_chk: if ( allocated(this%fC) .and. allocated(FC3in2%fC) .and. &
                      & allocated(this%atmNum) .and. allocated(FC3in2%atmNum) .and. &
                      & allocated(this%atmIndx) .and. allocated(FC3in2%atmIndx) ) then

            FC3out%fC = this%fC + FC3in2%fC
            FC3out%atmNum = this%atmNum
            FC3out%atmIndx = this%atmIndx

        else allct_chk

            write(*, 555)
            555 FORMAT('ERROR: One of the FC3type to be added is not already allocated')

        end if allct_chk

    end Function FC3_add

end module FC3_mod

