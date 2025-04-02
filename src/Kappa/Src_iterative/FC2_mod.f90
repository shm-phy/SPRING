
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

module FC2_mod

    use kinds,          only : dp
    use hdf5_wrap,      only : r2ndFC

    implicit none
    private

    type, public    ::  FC2type
        
        integer, dimension(:), allocatable                  :: atmNum
        integer, dimension(:,:,:), allocatable              :: atmIndx
        real(dp), dimension(:,:,:,:,:), allocatable         :: fC

    contains
        procedure, public, pass         :: set_FC2, set_Zero

        GENERIC     ::  ASSIGNMENT(=)   =>  FC2_assignment
        GENERIC     ::  OPERATOR(+)     =>  FC2_add
        procedure, private, pass        :: FC2_assignment, FC2_add

    end type FC2type

contains

    subroutine set_FC2(this, filename)

        implicit none

        class(FC2type)                                      :: this
        character(len=*), intent(in)                        :: filename
        integer, dimension(:), allocatable                  :: atmNum2
        integer, dimension(:,:,:), allocatable              :: atmIndx2
        real(dp), dimension(:,:,:,:,:), allocatable         :: FC2nd

        call r2ndFC(filename, atmNum2, atmIndx2, FC2nd)

        this%atmNum = atmNum2
        this%atmIndx = atmIndx2
        this%fC = FC2nd

        deallocate(atmNum2, atmIndx2, FC2nd)

    end subroutine set_FC2

    subroutine set_Zero(this)

        implicit none

        class(FC2type)              :: this

        allct_chk: if ( allocated(this%atmNum) .and. allocated(this%atmIndx) &
                      & .and. allocated(this%fC) ) then
            
            this%fC = 0.0_dp

        else allct_chk

            write(*, 355)
            355 FORMAT('ERROR in setting FC2 to zero: not allocated')

        end if allct_chk

    end subroutine set_Zero

    subroutine FC2_assignment(this, FC2_in)

        implicit none

        class(FC2type), intent(out)        :: this
        type(FC2type), intent(in)          :: FC2_in

        allct_chk: if ( allocated(FC2_in%atmNum) .and. allocated(FC2_in%atmIndx) &
                      & .and. allocated(FC2_in%fC) ) then

            this%atmNum = FC2_in%atmNum
            this%atmIndx = FC2_in%atmIndx
            this%fC = FC2_in%fC

        else allct_chk

            write(*, 455)
            455 FORMAT('ERROR: Assignment to a non-allocated FC2type')

        end if allct_chk

    end subroutine FC2_assignment

    Function FC2_add(this, FC2in2) RESULT(FC2out)

        implicit none

        class(FC2type), intent(in)      :: this
        type(FC2type), intent(in)       :: FC2in2
        type(FC2type)                   :: FC2out

        allct_chk: if ( allocated(this%fC) .and. allocated(FC2in2%fC) .and. &
                      & allocated(this%atmNum) .and. allocated(FC2in2%atmNum) .and. &
                      & allocated(this%atmIndx) .and. allocated(FC2in2%atmIndx) ) then

            FC2out%fC = this%fC + FC2in2%fC
            FC2out%atmNum = this%atmNum
            FC2out%atmIndx = this%atmIndx

        else allct_chk

            write(*, 555)
            555 FORMAT('ERROR: One of the FC2type to be added is not already allocated')

        end if allct_chk

    end Function FC2_add

end module FC2_mod

