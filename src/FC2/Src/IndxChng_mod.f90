
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


module IndxChng

    implicit none
    private

    public  :: num2indx, indx2num, fold_indx, PeriodicMap

contains

    subroutine num2indx(dimen, num, indx)
    
        implicit none
    
        integer, contiguous, dimension(:), intent(in)               :: dimen
        integer, intent(in)                                         :: num
        integer, dimension(:), allocatable, intent(out)             :: indx
    
        integer                                                     :: nd, j, prd, num2
    
        nd = size(dimen)
        allocate(indx(nd))
    
        num2 = num - 1
        do j = 1, nd-1
            
            prd = product(dimen(j+1:))
    
            indx(j) = (num2 / prd)
            num2 = mod(num2, prd)
    
        end do
    
        indx(nd) = num2
        indx = indx + 1
    
    end subroutine num2indx
    
    Function indx2num(dimen, indx) Result(num)
    
        implicit none
    
        integer, contiguous, dimension(:), intent(in)               :: dimen
        integer, contiguous, dimension(:), intent(in)               :: indx
        integer                                                     :: num ! Result
    
        integer, dimension(:), allocatable                          :: indx2
        integer                                                     :: nd, j
    
        nd = size(dimen)
        allocate(indx2(nd))
    
        indx2 = indx - 1
        num = 0
        do j = 1, nd-1
            
            num = num + indx2(j) * product(dimen(j+1:))
    
        end do
    
        num = num + indx2(nd) + 1
    
    end Function indx2num
    
    Function fold_indx(dimen, indx) Result(fold)
    
        implicit none
    
        integer, dimension(3), intent(in)       :: dimen
        integer, dimension(3), intent(in)       :: indx
        integer, dimension(3)                   :: fold ! Result
    
        ! ======================= Local variables ======================= !
        integer, dimension(3)                   :: dimen_2
        integer                                 :: nn
        ! ======================= Local variables ======================= !
    
        dimen_2 = dimen / 2
        fold = indx
    
        do nn = 1, 3
    
            if (fold(nn) > dimen_2(nn)) then
                fold(nn) = fold(nn) - dimen(nn)
            end if
    
        end do
    
    end Function fold_indx
    
    Function PeriodicMap(SupBound, Center, Cell) Result(CellFold)
    
        implicit none
        integer, dimension(3), intent(in)       :: SupBound
        integer, dimension(3), intent(in)       :: Center
        integer, dimension(3), intent(in)       :: Cell
    
        integer, dimension(3)                   :: CellFold !Result
    
        integer, dimension(3)                   :: CellEff
        integer                                 :: ii
    
        CellEff = (Center-1) + Cell
    
        loop_el: do ii = 1, 3
            CellFold(ii) = modulo( CellEff(ii), SupBound(ii) ) + 1
        end do loop_el
    
    end Function PeriodicMap

end module IndxChng

