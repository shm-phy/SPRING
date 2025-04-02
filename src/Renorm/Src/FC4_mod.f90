
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


module FC4_mod

    use kinds,          only : dp
    use hdf5_wrap,      only : r4thFC

    implicit none
    private

    type, public    ::  FC4type
        
        integer, dimension(:), allocatable                  :: atmNum
        integer, dimension(:,:,:), allocatable              :: atmIndx
        real(dp), dimension(:,:,:,:,:,:,:,:), allocatable   :: fC

    contains
        procedure, public, pass         :: set_FC4
        procedure, private, nopass      :: get_index_num

    end type FC4type

contains

    subroutine set_FC4(this, filename, Rcb2)

        implicit none

        class(FC4type)                                      :: this
        character(len=*), intent(in)                        :: filename
        integer, dimension(:,:,:), intent(in)               :: Rcb2

        integer, dimension(:), allocatable                  :: atmNum4
        integer, dimension(:,:,:), allocatable              :: Indx4
        integer, dimension(:,:,:), allocatable              :: IndxNum4
        real(dp), dimension(:,:,:,:,:,:,:,:), allocatable   :: FC4th

        integer, dimension(4)                               :: cb_indx
        integer                                             :: Nbasis, Natm4, bb, ii, N

        call r4thFC(filename, atmNum4, Indx4, FC4th)

        Nbasis = size(Indx4, 3)
        Natm4 = size(Indx4, 2)

        allocate( IndxNum4(5, Natm4, Nbasis) )
        IndxNum4 = 0

        IndxNum4(1:4, :, :) = Indx4

        basis_loop: do bb = 1, Nbasis
            atm4_loop: do ii = 1, Natm4

                cb_indx = Indx4(:, ii, bb)
                N = get_index_num(Rcb2, bb, cb_indx)
                IndxNum4(5, ii, bb) = N

            end do atm4_loop
        end do basis_loop

        this%atmNum = atmNum4
        this%atmIndx = IndxNum4
        this%fC = FC4th

        deallocate(atmNum4, Indx4, IndxNum4, FC4th)

    end subroutine set_FC4

    Function get_index_num(Rcb2, mu, indx) Result(N)
    
        implicit none
    
        integer, dimension(:,:,:), intent(in)       :: RCb2
        integer, intent(in)                         :: mu
        integer, dimension(4), intent(in)           :: indx

        integer                                     :: N    !Result
    
        integer, dimension(4)                       :: tmp_arr
        integer                                     :: nn, num_col
    
        num_col = size(RCb2, 2)
    
        N = 0
        loop_atm: do nn = 1, num_col
    
            tmp_arr = Rcb2(:, nn, mu)
    
            chk: if ( all((tmp_arr-indx) .eq. 0) ) then
                    N = nn
                    exit loop_atm
            end if chk
    
        end do loop_atm
    
        if ( N == 0 ) then
            write(*, 55) 
            55 FORMAT("Error: Index not found")
        end if
    
    end Function get_index_num


end module FC4_mod

