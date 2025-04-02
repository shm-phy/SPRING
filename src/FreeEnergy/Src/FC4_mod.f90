
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
    use unit_cell,      only : cell
    use hdf5_wrap,      only : r4thFC

    implicit none
    private

    type, public    ::  FC4type
        
        integer, dimension(:), allocatable                  :: atmNum
        integer, dimension(:,:,:), allocatable              :: atmIndx
        integer, dimension(:,:), allocatable                :: basisTyp
        real(dp), dimension(:,:,:), allocatable             :: atmPosR
        real(dp), dimension(:,:,:,:,:,:,:,:), allocatable   :: fC

    contains
        procedure, public, pass         :: set_FC4

    end type FC4type

contains


    subroutine set_FC4(this, sys, filename)

        implicit none

        class(FC4type)                                      :: this

        type(cell), intent(in)                              :: sys
        character(len=*), intent(in)                        :: filename

        ! =============================== Local Variables =============================== !
        integer, dimension(:), allocatable                  :: atmNum4
        integer, dimension(:,:,:), allocatable              :: Indx4
        real(dp), dimension(:,:,:,:,:,:,:,:), allocatable   :: FC4th

        integer, dimension(3)                               :: Natm4Indx
        integer                                             :: Nbasis, Natm4, Natm4max, &
                                                             & mu, atm
        ! =============================== Local Variables =============================== !

        call r4thFC(filename, atmNum4, Indx4, FC4th)

        this%atmNum = atmNum4
        this%atmIndx = Indx4
        this%fC = FC4th

        Nbasis = size( atmNum4 )
        Natm4max = maxval( atmNum4 )

        allocate( this%basisTyp( Natm4max, Nbasis ) )
        allocate( this%atmPosR( 3, Natm4max, Nbasis ) )

        this%basisTyp(:, :) = 0
        this%atmPosR(:, :, :) = 0.0_dp

        mu_loop: do mu = 1, Nbasis

            Natm4 = atmNum4( mu )

            atm_loop: do atm = 1, Natm4

                Natm4Indx = Indx4(1:3, atm, mu)
                this%atmPosR(:, atm, mu) = matmul( sys%latvec, dble(Natm4Indx) )

                this%basisTyp(atm, mu) = Indx4(4, atm, mu)

            end do atm_loop

        end do mu_loop

        deallocate(atmNum4, Indx4, FC4th)

    end subroutine set_FC4

end module FC4_mod

