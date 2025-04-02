
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


module phonon_m

    use kinds,          only : dp
    use DynaMat,        only : get_phonon_gv
    use unit_cell,      only : cell
    use FC2_mod,        only : FC2type
    use EwaldMod,       only : EwaldParam

    implicit none
    private

    type, public        :: Phon

        real(dp), dimension(:, :), allocatable          :: omega
        real(dp), dimension(:, :), allocatable          :: nBE
        complex(dp), dimension(:, :, :), allocatable    :: Evec
        real(dp), dimension(:, :, :), allocatable       :: grp_vel
        real(dp)                                        :: omega_min, omega_max

    contains
        procedure, public, pass                         :: set_Phon

    end type Phon
    
contains

    subroutine set_Phon(this, sys, FC2, qpoints, EwaldConst, T, LongEW)

        implicit none

        class(Phon)                                                 :: this
        type(cell), intent(in)                                      :: sys
        type(FC2type), intent(in)                                   :: FC2
        real(dp), dimension(:, :), intent(in)                       :: qpoints
        type(EwaldParam), intent(in)                                :: EwaldConst
        real(dp), intent(in)                                        :: T
        logical, intent(in)                                         :: LongEW

        ! ================================ Local variables ================================ !

        real(dp), dimension(:, :), allocatable         :: freq, nBE
        complex(dp), dimension(:, :, :), allocatable   :: Evec
        real(dp), dimension(:, :, :), allocatable      :: grp_vel
        real(dp)                                       :: freq_min, freq_max

        ! ================================ Local variables ================================ !

        call get_phonon_gv(sys, FC2, qpoints, EwaldConst, T, LongEw, &
                         & freq, nBE, Evec, grp_vel, freq_min, freq_max)

        this%omega = freq
        this%nBE = nBE
        this%Evec = Evec
        this%grp_vel = grp_vel
        this%omega_min = freq_min
        this%omega_max = freq_max

        deallocate( freq, nBE, Evec, grp_vel )

    end subroutine set_Phon

end module phonon_m

