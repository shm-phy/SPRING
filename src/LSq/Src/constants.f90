
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


module constants

    use kinds,      only : dp
    implicit none

    real(dp), parameter                     :: PI = 4.0_dp * datan(1.0_dp)

    real(dp), parameter                     :: cmp_prec = 1.0E-8_dp
    real(dp), parameter                     :: zero_prec = 8.0E-4_dp
    real(dp), parameter                     :: EPS = 1.0E-5_dp                        
    real(dp), parameter                     :: EPS_nzero = 1.0E-6_dp

    integer, parameter                      :: ordfc2 = 2
    integer, parameter                      :: FILLVAL = 1000000
    integer, parameter                      :: initMatsz2 = 2022

    real(dp), dimension(3, 3), parameter    :: Iden3 = reshape( (/1.0_dp, 0.0_dp, 0.0_dp, &
                                                                & 0.0_dp, 1.0_dp, 0.0_dp, &
                                                                & 0.0_dp, 0.0_dp, 1.0_dp/), (/3,3/) )
    integer, dimension(ordfc2), parameter   :: dummy_Comb = (/1, 2/)


    private
    public                          :: PI, cmp_prec, zero_prec, EPS, EPS_nzero, &
                                     & ordfc2, FILLVAL, initMatsz2, Iden3, dummy_Comb

end module constants

