
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

    real(dp), parameter             :: PI = 4.0_dp * datan(1.0_dp)
    real(dp), parameter             :: cmp_prec = 1.0E-8_dp
    real(dp), parameter             :: eV_A2 = ((1.602176634_dp * 2.99792458_dp)**2 &
                                             & * 6.241509074460763_dp * 0.1_dp)

    real(dp), parameter             :: THzConv = dsqrt(1.60217733_dp * 6.02214076_dp * 10.0_dp) * &
                                                       10.0_dp / (2.0_dp * PI)
    real(dp), parameter             :: EPS = 1.0E-8                        

    complex(dp), parameter          :: iu = dcmplx(0.0_dp, 1.0_dp)

    private
    public                          :: PI, cmp_prec, eV_A2, THzConv, EPS, iu

end module constants

