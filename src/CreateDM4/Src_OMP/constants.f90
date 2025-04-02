
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
    real(dp), parameter             :: zero_prec = 8.0E-4_dp
    real(dp), parameter             :: EPS = 1.0E-5                        
    real(dp), parameter             :: m1_6val = -0.16666666666666666_dp

    integer, parameter              :: MAX_ITER = 1000, ORDER = 4
    integer, parameter              :: FILLVAL = 1000000, MAX_DEF = 5


    real(dp), parameter             :: eV_A2 = ((1.602176634_dp * 2.99792458_dp)**2 &
                                             & * 6.241509074460763_dp * 0.1_dp)

    real(dp), parameter             :: THzConv = dsqrt(1.60217733_dp * 6.02214076_dp * 10.0_dp) * &
                                                       10.0_dp / (2.0_dp * PI)

    real(dp), parameter             :: A2Conv = dsqrt(1.60217733_dp * 6.02214076_dp * 10.0_dp) * &
                                              & 6.582119569_dp * 0.001_dp

    real(dp), parameter             :: BEConv = dsqrt(1.60217733_dp * 6.02214076_dp * 10.0_dp) * &
                                              & 6.582119569_dp * 100.0_dp / 8.617333262_dp

    complex(dp), parameter          :: iu = dcmplx(0.0_dp, 1.0_dp)

    private
    public                          :: PI, cmp_prec, zero_prec, EPS, MAX_ITER, ORDER, iu, &
                                     & eV_A2, THzConv, A2Conv, BEConv, FILLVAL, MAX_DEF, m1_6val

end module constants

