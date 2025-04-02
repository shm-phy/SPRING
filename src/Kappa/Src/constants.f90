
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

    integer, parameter              :: MAX_ITER = 1000, MAX_DEF = 5
    integer, parameter              :: print_step  = 500

    real(dp), parameter             :: PI = 4.0_dp * datan(1.0_dp)

    real(dp), parameter             :: Onebysqrt2PI = 1.0_dp / dsqrt(2.0_dp * PI) 

    real(dp), parameter             :: cmp_prec = 1.0E-8_dp

    real(dp), parameter             :: zero_prec = 5.0E-5_dp
    real(dp), parameter             :: OnebySix = 0.16666666666666666_dp
    real(dp), parameter             :: OneThird = 0.33333333333333331_dp

    real(dp), parameter             :: EPS = 1.0E-8                        

    real(dp), parameter             :: eV_A2 = ((1.602176634_dp * 2.99792458_dp)**2 &
                                             & * 6.241509074460763_dp * 0.1_dp)

    real(dp), parameter             :: wTHz = dsqrt(1.60217733_dp * 6.02214076_dp * 10.0_dp) * 10.0_dp 

    real(dp), parameter             :: THzConv = wTHz / (2.0_dp * PI)

    real(dp), parameter             :: vg_ms = wTHz * 100.0_dp

    real(dp), parameter             :: A2Conv = dsqrt(1.60217733_dp * 6.02214076_dp * 10.0_dp) * &
                                              & 6.582119569_dp * 0.001_dp

    real(dp), parameter             :: BEConv = dsqrt(1.60217733_dp * 6.02214076_dp * 10.0_dp) * &
                                              & 6.582119569_dp * 100.0_dp / 8.617333262_dp

    real(dp), parameter             :: BEConvTHz = 6.582119569_dp * 10.0_dp / 8.617333262_dp

    real(dp), parameter             :: Sctr3THz2 = 1.054571817_dp * ( (1.60217733_dp)**2 ) * &
                                                 & ( (6.02214076_dp)**3 ) * 1000000.0_dp * PI * 0.0625_dp

    real(dp), parameter             :: Sctr4THz2 = 0.125_dp * PI * ((1.054571817_dp * 1.602176565_dp)**2) * &
                                                 & (6.022140857_dp ** 4) * 1000000.0_dp

    real(dp), parameter             :: kW_m_k = (1.054571817_dp ** 2) * 0.001_dp / 1.380649_dp

    private
    public                          :: print_step, MaX_ITER, MAX_DEF, PI, Onebysqrt2PI, &
                                     & cmp_prec, zero_prec, OnebySix, OneThird, EPS, &
                                     & eV_A2, wTHz, THzConv, vg_ms, A2Conv, BEConv, BEConvTHz, &
                                     & Sctr3THz2, Sctr4THz2, kW_m_k

end module constants

