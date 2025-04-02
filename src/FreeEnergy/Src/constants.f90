
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

    integer, parameter              :: MAX_ITER = 1000
    integer, parameter              :: print_step  = 50, &
                                       print_step2 = 5

    real(dp), parameter             :: PI = 4.0_dp * datan(1.0_dp)
    real(dp), parameter             :: cmp_prec = 1.0E-8_dp

    real(dp), parameter             :: zero_prec = 5.0E-5_dp
    real(dp), parameter             :: OneSixth = 0.16666666666666666_dp
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

    real(dp), parameter             :: hbarTHz = 6.582119569_dp * 1.0E-4_dp
    real(dp), parameter             :: Boltzk = 8.617333262_dp * 1.0E-5_dp
    real(dp), parameter             :: BEConvTHz = 6.582119569_dp * 10.0_dp / 8.617333262_dp

    real(dp), parameter             :: eVConv3_t1 = 1.602176634_dp * (6.02214076_dp ** 3) * &
                                                  & (1.054571817_dp ** 2) * 1000.0_dp / 16.0_dp

    real(dp), parameter             :: eVConv3_t2 = eVConv3_t1 * 2.0_dp

    real(dp), parameter             :: eVConv4 = (( 6.02214076_dp * 1.054571817_dp ) ** 2) * 0.25_dp

    complex(dp), parameter          :: iu = dcmplx(0.0_dp, 1.0_dp)

    private
    public                          :: print_step, print_step2, MaX_ITER, PI, cmp_prec, zero_prec, OneSixth, OneThird, EPS, &
                                     & eV_A2, wTHz, THzConv, vg_ms, A2Conv, BEConv, hbarTHz, Boltzk, BEConvTHz, &
                                     & eVConv3_t1, eVConv3_t2, eVConv4, iu

end module constants

