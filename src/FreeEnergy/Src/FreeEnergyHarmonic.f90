
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


module HarmonicFreeEnergy

    use kinds,          only : dp
    use constants,      only : hbarTHz, Boltzk, BEConvTHz
    
    use Irr_q_point,    only : q_points_data
    use phonon_m,       only : Phon

    implicit none
    private

    public                  :: FreeEnergy2nd

contains

    subroutine FreeEnergy2nd( Qpoints, Nbasis, my_Qsize, my_offset, &
                            & phonon_dat, T, FreeEngHarmonic )

        implicit none

        type(q_points_data), intent(in)                             :: Qpoints
        integer, intent(in)                                         :: Nbasis, my_Qsize, my_offset
        type(Phon), intent(in)                                      :: phonon_dat
        real(dp), intent(in)                                        :: T

        real(dp), intent(out)                                       :: FreeEngHarmonic

        !================================= Local variable ===================================!

        real(dp), dimension(:), allocatable                         :: omegaq
        real(dp)                                                    :: FreeEng

        integer                                                     :: ii, i0, Ndof, strt, s
        integer, dimension(3)                                       :: q

        logical                                                     :: q0chk

        !================================= Local variable ===================================!

        FreeEngHarmonic = 0.0_dp

        Ndof = 3 * Nbasis
        allocate( omegaq(Ndof) )

        q_loop: do ii = 1, my_Qsize

            i0 = my_offset + ii
            q = Qpoints%q_pnt_int(1:3, i0)
            q0chk = all( q == 0 )

            omegaq = phonon_dat%omega(:, i0)

            strt = 1
            if ( q0chk ) strt = 4

            s_loop: do s = strt, Ndof

                FreeEng = 0.5_dp * hbarTHz * omegaq(s) + & 
                        & Boltzk * T * dlog( 1.0_dp - dexp(-1.0_dp * BEConvTHz * omegaq(s) / T) )

                FreeEngHarmonic = FreeEngHarmonic + FreeEng

                !-Debug-! write(*, *) q, "|", omegaq(s), "|",  (-1.0_dp * BEConvTHz * omegaq(s) / T), "|", &
                !-Debug-!           & dlog( 1.0_dp - dexp(-1.0_dp * BEConvTHz * omegaq(s) / T) ), "|", FreeEng

            end do s_loop

        end do q_loop

        deallocate( omegaq )

    end subroutine FreeEnergy2nd

end module HarmonicFreeEnergy

