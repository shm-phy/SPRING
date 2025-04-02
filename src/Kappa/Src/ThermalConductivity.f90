
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

module kappa_phonon

    use kinds,              only : dp
    use constants,          only : kW_m_k
    use Helper,             only : kron_prod
    use IndxChng,           only : indx2num
    use TetrahedronQ,       only : q_tetrahedron
    use phonon_m,           only : Phon

    implicit none
    private
    public                      :: kappa_noIter, kappa_noIter_symm

contains

    subroutine kappa_noIter( Nq1points, Nq0points, Ndof, my_QsizeAllq, my_offsetAllq, &
                           & tetraQ, ph_q1, vol, T, LW_Allq0, Polarization3rd, Polarization4th, &
                           & SpectralResolve, ModeResolve, FourPhonon, Spectral_k, Mode_k, kappa )

        implicit none

        integer, intent(in)                                                     :: Nq1points, Nq0points, Ndof, &
                                                                                 & my_QsizeAllq, my_offsetAllq

        type(q_tetrahedron), intent(in)                                         :: tetraQ
        type(Phon), intent(in)                                                  :: ph_q1

        real(dp), intent(in)                                                    :: vol, T
        real(dp), dimension(Ndof, Nq0points), intent(in)                        :: LW_Allq0
        integer, dimension(2, 3), intent(in)                                    :: Polarization3rd
        integer, dimension(2, 4), intent(in)                                    :: Polarization4th

        logical, intent(in)                                                     :: SpectralResolve, ModeResolve, &
                                                                                 & FourPhonon

        real(dp), allocatable, dimension(:,:), codimension[:], intent(inout)    :: Spectral_k, Mode_k
        real(dp), dimension(9), intent(out)                                     :: kappa 

        ! ================================== Local Variables ================================== !

        character(len=128)                  :: msg

        real(dp)                            :: omega_qs, nBE_qs, Lw_qs, tau_qs, multiply_fact1, &
                                             & multiply_fact2
        real(dp), dimension(3)              :: grpVel_qs
        real(dp), dimension(9)              :: kronProd_vel, kappa_qs

        integer                             :: ii, qi, s, strt, eqvl_irrQpos, current_indx, istat, &
                                             & default_strt, default_end
        integer, dimension(3)               :: q_int
        integer, dimension(2)               :: bound, current_address

        ! ================================== Local Variables ================================== !

        ! ..............:::::::::::::: Coindexed Objects Allocation ::::::::::::::............. !
        Spectral: if ( SpectralResolve ) then

            allocate( Spectral_k(11, (Ndof*Nq0points))[*], STAT=istat, ERRMSG=msg )
            if ( istat /= 0 ) then
                write(*, 24) this_image(), msg
                24 FORMAT( " Memory allocation ERROR in image : ", I5, ". Error message: ", A128 )
                ERROR STOP
            end if
            Spectral_k = 0.0_dp
            SYNC ALL

            bound = (/Nq0points, Ndof/)

        end if Spectral

        ModeKappa: if ( ModeResolve ) then

            allocate( Mode_k(9, Ndof)[*], STAT=istat, ERRMSG=msg )
            if ( istat /= 0 ) then
                write(*, 22) this_image(), msg
                22 FORMAT( " Memory allocation ERROR in image : ", I5, ". Error message: ", A128 )
                ERROR STOP
            end if
            Mode_k = 0.0_dp
            SYNC ALL

        end if ModeKappa
        ! ..............:::::::::::::: Coindexed Objects Allocation ::::::::::::::............. !

        FindStrtEnd: if ( FourPhonon ) then
            default_strt = min( Polarization3rd(1, 1), Polarization4th(1, 1) ) 
            default_end = max( Polarization3rd(2, 1), Polarization4th(2, 1) )

        else FindStrtEnd
            default_strt = Polarization3rd(1, 1)
            default_end = Polarization3rd(2, 1)

        end if FindStrtEnd

        multiply_fact1 = kW_m_k / ( (T**2) * vol )
        multiply_fact2 = multiply_fact1 / dble(Nq1points)

        kappa(:) = 0.0_dp

        q_loop: do ii = 1, my_QsizeAllq

            qi = my_offsetAllq + ii
            q_int = tetraQ%q_pnt_int(1:3, qi)

            strt = default_strt
            if( all(q_int == 0) ) strt = 4

            eqvl_irrQpos = tetraQ%pos(qi)

            polarization_loop: do s = strt, default_end

                omega_qs = ph_q1%omega(s, qi)
                nBE_qs = ph_q1%nBE(s, qi)
                grpVel_qs = ph_q1%grp_vel(:, s, qi)

                Lw_qs = LW_Allq0(s, eqvl_irrQpos)
                tau_qs = 0.5_dp / Lw_qs

                kronProd_vel = kron_prod( grpVel_qs, grpVel_qs )

                kappa_qs = ( (omega_qs**2) * nBE_qs*(nBE_qs + 1.0_dp) * tau_qs ) * kronProd_vel

                kappa = kappa + kappa_qs

                if ( ModeResolve ) Mode_k(:, s) = Mode_k(:, s) + kappa_qs

                Spectral2: if ( SpectralResolve ) then

                    current_address = (/eqvl_irrQpos, s/)
                    current_indx = indx2num( bound, current_address )
                    Spectral_k(1, current_indx) = dble( s )
                    Spectral_k(2, current_indx) = omega_qs
                    Spectral_k(3:11, current_indx) = Spectral_k(3:11, current_indx) + (kappa_qs * multiply_fact1)

                end if Spectral2

            end do polarization_loop

        end do q_loop

        kappa = kappa * multiply_fact2
        if ( ModeResolve ) Mode_k = Mode_k * multiply_fact2

    end subroutine kappa_noIter


    subroutine kappa_noIter_symm( Nq1points, Nq0points, Ndof, my_Qsizeq0, my_offsetq0, &
                                & tetraQ, ph_q0, vol, T, LW_Allq0, KrnGrpVelSym, Polarization3rd, &
                                & Polarization4th, SpectralResolve, ModeResolve, FourPhonon, &
                                & Spectral_k, Mode_k, kappa )

        implicit none

        integer, intent(in)                                                     :: Nq1points, Nq0points, Ndof, &
                                                                                 & my_Qsizeq0, my_offsetq0

        type(q_tetrahedron), intent(in)                                         :: tetraQ
        type(Phon), intent(in)                                                  :: ph_q0

        real(dp), intent(in)                                                    :: vol, T
        real(dp), dimension(Ndof, Nq0points), intent(in)                        :: LW_Allq0
        real(dp), dimension(9, Ndof, Nq0points), intent(in)                     :: KrnGrpVelSym
        integer, dimension(2, 3), intent(in)                                    :: Polarization3rd
        integer, dimension(2, 4), intent(in)                                    :: Polarization4th

        logical, intent(in)                                                     :: SpectralResolve, ModeResolve, &
                                                                                 & FourPhonon

        real(dp), allocatable, dimension(:,:), codimension[:], intent(inout)    :: Spectral_k, Mode_k
        real(dp), dimension(9), intent(out)                                     :: kappa 

        ! ================================== Local Variables ================================== !

        character(len=128)                  :: msg

        real(dp)                            :: omega_qs, nBE_qs, Lw_qs, tau_qs, multiply_fact1, &
                                             & multiply_fact2
        real(dp), dimension(9)              :: kronProd_vel_sym, kappa_qs

        integer                             :: ii, qi, s, strt, current_indx, istat, &
                                             & default_strt, default_end
        integer, dimension(3)               :: q_int
        integer, dimension(2)               :: bound, current_address

        ! ================================== Local Variables ================================== !

        ! ..............:::::::::::::: Coindexed Objects Allocation ::::::::::::::............. !
        Spectral: if ( SpectralResolve ) then

            allocate( Spectral_k(11, (Ndof*Nq0points))[*], STAT=istat, ERRMSG=msg )
            if ( istat /= 0 ) then
                write(*, 24) this_image(), msg
                24 FORMAT( " Memory allocation ERROR in image : ", I5, ". Error message: ", A128 )
                ERROR STOP
            end if
            Spectral_k = 0.0_dp
            SYNC ALL

            bound = (/Nq0points, Ndof/)

        end if Spectral

        ModeKappa: if ( ModeResolve ) then

            allocate( Mode_k(9, Ndof)[*], STAT=istat, ERRMSG=msg )
            if ( istat /= 0 ) then
                write(*, 22) this_image(), msg
                22 FORMAT( " Memory allocation ERROR in image : ", I5, ". Error message: ", A128 )
                ERROR STOP
            end if
            Mode_k = 0.0_dp
            SYNC ALL

        end if ModeKappa
        ! ..............:::::::::::::: Coindexed Objects Allocation ::::::::::::::............. !

        FindStrtEnd: if ( FourPhonon ) then
            default_strt = min( Polarization3rd(1, 1), Polarization4th(1, 1) ) 
            default_end = max( Polarization3rd(2, 1), Polarization4th(2, 1) )

        else FindStrtEnd
            default_strt = Polarization3rd(1, 1)
            default_end = Polarization3rd(2, 1)

        end if FindStrtEnd

        multiply_fact1 = kW_m_k / ( (T**2) * vol )
        multiply_fact2 = multiply_fact1 / dble(Nq1points)

        kappa(:) = 0.0_dp

        q_loop: do ii = 1, my_Qsizeq0

            qi = my_offsetq0 + ii
            q_int = tetraQ%irr_q_pnt_int(1:3, qi)

            strt = default_strt
            if( all(q_int == 0) ) strt = 4

            polarization_loop: do s = strt, default_end

                omega_qs = ph_q0%omega(s, qi)
                nBE_qs = ph_q0%nBE(s, qi)

                Lw_qs = LW_Allq0(s, qi)
                tau_qs = 0.5_dp / Lw_qs

                kronProd_vel_sym = KrnGrpVelSym(:, s, qi)

                kappa_qs = ( (omega_qs**2) * nBE_qs*(nBE_qs + 1.0_dp) * tau_qs ) * kronProd_vel_sym

                kappa = kappa + kappa_qs

                if ( ModeResolve ) Mode_k(:, s) = Mode_k(:, s) + kappa_qs

                Spectral2: if ( SpectralResolve ) then

                    current_address = (/qi, s/)
                    current_indx = indx2num( bound, current_address )
                    Spectral_k(1, current_indx) = dble( s )
                    Spectral_k(2, current_indx) = omega_qs
                    Spectral_k(3:11, current_indx) = Spectral_k(3:11, current_indx) + (kappa_qs * multiply_fact1)

                end if Spectral2

            end do polarization_loop

        end do q_loop

        kappa = kappa * multiply_fact2
        if ( ModeResolve ) Mode_k = Mode_k * multiply_fact2

    end subroutine kappa_noIter_symm

end module kappa_phonon

