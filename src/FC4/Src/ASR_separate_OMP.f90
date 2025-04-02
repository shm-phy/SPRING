
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

subroutine ApplyASR(this, Cut_Off, Nbasis)

    implicit none

    class(FCtyp)                                :: this

    type(CutOff), intent(in)                    :: Cut_Off
    integer, intent(in)                         :: Nbasis

    ! ====================================== Local Variables ===================================== !

    real(dp), dimension(:), allocatable     :: asr_row

    integer                                 :: NumASR, Nrow

    integer                                 :: NatmMax, chunk

    integer                                 :: mu, Natm_mu, N2, N3, N4, alpha, beta, gama, delta, &
                                             & asr_count

    ! ====================================== Local Variables ===================================== !

    ! --------------------------------- Initialize ASR variables --------------------------------- !

    NatmMax = maxval( Cut_Off%numAtom )
    NumASR = Nbasis * ( (NatmMax ** (ordfc4-2)) * (3**ordfc4) )

    Nrow = this%ind_f

    allocate( this%asr_mat(Nrow, NumASR) )
    this%asr_mat = 0.0_dp

    this%asr_zero = 0

    this%N_ASR = NumASR
    ! --------------------------------- Initialize ASR variables --------------------------------- !

    ! -------------------------------- Iter Over N1, N2, N3, cart -------------------------------- !
    write(*, 500) NumASR

    allocate( asr_row(Nrow) )
    asr_row = 0.0_dp

    chunk = 1 * 81

    !$omp parallel do default(shared) &
                !$omp schedule(dynamic, chunk) &
                !$omp private(mu, N2, N3, alpha, beta, gama, delta, Natm_mu, N4, asr_count) &
                !$omp firstprivate(asr_row) &
                !$omp shared(Nbasis, NatmMax, chunk, Cut_Off) &
                !$omp collapse(7)

    mu_loop: do mu = 1, Nbasis

        N2_loop: do N2 = 1, NatmMax
            N3_loop: do N3 = 1, NatmMax
                alpha_loop: do alpha = 1, 3
                    beta_loop: do beta = 1, 3
                        gama_loop: do gama = 1, 3
                            delta_loop: do delta = 1, 3

                                Natm_mu = Cut_Off%numAtom(mu)

                                NatmBound: if ( ( N2 <= Natm_mu) .and. (N3 <= Natm_mu) ) then
                                    ! *** !
                                    asr_row = 0.0_dp
                                    ! *** !

                                    N4_loop: do N4 = 1, Natm_mu

                                        call this%FindLC_ASR(mu, N2, N3, N4, alpha, beta, gama, delta, asr_row)

                                    end do N4_loop

                                    !^! Not necessary !^!
                                    WHERE( dabs(asr_row) < EPS )  asr_row = 0.0_dp
                                    !^! Not necessary !^!

                                    asr_count = (mu - 1) * ( (NatmMax**2) * 81 ) + &
                                              & (N2 - 1) * (NatmMax * 81) + &
                                              & (N3 - 1) * 81 + &
                                              & (alpha - 1) * 27 + &
                                              & (beta - 1) * 9 + &
                                              & (gama - 1) * 3 + &
                                              &  delta

                                    this%asr_mat(:, asr_count) = asr_row

                                    !*Off*! NonAdvancingOut: if ( this%NonAdvancWrite ) then
                                    !*Off*!     write(6, 520, advance="no") char(13), asr_count, NumASR
                                    !*Off*!     FLUSH(6)
                                    !*Off*! 
                                    !*Off*! else NonAdvancingOut
                                    !*Off*!     write(*, 550) asr_count, NumASR
                                    !*Off*! 
                                    !*Off*! end if NonAdvancingOut

                                end if NatmBound

                            end do delta_loop
                        end do gama_loop
                    end do beta_loop
                end do alpha_loop
            end do N3_loop
        end do N2_loop

    end do mu_loop

    !$omp end parallel do
    ! -------------------------------- Iter Over N1, N2, N3, cart -------------------------------- !

    deallocate( asr_row )
    write(*, *)

    500 FORMAT("Applying Translational Acoustic Sum Rule. Total Number of ASR conditions: ", I6)
    !*Off*! 520 FORMAT(1a1, "ASR Condition: ", I6, "/ ", I6)
    !*Off*! 550 FORMAT("ASR Condition: ", I6, "/ ", I6)

end subroutine ApplyASR


subroutine FindLC_ASR(this, mu, N2, N3, N4, alpha, beta, gama, delta, asr_row)

    implicit none

    class(FCtyp)                                :: this

    integer, intent(in)                         :: mu, N2, N3, N4, &
                                                 & alpha, beta, gama, delta
    real(dp), dimension(:), intent(inout)       :: asr_row

    ! ====================================== Local Variables ===================================== !

    real(dp)                                    :: dep_var1
    integer                                     :: colno

    ! ====================================== Local Variables ===================================== !

    dep_var1 = this%F(mu)%FCmu(1, delta, gama, beta, alpha, N4, N3, N2)

    depVarChk: if ( (dabs(dep_var1+1.0_dp) < EPS_nzero) .or. &
                  & (dabs(dep_var1+7.0_dp) < EPS_nzero) ) then

        call this%PntPerPreReduced(mu, N2, N3, N4, alpha, beta, gama, delta, asr_row)

    else if ( dabs(dep_var1+5.0_dp) < EPS_nzero ) then depVarChk

        call this%DepIndx5_ASR(mu, N2, N3, N4, alpha, beta, gama, delta, 1.0_dp, asr_row)

    else if ( dabs(dep_var1-1.0_dp) < EPS_nzero ) then depVarChk

        call this%find_col_ind(mu, N2, N3, N4, alpha, beta, gama, delta, colno)
        asr_row(colno) = asr_row(colno) + 1.0_dp

    else depVarChk

        not8: if ( dabs(dep_var1-8.0_dp) > EPS_nzero ) then
            write(*, 325) dep_var1
            ERROR STOP

        end if not8

    end if depVarChk

    325 FORMAT("ERROR( in FindLC_ASR ): Illegal value of dep var index ", F6.3)

end subroutine FindLC_ASR


subroutine PntPerPreReduced(this, mu, N2, N3, N4, alpha, beta, gama, delta, asr_row)

    implicit none

    class(FCtyp)                                :: this

    integer, intent(in)                         :: mu, N2, N3, N4, alpha, beta, gama, delta
    real(dp), dimension(:), intent(inout)       :: asr_row

    ! ====================================== Local Variables ===================================== !

    real(dp), dimension(3**ordfc4)              :: coef
    real(dp)                                    :: dep_var2, cf
    integer, dimension(ordfc4)                  :: ind_Indx
    integer                                     :: a_ind, res1, b_ind, res2, c_ind, d_ind, colno
    integer                                     :: i, ip

    ! ====================================== Local Variables ===================================== !

    ind_Indx = int( this%F(mu)%FCmu(3:6, delta, gama, beta, alpha, N4, N3, N2) )

    coef = this%F(mu)%FCmu(7:, delta, gama, beta, alpha, N4, N3, N2)

    zro_count: if ( all(dabs(coef) < EPS) ) then

        !$omp atomic
            this%asr_zero = this%asr_zero + 1

    end if zro_count

    coef_loop: do i = 1, (3**ordfc4)

        cf = coef(i)
        nonZerocf: if ( dabs(cf) > EPS ) then

            ip = (i - 1)

            a_ind = (ip / 27) + 1       !**!
            res1 = mod( ip, 27 )
            b_ind = (res1 / 9) + 1      !**!
            res2 = mod( res1, 9 )
            c_ind = (res2 / 3) + 1      !**!
            d_ind = mod( res2, 3 ) + 1  !**!

            dep_var2 = this%F(ind_Indx(1))%FCmu(1, d_ind, c_ind, b_ind, a_ind, &
                                         & ind_Indx(4), ind_Indx(3), ind_Indx(2))

            depVar2Chk: if ( dabs(dep_var2-1.0_dp) < EPS ) then

                call this%find_col_ind(ind_Indx(1), ind_Indx(2), ind_Indx(3), ind_Indx(4), &
                                     & a_ind, b_ind, c_ind, d_ind, colno)

                asr_row(colno) = asr_row(colno) + (1.0_dp * cf)

            else depVar2Chk
                write(*, 350) dep_var2
                ERROR STOP

            end if depVar2Chk

        end if nonZerocf

    end do coef_loop

    350 FORMAT("ERROR ( in PntPerPreReduced ): dep var index can not be other than 1 for pre-reduced ", F6.3)

end subroutine PntPerPreReduced


subroutine DepIndx5_ASR(this, mu, N2, N3, N4, alpha, beta, gama, delta, cof, asr_row)

    implicit none

    class(FCtyp)                                :: this

    integer, intent(in)                         :: mu, N2, N3, N4, alpha, beta, gama, delta
    real(dp), intent(in)                        :: cof
    real(dp), dimension(:), intent(inout)       :: asr_row

    ! ====================================== Local Variables ===================================== !

    real(dp), dimension(3**ordfc4)              :: coef_dep
    real(dp)                                    :: cf, check_var1
    integer                                     :: a_ind, res1, b_ind, res2, c_ind, d_ind, colno
    integer                                     :: i, ip

    ! ====================================== Local Variables ===================================== !

    coef_dep = this%F(mu)%FCmu(7:, delta, gama, beta, alpha, N4, N3, N2)

    debug: if ( all( dabs(coef_dep) < EPS_nzero ) ) then
        write(*, 375)
        STOP

    end if debug

    coef_depLoop: do i = 1, 3**ordfc4

        cf = coef_dep(i)

        nonZeroChk: if ( dabs(cf) > EPS ) then

            ip = (i - 1)

            a_ind = (ip / 27) + 1       !**!
            res1 = mod( ip, 27 )
            b_ind = (res1 / 9) + 1      !**!
            res2 = mod( res1, 9 )
            c_ind = (res2 / 3) + 1      !**!
            d_ind = mod( res2, 3 ) + 1  !**!

            check_var1 = this%F(mu)%FCmu(1, d_ind, c_ind, b_ind, a_ind, N4, N3, N2)

            CheckVar1Chk: if ( dabs(check_var1-1.0_dp) > EPS_nzero ) then
                write(*, 390) check_var1
                ERROR STOP

            end if CheckVar1Chk

            call this%find_col_ind(mu, N2, N3, N4, a_ind, b_ind, c_ind, d_ind, colno)
            asr_row(colno) = asr_row(colno) + (cof*cf)

        end if nonZeroChk

    end do coef_depLoop

    375 FORMAT("ERROR( in DepIndx5_ASR ): All coefficients can not be zero for dep_var=-5")
    390 FORMAT("ERROR( in DepIndx5_ASR ): dep_var for independent FC must be 1 ", F6.3)

end subroutine DepIndx5_ASR

