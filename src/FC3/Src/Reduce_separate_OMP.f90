
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

subroutine IterFC(this, Cut_Off, Nbasis)

    implicit none

    class(FCtyp)                                :: this

    type(CutOff), intent(in)                    :: Cut_Off
    integer, intent(in)                         :: Nbasis

    ! ================================ Local Variables ================================ !

    integer                             :: mu, N2, N3, alpha, beta, gama, Natm_mu
    integer                             :: chunk, NatmMax
    real(dp)                            :: dep_var1
    integer, dimension(ordfc3)          :: ind_Indx
    real(dp), dimension(3**ordfc3)      :: coef, coef_red

    ! ================================ Local Variables ================================ !

    write(*, 200)

    coef_red = 0.0_dp

    chunk = 2 * 27
    NatmMax = maxval( Cut_Off%numAtom )

    !$omp parallel do default(shared) &
                !$omp schedule(dynamic, chunk) &
                !$omp private(mu, N2, N3, alpha, beta, gama, Natm_mu, ind_Indx) &
                !$omp private(dep_var1, coef) &
                !$omp firstprivate(coef_red) &
                !$omp shared(Nbasis, NatmMax, chunk, Cut_Off) &
                !$omp collapse(6)

    mu_loop: do mu = 1, Nbasis

        N2_loop: do N2 = 1, NatmMax
            N3_loop: do N3 = 1, NatmMax
                alpha_loop: do alpha = 1, 3
                    beta_loop: do beta = 1, 3
                        gama_loop: do gama = 1, 3

                            Natm_mu = Cut_Off%numAtom(mu)

                            NatmBound: if ( (N2 <= Natm_mu) .and. (N3 <= Natm_mu) ) then

                                dep_var1 = this%F(mu)%FCmu(1, gama, beta, alpha, N3, N2)

                                depVarChk: if ( (dabs(dep_var1+1.0_dp) < EPS_nzero) .or. &
                                              & (dabs(dep_var1+7.0_dp) < EPS_nzero) ) then

                                    ind_Indx = int( this%F(mu)%FCmu(3:5, gama, beta, alpha, N3, N2) )
                                    coef = this%F(mu)%FCmu(6:, gama, beta, alpha, N3, N2)

                                    debug: if ( all( dabs(coef) < EPS ) ) then
                                        write(*, 225)
                                        ERROR STOP
                                    end if debug

                                    call this%find_ind2(ind_Indx(1), ind_Indx(2), ind_Indx(3), &
                                                      & coef, coef_red)

                                    this%F(mu)%FCmu(6:, gama, beta, alpha, N3, N2) = coef_red

                                else depVarChk

                                    not158: if ( .not. ( (dabs(dep_var1-1.0_dp) < EPS_nzero) .or. &
                                                       & (dabs(dep_var1+5.0_dp) < EPS_nzero) .or. &
                                                       & (dabs(dep_var1-8.0_dp) < EPS_nzero) ) ) then
                                        write(*, 265) dep_var1
                                        ERROR STOP

                                    end if not158

                                end if depVarChk

                            end if NatmBound

                        end do gama_loop
                    end do beta_loop
                end do alpha_loop
            end do N3_loop
        end do N2_loop
    end do mu_loop

    !$omp end parallel do

    200 FORMAT("Reducing all the dependent FC to linear combination of independent FC")
    225 FORMAT("ERROR( in IterFC ): All coefficients can not be zero for dep_var -1 or -7")
    265 FORMAT("ERROR( in IterFC ): Illegal value of dep_var", F6.3)

end subroutine IterFC


subroutine find_ind2(this, mu_ind, ind_N2, ind_N3, coef, coef_red)

    implicit none

    class(FCtyp)                                :: this

    integer, intent(in)                         :: mu_ind, ind_N2, ind_N3
    real(dp), dimension(3**ordfc3), intent(in)  :: coef

    real(dp), dimension(3**ordfc3), intent(out) :: coef_red

    ! ======================================== Local Variables ======================================== !

    real(dp)                                    :: ind_var, cf
    integer                                     :: a, res, b, c
    integer                                     :: i, ip

    ! ======================================== Local Variables ======================================== !

    coef_red = 0.0_dp

    coefLoop: do i = 1, (3**ordfc3)
        cf = coef(i)

        NonZero: if ( dabs(cf) > EPS ) then

            ip = (i - 1)

            a = (ip / 9) + 1     !**!
            res = mod( ip, 9 )
            b = (res / 3) + 1    !**!
            c = mod(res, 3) + 1  !**!

            ind_var = this%F(mu_ind)%FCmu(1, c, b, a, ind_N3, ind_N2)

            indVarChk: if ( dabs(ind_var-1.0_dp) < EPS_nzero ) then
                coef_red(i) = coef_red(i) + cf

            else if ( dabs(ind_var+5.0_dp) < EPS_nzero ) then indVarChk
                call this%find_ind5Red(mu_ind, ind_N2, ind_N3, a, b, c, cf, coef_red)

            else if ( dabs(ind_var-8.0_dp) < EPS_nzero ) then indVarChk
                coef_red(i) = coef_red(i) + (cf*0.0_dp)

            else indVarChk
                write(*, 285) ind_var
                ERROR STOP

            end if indVarChk

        end if NonZero

    end do coefLoop

    285 FORMAT("ERROR( in find_ind2 ): Illegal value of ind_var ", F6.3)

end subroutine find_ind2


subroutine find_ind5Red(this, mu_ind, ind_N2, ind_N3, &
                      & a, b, c, cof, coef_red)

    implicit none

    class(FCtyp)                                    :: this

    integer, intent(in)                             :: mu_ind, ind_N2, ind_N3, &
                                                     & a, b, c
    real(dp), intent(in)                            :: cof

    real(dp), dimension(3**ordfc3), intent(inout)   :: coef_red

    ! ================================== Local Variables ================================== !

    real(dp), dimension(3**ordfc3)              :: coef_dep
    real(dp)                                    :: cf, check_var1
    integer                                     :: a_ind, res, b_ind, c_ind
    integer                                     :: i, ip

    ! ================================== Local Variables ================================== !

    coef_dep = this%F(mu_ind)%FCmu(6:, c, b, a, ind_N3, ind_N2)

    debug: if ( all( dabs(coef_dep) < EPS ) ) then
        write(*, 300)
        ERROR STOP

    end if debug

    coef_depLoop: do i = 1, (3**ordfc3)

        cf = coef_dep(i)

        nonZeroChk: if ( dabs(cf) > EPS ) then
            ip = (i - 1)

            a_ind = (ip / 9) + 1     !**!
            res = mod( ip, 9 )
            b_ind = (res / 3) + 1    !**!
            c_ind = mod(res, 3) + 1  !**!

            check_var1 = this%F(mu_ind)%FCmu(1, c_ind, b_ind, a_ind, ind_N3, ind_N2)

            CheckVar1Chk: if ( dabs(check_var1-1.0_dp) > EPS_nzero ) then
                write(*, 320) check_var1
                ERROR STOP

            end if CheckVar1Chk

            coef_red(i) = coef_red(i) + (cof*cf)

        end if nonZeroChk

    end do coef_depLoop

    300 FORMAT("ERROR( in find_ind5 ): All coefficients can not be zero for dep_var=-5")
    320 FORMAT("ERROR( in find_ind5 ): dep_var for independent FC must be 1 ", F6.3)

end subroutine find_ind5Red

