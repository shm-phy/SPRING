
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

    integer                             :: mu, Natm, N2, alpha, beta
    real(dp)                            :: dep_var1
    integer, dimension(ordfc2)          :: ind_Indx
    real(dp), dimension(3**ordfc2)      :: coef, coef_red

    ! ================================ Local Variables ================================ !

    write(*, 200)

    coef_red = 0.0_dp

    mu_loop: do mu = 1, Nbasis
        Natm = Cut_Off%numAtom(mu)

        N2_loop: do N2 = 1, Natm
            alpha_loop: do alpha = 1, 3
                beta_loop: do beta = 1, 3

                    dep_var1 = this%F(mu)%FCmu(1, beta, alpha, N2)

                    depVarChk: if ( (dabs(dep_var1+1.0_dp) < EPS_nzero) .or. &
                                  & (dabs(dep_var1+7.0_dp) < EPS_nzero) ) then

                        ind_Indx = int( this%F(mu)%FCmu(3:4, beta, alpha, N2) )

                        coef = this%F(mu)%FCmu(5:, beta, alpha, N2)

                        debug: if ( all( dabs(coef) < EPS ) ) then
                            write(*, 225)
                            STOP
                        end if debug

                        call this%find_ind2(ind_Indx(1), ind_Indx(2), coef, coef_red)

                        this%F(mu)%FCmu(5:, beta, alpha, N2) = coef_red

                    else depVarChk

                        not158: if ( .not. ( (dabs(dep_var1-1.0_dp) < EPS_nzero) .or. &
                                           & (dabs(dep_var1+5.0_dp) < EPS_nzero) .or. &
                                           & (dabs(dep_var1-8.0_dp) < EPS_nzero) ) ) then
                            write(*, 265) dep_var1
                            STOP

                        end if not158

                    end if depVarChk

                end do beta_loop
            end do alpha_loop
        end do N2_loop
    end do mu_loop

    200 FORMAT("Reducing all the dependent FC to linear combination of independent FC")
    225 FORMAT("ERROR( in IterFC ): All coefficients can not be zero for dep_var -1 or -7")
    265 FORMAT("ERROR( in IterFC ): Illegal value of dep_var", F6.3)

end subroutine IterFC


subroutine find_ind2(this, mu_ind, ind_N2, coef, coef_red)

    implicit none

    class(FCtyp)                                :: this

    integer, intent(in)                         :: mu_ind, ind_N2
    real(dp), dimension(3**ordfc2), intent(in)  :: coef

    real(dp), dimension(3**ordfc2), intent(out) :: coef_red

    ! ======================================== Local Variables ======================================== !

    real(dp)                                    :: ind_var, cf
    integer                                     :: a, b
    integer                                     :: i, ip

    ! ======================================== Local Variables ======================================== !

    coef_red = 0.0_dp

    coefLoop: do i = 1, (3**ordfc2)
        cf = coef(i)

        NonZero: if ( dabs(cf) > EPS ) then

            ip = (i - 1)

            a = (ip / 3) + 1        !**!
            b = mod( ip, 3 ) + 1    !**!

            ind_var = this%F(mu_ind)%FCmu(1, b, a, ind_N2)

            indVarChk: if ( dabs(ind_var-1.0_dp) < EPS_nzero ) then
                coef_red(i) = coef_red(i) + cf

            else if ( dabs(ind_var+5.0_dp) < EPS_nzero ) then indVarChk
                call this%find_ind5Red(mu_ind, ind_N2, a, b, cf, coef_red)

            else if ( dabs(ind_var-8.0_dp) < EPS_nzero ) then indVarChk
                coef_red(i) = coef_red(i) + (cf*0.0_dp)

            else indVarChk
                write(*, 285) ind_var
                STOP

            end if indVarChk

        end if NonZero

    end do coefLoop

    285 FORMAT("ERROR( in find_ind2 ): Illegal value of ind_var ", F6.3)

end subroutine find_ind2


subroutine find_ind5Red(this, mu_ind, ind_N2, a, b, cof, coef_red)

    implicit none

    class(FCtyp)                                    :: this

    integer, intent(in)                             :: mu_ind, ind_N2, a, b

    real(dp), intent(in)                            :: cof

    real(dp), dimension(3**ordfc2), intent(inout)   :: coef_red

    ! ================================== Local Variables ================================== !

    real(dp), dimension(3**ordfc2)              :: coef_dep
    real(dp)                                    :: cf, check_var1
    integer                                     :: a_ind, b_ind
    integer                                     :: i, ip

    ! ================================== Local Variables ================================== !

    coef_dep = this%F(mu_ind)%FCmu(5:, b, a, ind_N2)

    debug: if ( all( dabs(coef_dep) < EPS ) ) then
        write(*, 300)
        STOP

    end if debug

    coef_depLoop: do i = 1, (3**ordfc2)

        cf = coef_dep(i)

        nonZeroChk: if ( dabs(cf) > EPS ) then
            ip = (i - 1)

            a_ind = (ip / 3) + 1        !**!
            b_ind = mod( ip, 3 ) + 1    !**!

            check_var1 = this%F(mu_ind)%FCmu(1, b_ind, a_ind, ind_N2)

            CheckVar1Chk: if ( dabs(check_var1-1.0_dp) > EPS_nzero ) then
                write(*, 320) check_var1
                STOP

            end if CheckVar1Chk

            coef_red(i) = coef_red(i) + (cof*cf)

        end if nonZeroChk

    end do coef_depLoop

    300 FORMAT("ERROR( in find_ind5 ): All coefficients can not be zero for dep_var=-5")
    320 FORMAT("ERROR( in find_ind5 ): dep_var for independent FC must be 1 ", F6.3)

end subroutine find_ind5Red

