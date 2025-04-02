
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

module ReconstFC3

    use kinds,          only : dp
    use constants,      only : EPS, FILLVAL
    use FC3_data_mod,   only : FC3_dat, RedInFCsp

    implicit none
    private

    public      :: IterOver_N2_N3

contains

    subroutine  IterOver_N2_N3(FC3, Nbasis, FC_full)

        implicit none

        type(FC3_dat), intent(in)                           :: FC3
        integer, intent(in)                                 :: Nbasis

        real(dp), dimension(:,:,:,:,:,:), intent(inout)     :: FC_full !Assumed-shape dummy Array

        !.......................  Local variables     .........................!

        real(dp)                                :: fc_val
        integer                                 :: Nmax, Natm_mu, total, cnt, chunk
        integer                                 :: mu, N2, N3, alpha, &
                                                 & beta, gama

        !.......................  Local variables     .........................!

        Nmax = maxval( FC3%atmNum )
        total = Nbasis * (Nmax**2) * (3**3)
        cnt = 1

        chunk = 1 * 27

        !$omp parallel do default(shared) &
                    !$omp schedule(dynamic, chunk) &
                    !$omp private(fc_val, mu, N2, N3, alpha, beta, gama, Natm_mu) &
                    !$omp shared(FC3, FC_full, Nbasis, Nmax, total, cnt, chunk) &
                    !$omp collapse(6)

        mu_loop: do mu = 1, Nbasis

            N2_loop: do N2 = 1, Nmax
                N3_loop: do N3 = 1, Nmax
                    alpha_loop: do alpha = 1, 3
                        beta_loop: do beta = 1, 3
                            gama_loop: do gama = 1, 3

                                Natm_mu = FC3%atmNum(mu)

                                NatmBound: if ( (N2 <= Natm_mu) .and. (N3 <= Natm_mu) ) then

                                    call FindLC(FC3, mu, N2, N3, alpha, beta, gama, fc_val)

                                    FC_full( gama, beta, alpha, N3, N2, mu ) = fc_val

                                    !Off! write(*, 32) cnt, total
                                    !Off! 32 FORMAT(I6, ' /', I6)

                                    !Off! cnt = cnt + 1

                                end if NatmBound

                            end do gama_loop
                        end do beta_loop
                    end do alpha_loop
                end do N3_loop
            end do N2_loop

        end do mu_loop

        !$omp end parallel do

    end subroutine IterOver_N2_N3


    subroutine FindLC(FC3, mu, N2, N3, alpha, beta, gama, fc)

        implicit none

        type(FC3_dat), intent(in)                           :: FC3
        integer, intent(in)                                 :: mu, N2, N3, &
                                                             & alpha, beta, gama

        real(dp), intent(out)                               :: fc
        !.......................  Local variables     .........................!

        real(dp)                                    :: dv
        integer                                     :: col_no

        !.......................  Local variables     .........................!

        dv = FC3%F(1, gama, beta, alpha, N3, N2, mu)

        fc = 0.0_dp
        dv_chk: if ( (dabs(dv+1.0_dp) < EPS) .or. (dabs(dv+7.0_dp) < EPS) ) then

            call find_ind17(FC3, mu, N2, N3, alpha, beta, gama, fc)

        else if ( dabs(dv+5.0_dp) < EPS ) then dv_chk

            call find_ind5(FC3, mu, N2, N3, alpha, beta, gama, fc)

        else if ( dabs(dv+9.0_dp) < EPS ) then dv_chk

            call find_ind9(FC3, mu, N2, N3, alpha, beta, gama, 1.0_dp, fc)

        else if ( dabs(dv-1.0_dp) < EPS ) then dv_chk

            col_no = FC3%PosIndx(gama, beta, alpha, N3, N2, mu)

            col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC3%ind_fc) ) then
                write(*, 45)
                45 FORMAT('Something going wrong, cannot be FILLVAL')
                ERROR STOP

            else col_no_chk
                fc = FC3%FC3_ind(col_no)

            end if col_no_chk

        else if ( dabs(dv-8.0_dp) < EPS ) then dv_chk
            fc = 0.0_dp


        else dv_chk
            write(*, 23)
            write(*, 13) dv
            write(*, 23)
            23 FORMAT('*********************** Something going wrong **************************')
            13 FORMAT('                     Illegal value of dv:', E14.6)
            ERROR STOP

        end if dv_chk
            
    end subroutine FindLC


    subroutine find_ind17(FC3, mu, N2, N3, &
                        & alpha, beta, gama, fc)

        implicit none
        integer, parameter                                  :: MAX_CF = 27

        type(FC3_dat), intent(in)                           :: FC3
        integer, intent(in)                                 :: mu, N2, N3, &
                                                             & alpha, beta, gama

        real(dp), intent(out)                               :: fc

        !.......................  Local variables     .........................!

        real(dp), dimension(MAX_CF)             :: coef_lc
        real(dp)                                :: cf, dep_var3
        real(dp)                                :: my_fc17, my_fc9
        integer, dimension(3)                   :: fc_lc
        integer                                 :: jj, jj_p, res1, col_no
        integer                                 :: a_lc, b_lc, c_lc

        !.......................  Local variables     .........................!

        fc_lc = int( FC3%F(3:5, gama, beta, alpha, N3, N2, mu) )

        coef_lc = FC3%F(6:, gama, beta, alpha, N3, N2, mu)

        my_fc17 = 0.0_dp
        coef_lc_loop: do jj = 1, MAX_CF

            cf = coef_lc(jj)
            cf0_chk: if ( dabs(cf) > EPS ) then

                jj_p = jj - 1

                a_lc = (jj_p / 9) + 1       !**!
                res1 = mod( jj_p, 9 )
                b_lc = (res1 / 3) + 1       !**!
                c_lc = mod( res1, 3 ) + 1   !**!

                dep_var3 = FC3%F(1, c_lc, b_lc, a_lc, fc_lc(3), fc_lc(2), fc_lc(1))

                dv3_chk: if ( dabs(dep_var3+9.0_dp) < EPS ) then

                    call find_ind9(FC3, fc_lc(1), fc_lc(2), fc_lc(3), &
                                 & a_lc, b_lc, c_lc, cf, my_fc9)

                    my_fc17 = my_fc17 + my_fc9

                else if ( dabs(dep_var3-1.0_dp) < EPS ) then dv3_chk
                    col_no = FC3%PosIndx(c_lc, b_lc, a_lc, fc_lc(3), fc_lc(2), fc_lc(1))

                    col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC3%ind_fc) ) then
                        write(*, 45)
                        45 FORMAT('Something going wrong, cannot be FILLVAL')
                        ERROR STOP

                    else col_no_chk
                        my_fc17 = my_fc17 + ( cf * FC3%FC3_ind(col_no) )

                    end if col_no_chk

                else if ( (dabs(dep_var3+1.0_dp) < EPS) .or. (dabs(dep_var3+7.0_dp) < EPS) .or. &
                          (dabs(dep_var3+5.0_dp) < EPS) .or. (dabs(dep_var3-8.0_dp) < EPS) )  then dv3_chk
                    write(*, 23)
                    write(*, 13) dep_var3
                    write(*, 23)
                    23 FORMAT('*********************** Something going wrong **************************')
                    13 FORMAT('                     Illegal value of dv:', E14.6)
                    ERROR STOP

                else dv3_chk
                    write(*, 62) dep_var3
                    62 FORMAT('Illelegal value of dep_var3 encountered', E14.6)
                    ERROR STOP

                end if dv3_chk

            end if cf0_chk

        end do coef_lc_loop

        fc = my_fc17

    end subroutine find_ind17


    subroutine find_ind5(FC3, mu, N2, N3, &
                       & alpha, beta, gama, fc)

        implicit none
        integer, parameter                                  :: MAX_CF = 27

        type(FC3_dat), intent(in)                           :: FC3
        integer, intent(in)                                 :: mu, N2, N3, &
                                                             & alpha, beta, gama

        real(dp), intent(out)                               :: fc

        !.......................  Local variables     .........................!

        real(dp), dimension(MAX_CF)             :: coef_lc
        real(dp)                                :: cf, dep_var2
        real(dp)                                :: my_fc5, my_fc9
        integer                                 :: jj, jj_p, res1, col_no
        integer                                 :: a_lc, b_lc, c_lc

        !.......................  Local variables     .........................!

        coef_lc = FC3%F(6:, gama, beta, alpha, N3, N2, mu)

        my_fc5 = 0.0_dp
        coef_lc_loop: do jj = 1, MAX_CF

            cf = coef_lc(jj)
            cf0_chk: if ( dabs(cf) > EPS ) then

                jj_p = jj - 1

                a_lc = (jj_p / 9) + 1       !**!
                res1 = mod( jj_p, 9 )
                b_lc = (res1 / 3) + 1       !**!
                c_lc = mod( res1, 3 ) + 1   !**!

                dep_var2 = FC3%F(1, c_lc, b_lc, a_lc, N3, N2, mu)

                dv2_chk: if ( dabs(dep_var2+9.0_dp) < EPS ) then

                    call find_ind9(FC3, mu, N2, N3, a_lc, b_lc, c_lc, cf, my_fc9)
                    my_fc5 = my_fc5 + my_fc9

                else if ( dabs(dep_var2-1.0_dp) < EPS ) then dv2_chk
                    col_no = FC3%PosIndx(c_lc, b_lc, a_lc, N3, N2, mu)

                    col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC3%ind_fc) ) then
                        write(*, 45)
                        45 FORMAT('Something going wrong, cannot be FILLVAL')
                        ERROR STOP

                    else col_no_chk
                        my_fc5 = my_fc5 + ( cf * FC3%FC3_ind(col_no) )

                    end if col_no_chk

                else if ( (dabs(dep_var2+1.0_dp) < EPS) .or. (dabs(dep_var2+7.0_dp) < EPS) .or. &
                        & (dabs(dep_var2+5.0_dp) < EPS) .or. (dabs(dep_var2-8.0_dp) < EPS) )  then dv2_chk
                    write(*, 23)
                    write(*, 13) dep_var2
                    write(*, 23)
                    23 FORMAT('*********************** Something going wrong **************************')
                    13 FORMAT('                     Illegal value of dv:', E14.6)
                    ERROR STOP

                else dv2_chk
                    write(*, 62) dep_var2
                    62 FORMAT('Illelegal value of dep_var2 encountered', E14.6)
                    ERROR STOP

                end if dv2_chk

            end if cf0_chk

        end do coef_lc_loop

        fc = my_fc5

    end subroutine find_ind5


    subroutine find_ind9(FC3, mu, N2, N3, &
                       & a, b, c, cf_in, fc)

        implicit none
        type(FC3_dat), intent(in)                           :: FC3
        integer, intent(in)                                 :: mu, N2, N3, &
                                                             & a, b, c
        real(dp), intent(in)                                :: cf_in

        real(dp), intent(out)                               :: fc

        !.......................  Local variables     .........................!

        real(dp)                                :: mat_el
        real(dp)                                :: my_fc9
        integer                                 :: N, num_el, jj, N_ind, col_no
        integer, dimension(6)                   :: fc_ind

        integer                                 :: pivot, Indx_el
        integer, dimension(:), allocatable      :: afterI
        real(dp), dimension(:), allocatable     :: afterC

        !.......................  Local variables     .........................!

        N = FC3%PosIndx_old(c, b, a, N3, N2, mu)
        N_chk: if ( N == FILLVAL ) then
            write(*, 45)
            45 FORMAT('Something going wrong, cannot be FILLVAL')
            ERROR STOP
        end if N_chk

        pivot = FC3%DepIndxMap(N)

        pvt_chk: if ( (pivot == -2) .or. (pivot == -8) .or. (pivot == FILLVAL) ) then
            write(*, 786) pivot
            786 FORMAT('Something going wrong in FC3%DepIndxMap', I6)
            ERROR STOP
        end if pvt_chk

        num_el = size( FC3%SprMat(pivot)%Indx )
        allocate( afterI(num_el) )
        allocate( afterC(num_el) )

        afterI = FC3%SprMat(pivot)%Indx
        afterC = FC3%SprMat(pivot)%Coeff

        !call dcopy(num_el, FC3%SprMat(pivot)%Indx, 1, afterI, 1)
        !call dcopy(num_el, FC3%SprMat(pivot)%Coeff, 1, afterC, 1)

        my_fc9 = 0.0_dp
        mat_el_loop: do jj = 1, num_el

            Indx_el = afterI(jj)
            mat_el = afterC(jj)

            N_ind = N + Indx_el
            fc_ind = FC3%PosIndx_flat_old(:, N_ind)

            col_no = FC3%PosIndx(fc_ind(6), fc_ind(5), fc_ind(4), fc_ind(3), fc_ind(2), fc_ind(1))

            col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC3%ind_fc) ) then
                write(*, 45)
                ERROR STOP

            else col_no_chk
                my_fc9 = my_fc9 + ( cf_in * mat_el * FC3%FC3_ind(col_no) )

            end if col_no_chk
            
        end do mat_el_loop

        fc = my_fc9

        deallocate( afterI )
        deallocate( afterC )

    end subroutine find_ind9

end module ReconstFC3

