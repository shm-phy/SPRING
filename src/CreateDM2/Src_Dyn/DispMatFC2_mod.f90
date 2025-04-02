
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


module DispMatFC2

    use kinds,          only : dp
    use constants,      only : EPS, FILLVAL, ORDER
    use IndxChng,       only : PeriodicMap
    use FC2_data_mod,   only : FC2_dat, RedInFCsp

    implicit none

#ifdef _USE_DAXPY
    EXTERNAL            :: daxpy
#endif
    private

    public      :: IterOver_N2_mp2, find_force

contains

    subroutine  IterOver_N2_mp2(mu, alpha, num_row, nstep, NumDef_sup, PointDef, &
                              & cell_atm_def, my_FC_indx, sup_dim, FC2, disp, flags, Mat)

        implicit none

        integer, intent(in)                                 :: mu, alpha, num_row, nstep
        integer, intent(in)                                 :: NumDef_sup !** Variable for point-defect calculation **!
        logical, intent(in)                                 :: PointDef !** Variable for point-defect calculation **!
        integer, dimension(4, NumDef_sup), intent(in)       :: cell_atm_def !** Variable for point-defect calculation **!
        integer, dimension(3), intent(in)                   :: sup_dim
        integer, dimension(:, :), intent(in)                :: my_FC_indx

        type(FC2_dat), intent(in)                           :: FC2
        real(dp), dimension(:,:,:,:,:,:), intent(in)        :: disp !Assumed-shape dummy Array

        integer, dimension(num_row), intent(inout)          :: flags !** Variable for point-defect calculation **!
        real(dp), dimension(:,:), intent(inout)             :: Mat !Assumed-shape dummy Array

        !.......................  Local variables     .........................!

        real(dp)                                :: dv
        real(dp), dimension(num_row)            :: dsp_col !Automatic Array
        integer                                 :: ii, my_FC_num, cnt
        integer, dimension(2)                   :: FCi

        !.......................  Local variables     .........................!

        my_FC_num = size(my_FC_indx, 2)
        cnt = 1

        FC_part_loop: do ii = 1, my_FC_num

            FCi = my_FC_indx(:, ii)

            dv = FC2%F(1, FCi(2), alpha, FCi(1), mu)

            dv8_chk: if ( dabs(dv - 8.0_dp) > EPS ) then
                call FindLC(Mat, disp, num_row, nstep, sup_dim, FC2, &
                          & mu, FCi(1), alpha, FCi(2), NumDef_sup, PointDef, &
                          & cell_atm_def, dv, flags, dsp_col)

            end if dv8_chk

            cnt = cnt + 1

        end do FC_part_loop

    end subroutine IterOver_N2_mp2


    subroutine FindLC(Mat, disp, num_row, nstep, sup_dim, FC2, &
                    & mu, N2, alpha, beta, NumDef_sup, PointDef, &
                    & cell_atm_def, dv, flags, dsp_col)

        implicit none

        real(dp), dimension(:, :), intent(inout)            :: Mat !Assumed-shape dummy Array
        real(dp), dimension(:,:,:,:,:,:), intent(in)        :: disp !Assumed-shape dummy Array
        type(FC2_dat), intent(in)                           :: FC2
        integer, intent(in)                                 :: num_row, nstep, mu, N2, &
                                                             & alpha, beta
        integer, intent(in)                                 :: NumDef_sup !** Variable for point-defect calculation **!
        logical, intent(in)                                 :: PointDef !** Variable for point-defect calculation **!
        integer, dimension(4, NumDef_sup), intent(in)       :: cell_atm_def !** Variable for point-defect calculation **!
        integer, dimension(3)                               :: sup_dim
        real(dp), intent(in)                                :: dv

        integer, dimension(num_row), intent(inout)          :: flags !** Variable for point-defect calculation **!
        real(dp), dimension(num_row), intent(inout)         :: dsp_col !Explicit-shape dummy Array

        !.......................  Local variables     .........................!

        integer                                     :: col_no

        !.......................  Local variables     .........................!

        call find_disp_coeff(FC2, sup_dim, nstep, num_row, mu, N2, beta, &
                           & NumDef_sup, PointDef, cell_atm_def, disp, flags, dsp_col)

        dv_chk: if ( (dabs(dv+1.0_dp) < EPS) .or. (dabs(dv+7.0_dp) < EPS) ) then

            call find_ind17(Mat, FC2, num_row, mu, N2, alpha, beta, dsp_col)

        else if ( dabs(dv+5.0_dp) < EPS ) then dv_chk

            call find_ind5(Mat, FC2, num_row, mu, N2, alpha, beta, dsp_col)

        else if ( dabs(dv+9.0_dp) < EPS ) then dv_chk

            call find_ind9(Mat, FC2, num_row, mu, N2, alpha, beta, 1.0_dp, dsp_col)

        else if ( dabs(dv-1.0_dp) < EPS ) then dv_chk

            col_no = FC2%PosIndx(beta, alpha, N2, mu)

            col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC2%ind_fc) ) then
                write(*, 45)
                45 FORMAT('Something going wrong, cannot be FILLVAL')
                ERROR STOP

            else col_no_chk

#ifdef _USE_DAXPY
                call daxpy(num_row, 1.0_dp, dsp_col, 1, Mat(:, col_no), 1)
#else
                Mat(:, col_no) = Mat(:, col_no) + dsp_col
#endif

            end if col_no_chk

        else dv_chk
            write(*, 23)
            write(*, 13) dv
            write(*, 23)
            23 FORMAT('*********************** Something going wrong **************************')
            13 FORMAT('                     Illegal value of dv:', E14.6)
            ERROR STOP

        end if dv_chk
            
    end subroutine FindLC


    subroutine find_ind17(Mat, FC2, num_row, mu, N2, &
                        & alpha, beta, dsp_col)

        implicit none
        integer, parameter                                  :: MAX_CF = 9

        real(dp), dimension(:, :), intent(inout)            :: Mat !Assumed-shape dummy Array
        type(FC2_dat), intent(in)                           :: FC2
        integer, intent(in)                                 :: num_row, mu, N2, &
                                                             & alpha, beta
        real(dp), dimension(num_row), intent(in)            :: dsp_col !Explicit-shape dummy Array

        !.......................  Local variables     .........................!

        real(dp)                                :: cf, dep_var3
        real(dp), dimension(MAX_CF)             :: coef_lc
        integer, dimension(2)                   :: fc_lc
        integer                                 :: jj, jj_p, col_no
        integer                                 :: a_lc, b_lc

        !.......................  Local variables     .........................!

        fc_lc = int( FC2%F(3:4, beta, alpha, N2, mu) )

        coef_lc = FC2%F(5:, beta, alpha, N2, mu)

        coef_lc_loop: do jj = 1, MAX_CF

            cf = coef_lc(jj)
            cf0_chk: if ( dabs(cf) > EPS ) then

                jj_p = jj - 1

                a_lc = (jj_p / 3) + 1       !**!
                b_lc = mod( jj_p, 3 ) + 1   !**!

                dep_var3 = FC2%F(1, b_lc, a_lc, fc_lc(2), fc_lc(1))

                dv3_chk: if ( dabs(dep_var3+9.0_dp) < EPS ) then

                    call find_ind9(Mat, FC2, num_row, fc_lc(1), fc_lc(2), a_lc, b_lc, cf, dsp_col)

                else if ( dabs(dep_var3-1.0_dp) < EPS ) then dv3_chk
                    col_no = FC2%PosIndx(b_lc, a_lc, fc_lc(2), fc_lc(1))

                    col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC2%ind_fc) ) then
                        write(*, 45)
                        45 FORMAT('Something going wrong, cannot be FILLVAL')
                        ERROR STOP

                    else col_no_chk

#ifdef _USE_DAXPY
                        call daxpy(num_row, cf, dsp_col, 1, Mat(:, col_no), 1)
#else
                        Mat(:, col_no) = Mat(:, col_no) + cf * dsp_col
#endif

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

    end subroutine find_ind17


    subroutine find_ind5(Mat, FC2, num_row, mu, N2, &
                       & alpha, beta, dsp_col)

        implicit none
        integer, parameter                                  :: MAX_CF = 9

        real(dp), dimension(:, :), intent(inout)            :: Mat !Assumed-shape dummy Array
        type(FC2_dat), intent(in)                           :: FC2
        integer, intent(in)                                 :: num_row, mu, N2, &
                                                             & alpha, beta
        real(dp), dimension(num_row), intent(in)            :: dsp_col !Explicit-shape dummy Array

        !.......................  Local variables     .........................!

        real(dp)                                :: cf, dep_var2

        real(dp), dimension(MAX_CF)             :: coef_lc
        integer                                 :: jj, jj_p, col_no
        integer                                 :: a_lc, b_lc

        !.......................  Local variables     .........................!

        coef_lc = FC2%F(5:, beta, alpha, N2, mu)

        coef_lc_loop: do jj = 1, MAX_CF

            cf = coef_lc(jj)
            cf0_chk: if ( dabs(cf) > EPS ) then

                jj_p = jj - 1

                a_lc = (jj_p / 3) + 1       !**!
                b_lc = mod( jj_p, 3 ) + 1   !**!

                dep_var2 = FC2%F(1, b_lc, a_lc, N2, mu)

                dv2_chk: if ( dabs(dep_var2+9.0_dp) < EPS ) then

                    call find_ind9(Mat, FC2, num_row, mu, N2, a_lc, b_lc, cf, dsp_col)

                else if ( dabs(dep_var2-1.0_dp) < EPS ) then dv2_chk
                    col_no = FC2%PosIndx(b_lc, a_lc, N2, mu)

                    col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC2%ind_fc) ) then
                        write(*, 45)
                        45 FORMAT('Something going wrong, cannot be FILLVAL')
                        ERROR STOP

                    else col_no_chk

#ifdef _USE_DAXPY
                        call daxpy(num_row, cf, dsp_col, 1, Mat(:, col_no), 1)
#else
                        Mat(:, col_no) = Mat(:, col_no) + cf * dsp_col
#endif

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

    end subroutine find_ind5


    subroutine find_ind9(Mat, FC2, num_row, mu, N2, &
                       & a, b, cf, dsp_col)

        implicit none
        real(dp), dimension(:, :), intent(inout)            :: Mat !Assumed-shape dummy Array
        type(FC2_dat), intent(in)                           :: FC2
        integer, intent(in)                                 :: num_row, mu, N2, &
                                                             & a, b
        real(dp), intent(in)                                :: cf
        real(dp), dimension(num_row), intent(in)            :: dsp_col !Explicit-shape dummy Array

        !.......................  Local variables     .........................!

        real(dp)                                :: mat_el, coeff
        integer                                 :: N, num_el, jj, N_ind, col_no
        integer, dimension(4)                   :: fc_ind

        integer                                 :: pivot, Indx_el
        !integer, dimension(:), allocatable      :: afterI
        !real(dp), dimension(:), allocatable     :: afterC

        !.......................  Local variables     .........................!

        N = FC2%PosIndx_old(b, a, N2, mu)
        N_chk: if ( N == FILLVAL ) then
            write(*, 45)
            45 FORMAT('Something going wrong, cannot be FILLVAL')
            ERROR STOP
        end if N_chk

        pivot = FC2%DepIndxMap(N)

        pvt_chk: if ( (pivot == -2) .or. (pivot == -8) .or. (pivot == FILLVAL) ) then
            write(*, 786) pivot
            786 FORMAT('Something going wrong in FC2%DepIndxMap', I6)
            ERROR STOP
        end if pvt_chk

        num_el = size( FC2%SprMat(pivot)%Indx )

        ! ** This copy is Unnecessary ** !
        !allocate( afterI(num_el) )
        !allocate( afterC(num_el) )
        !afterI = FC2%SprMat(pivot)%Indx
        !afterC = FC2%SprMat(pivot)%Coeff

        !call dcopy(num_el, FC2%SprMat(pivot)%Indx, 1, afterI, 1)
        !call dcopy(num_el, FC2%SprMat(pivot)%Coeff, 1, afterC, 1)
        ! ** This copy is Unnecessary ** !

        mat_el_loop: do jj = 1, num_el

            Indx_el = FC2%SprMat(pivot)%Indx(jj)
            mat_el = FC2%SprMat(pivot)%Coeff(jj)

            N_ind = N + Indx_el
            fc_ind = FC2%PosIndx_flat_old(:, N_ind)

            col_no = FC2%PosIndx(fc_ind(4), fc_ind(3), fc_ind(2), fc_ind(1))

            col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC2%ind_fc) ) then
                write(*, 45)
                ERROR STOP

            else col_no_chk

                coeff = cf * mat_el
#ifdef _USE_DAXPY
                call daxpy(num_row, coeff, dsp_col, 1, Mat(:, col_no), 1)
#else
                Mat(:, col_no) = Mat(:, col_no) + coeff * dsp_col
#endif
            end if col_no_chk
            
        end do mat_el_loop

        !deallocate( afterI )
        !deallocate( afterC )

    end subroutine find_ind9


    subroutine find_disp_coeff(FC2, sup_dim, nstep, num_row, mu, N2, beta, &
                             & NumDef_sup, PointDef, cell_atm_def, disp, flags, dsp_col)

        implicit none

        type(FC2_dat), intent(in)                               :: FC2
        integer, intent(in)                                     :: nstep, num_row, mu, N2, beta
        integer, intent(in)                                     :: NumDef_sup !** Variable for point-defect calculation **!
        logical, intent(in)                                     :: PointDef !** Variable for point-defect calculation **!
        integer, dimension(4, NumDef_sup), intent(in)           :: cell_atm_def !** Variable for point-defect calculation **!
        integer, dimension(3), intent(in)                       :: sup_dim
        real(dp), dimension(:,:,:,:,:,:), intent(in)            :: disp !Assumed-shape dummy Array

        integer, dimension(num_row), intent(inout)              :: flags !** Variable for point-defect calculation **!
        real(dp), dimension(num_row), intent(out)               :: dsp_col !Explicit-shape dummy Array

        !.......................  Local variables     .........................!

        !-! real(dp)                                    :: u2
        integer, dimension(4)                       :: N2crt_nu
        integer, dimension(3)                       :: N2crt, &
                                                     & cell, cell2
        integer                                     :: nu, numd, &
                                                     & cx, cy, cz, ts
        integer, dimension(4, ORDER)                :: cell_atm_ifc !** Variable for point-defect calculation **!
        integer                                     :: flag !** Variable for point-defect calculation **!

        !.......................  Local variables     .........................!

        N2crt_nu = FC2%cb_Indx(:, N2, mu)
        N2crt = N2crt_nu(1:3)
        nu = N2crt_nu(4)

        !_! allocate( dsp_col(num_row) )
        !_! dsp_col = 0.0_dp

        numd = 0
        !-! numd = 1
        flag = 0

        cx_loop: do cx = 1, sup_dim(1)
            cy_loop: do cy = 1, sup_dim(2)
                cz_loop: do cz = 1, sup_dim(3)

                    cell = (/cx, cy, cz/)
                    cell2 = PeriodicMap(sup_dim, cell, N2crt)

                    !** Point defect realted calculation **!
                    if ( PointDef ) then
                        cell_atm_ifc(1:3, 1) = cell
                        cell_atm_ifc(4, 1) = mu
                        cell_atm_ifc(1:3, 2) = cell2
                        cell_atm_ifc(4, 2) = nu
                    
                        flag = CheckDefCellAtom( NumDef_sup, cell_atm_def, cell_atm_ifc )
                    end if
                    !** Point defect realted calculation **!

                    ts_loop: do ts = 1, nstep
                        !-! u2 = disp(beta, nu, cell2(3), cell2(2), cell2(1), ts)

                        dsp_col(numd+ts) = ( disp(ts, beta, nu, cell2(3), cell2(2), cell2(1)) * (-1.0_dp) )
 
                        !-! dsp_col(numd) = ( u2 * (-1.0_dp) )
                        !-! numd = numd + 1

                        if ( PointDef ) flags(numd+ts) = flags(numd+ts) + flag

                    end do ts_loop

                    numd = numd + nstep

                end do cz_loop
            end do cy_loop
        end do cx_loop

        !** For Debug Purpose **!
        !chk_num: if ( numd /= num_row ) then
        !    write(*, 99) 
        !    99 FORMAT("Error: All the timesteps and cells are not covered")
        !    ERROR STOP
        !end if chk_num
        !** For Debug Purpose **!

    end subroutine find_disp_coeff


    subroutine find_force(sup_dim, nstep, num_row, mu, alpha, force, force_mua)

        implicit none

        integer, dimension(3), intent(in)                       :: sup_dim
        integer, intent(in)                                     :: nstep, num_row, mu, alpha
        real(dp), dimension(:,:,:,:,:,:), intent(in)            :: force !Assumed-shape dummy Array
        real(dp), dimension(:), allocatable, intent(out)        :: force_mua !Deferred-shape Array

        !.......................  Local variables     .........................!
        integer                                     :: numf, cx, cy, cz, ts
        !.......................  Local variables     .........................!

        allocate( force_mua(num_row) )
        force_mua = 0.0_dp

        numf = 1
        cx_loop: do cx = 1, sup_dim(1)
            cy_loop: do cy = 1, sup_dim(2)
                cz_loop: do cz = 1, sup_dim(3)

                    ts_loop: do ts = 1, nstep
                        force_mua(numf) = force(ts, alpha, mu, cz, cy, cx)
                        numf = numf + 1

                    end do ts_loop

                end do cz_loop
            end do cy_loop
        end do cx_loop

        chk_num: if ( (numf-1) /= num_row ) then
            write(*, 99) 
            99 FORMAT("Error: All the timesteps and cells are not covered")
            ERROR STOP
        end if chk_num

    end subroutine find_force


    Pure Function CheckDefCellAtom( NumDef_sup, cell_atm_def, cell_atm_ifc ) Result( flag )

        implicit none

        integer, intent(in)                                 :: NumDef_sup
        integer, dimension(4, NumDef_sup), intent(in)       :: cell_atm_def
        integer, dimension(4, ORDER), intent(in)            :: cell_atm_ifc

        integer                                             :: flag !Result

        !.......................  Local variables     .........................!

        integer, dimension(4)               :: cell_atm
        integer                             :: i_f, i_d

        !.......................  Local variables     .........................!

        flag = 0

        ifc_cell: do i_f = 1, ORDER

            cell_atm = cell_atm_ifc(:, i_f)

            def_cell: do i_d = 1, NumDef_sup

                if ( all((cell_atm - cell_atm_def(:, i_d)) == 0) ) flag = i_f

            end do def_cell

        end do ifc_cell

    end Function CheckDefCellAtom

end module DispMatFC2

