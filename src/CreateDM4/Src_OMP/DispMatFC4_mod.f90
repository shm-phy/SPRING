
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


module DispMatFC4

    use kinds,          only : dp
    use constants,      only : EPS, FILLVAL, m1_6val, ORDER 
    use IndxChng,       only : PeriodicMap
    use FC4_data_mod,   only : FC4_dat, RedInFCsp, Chunk7d
    use FCType,         only : FChunk7d
    use unit_cell,      only : cell

    use omp_lib,        only : omp_lock_kind, OMP_INIT_LOCK, &
                             & OMP_SET_LOCK, OMP_UNSET_LOCK
    implicit none
#ifdef _USE_DAXPY
    EXTERNAL                :: daxpy
#endif
    private

    public      :: IterOver_N2_N3_N4, find_force

contains

    subroutine  IterOver_N2_N3_N4( mu, alpha, num_row, nstep, NIndFC, Nbasis, NumDef_sup, PointDef, & 
                                 & cell_atm_def, sup_dim, FC4, F, disp, flags, Mat)

        implicit none

        integer, intent(in)                                 :: mu, alpha, num_row, &
                                                             & nstep, NIndFC, Nbasis
        integer, intent(in)                                 :: NumDef_sup !** Variable for point-defect calculation **!
        logical, intent(in)                                 :: PointDef !** Variable for point-defect calculation **!
        integer, dimension(4, NumDef_sup), intent(in)       :: cell_atm_def !** Variable for point-defect calculation **!
        integer, dimension(3), intent(in)                   :: sup_dim

        type(FC4_dat), intent(in)                           :: FC4
        type(FChunk7d), dimension(3,Nbasis), intent(in)     :: F

        real(dp), dimension(:,:,:,:,:,:), intent(in)        :: disp !Assumed-Shape Dummy Array

        integer, dimension(num_row), intent(inout)          :: flags !** Variable for point-defect calculation **!
        real(dp), dimension(:,:), intent(inout)             :: Mat !Assumed-Shape Dummy Array

        !.......................  Local variables     .........................!

        real(dp)                                    :: dv
        integer                                     :: Natm, chunk
        integer                                     :: N2, N3, N4, &
                                                     & beta, gama, delta, ll

        integer(omp_lock_kind), dimension(NIndFC)   :: ColLock !Automatic Array

        !.......................  Local variables     .........................!

        ! * Initialize ColLock * !
        lockVarLoop: do ll = 1, NIndFC

            call OMP_INIT_LOCK( ColLock(ll) )

        end do lockVarLoop

        chunk = 9 * 50

        Natm = FC4%atmNum(mu)

        !$omp parallel do default(shared) &
                    !$omp schedule(dynamic, chunk) &
                    !$omp private(N2, N3, N4, beta, gama, delta, dv) &
                    !$omp shared(Mat, flags, ColLock, disp, sup_dim, FC4, F) &
                    !$omp shared(mu, alpha, Nbasis, Natm, num_row, nstep, chunk) &
                    !$omp shared(NumDef_sup, PointDef, cell_atm_def) &
                    !$omp collapse(6)

        N2_loop: do N2 = 1, Natm
            N3_loop: do N3 = 1, Natm
                N4_loop: do N4 = 1, Natm
                    beta_loop: do beta = 1, 3
                        gama_loop: do gama = 1, 3
                            delta_loop: do delta = 1, 3

                                dv = F(alpha,mu)%FArr7d(1, delta, gama, beta, N4, N3, N2)

                                dv8_chk: if ( dabs(dv - 8.0_dp) > EPS ) then

                                    call FindLC( disp, Nbasis, num_row, nstep, sup_dim, FC4, F, mu, N2, N3, &
                                               & N4, alpha, beta, gama, delta, NumDef_sup, PointDef, cell_atm_def, &
                                               & dv, Mat, flags, ColLock )

                                end if dv8_chk

                            end do delta_loop
                        end do gama_loop
                    end do beta_loop
                end do N4_loop
            end do N3_loop
        end do N2_loop

        !$omp end parallel do

    end subroutine IterOver_N2_N3_N4


    subroutine FindLC( disp, Nbasis, num_row, nstep, sup_dim, FC4, F, mu, N2, N3, &
                     & N4, alpha, beta, gama, delta, NumDef_sup, PointDef, cell_atm_def, &
                     & dv, Mat, flags, ColLock )

        implicit none

        real(dp), dimension(:,:,:,:,:,:), intent(in)        :: disp !Assumed-Shape Dummy Array
        integer, intent(in)                                 :: Nbasis, num_row, nstep, mu, N2, N3, N4, &
                                                             & alpha, beta, gama, delta
        integer, intent(in)                                 :: NumDef_sup !** Variable for point-defect calculation **!
        logical, intent(in)                                 :: PointDef !** Variable for point-defect calculation **!
        integer, dimension(4, NumDef_sup), intent(in)       :: cell_atm_def !** Variable for point-defect calculation **!
        integer, dimension(3)                               :: sup_dim

        type(FC4_dat), intent(in)                           :: FC4
        type(FChunk7d), dimension(3,Nbasis), intent(in)     :: F

        real(dp), intent(in)                                :: dv

        real(dp), dimension(:, :), intent(inout)            :: Mat !Assumed-Shape Dummy Array
        integer, dimension(num_row), intent(inout)          :: flags !** Variable for point-defect calculation **!
        integer(omp_lock_kind), dimension(:), intent(inout) :: ColLock !Assumed-Shape Dummy Array

        !.......................  Local variables     .........................!

        real(dp), dimension(num_row)                :: dsp_col !Automatic Array
        integer                                     :: col_no

        !.......................  Local variables     .........................!

        call find_disp_coeff( FC4, sup_dim, nstep, num_row, mu, N2, beta, N3, gama, &
                            & N4, delta, NumDef_sup, PointDef, cell_atm_def, disp, flags, dsp_col )

        dv_chk: if ( (dabs(dv+1.0_dp) < EPS) .or. (dabs(dv+7.0_dp) < EPS) ) then

            call find_ind17(Mat, ColLock, FC4, num_row, mu, N2, N3, N4, &
                          & alpha, beta, gama, delta, Nbasis, F, dsp_col)

        else if ( dabs(dv+5.0_dp) < EPS ) then dv_chk

            call find_ind5(Mat, ColLock, FC4, num_row, mu, N2, N3, N4, &
                         & alpha, beta, gama, delta, Nbasis, F, dsp_col)

        else if ( dabs(dv+9.0_dp) < EPS ) then dv_chk

            call find_ind9(Mat, ColLock, FC4, num_row, mu, N2, N3, N4, &
                         & alpha, beta, gama, delta, 1.0_dp, dsp_col)

        else if ( dabs(dv-1.0_dp) < EPS ) then dv_chk

            col_no = FC4%PosIndx(mu)%Arr7d(delta, gama, beta, alpha, N4, N3, N2)

            col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC4%ind_fc) ) then
                write(*, 45)
                45 FORMAT('Something going wrong, cannot be FILLVAL')
                ERROR STOP

            else col_no_chk

                call OMP_SET_LOCK( ColLock(col_no) )
#ifdef _USE_DAXPY
                    call daxpy(num_row, 1.0_dp, dsp_col, 1, Mat(:, col_no), 1)
#else
                    Mat(:, col_no) = Mat(:, col_no) + dsp_col
#endif
                call OMP_UNSET_LOCK( ColLock(col_no) )

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


    subroutine find_ind17(Mat, ColLock, FC4, num_row, mu, N2, N3, N4, &
                        & alpha, beta, gama, delta, Nbasis, F, dsp_col)

        implicit none
        integer, parameter                                  :: MAX_CF = 81

        real(dp), dimension(:, :), intent(inout)            :: Mat !Assumed-Shape Dummy Array
        integer(omp_lock_kind), dimension(:), intent(inout) :: ColLock !Assumed-Shape Dummy Array

        type(FC4_dat), intent(in)                           :: FC4
        integer, intent(in)                                 :: num_row, mu, N2, N3, N4, &
                                                             & alpha, beta, gama, delta, Nbasis

        type(FChunk7d), dimension(3,Nbasis), intent(in)     :: F

        real(dp), dimension(num_row), intent(in)            :: dsp_col !Explicit-Shape Dummy Array

        !.......................  Local variables     .........................!
        real(dp), dimension(MAX_CF)             :: coef_lc
        real(dp)                                :: cf, dep_var3
        integer, dimension(4)                   :: fc_lc
        integer                                 :: jj, jj_p, res1, res2, col_no
        integer                                 :: a_lc, b_lc, c_lc, d_lc
        !.......................  Local variables     .........................!

        fc_lc = int( F(alpha,mu)%FArr7d(3:6, delta, gama, beta, N4, N3, N2) )

        coef_lc = F(alpha,mu)%FArr7d(7:, delta, gama, beta, N4, N3, N2)

        coef_lc_loop: do jj = 1, MAX_CF

            cf = coef_lc(jj)
            cf0_chk: if ( dabs(cf) > EPS ) then

                jj_p = jj - 1

                a_lc = (jj_p / 27) + 1      !**!
                res1 = mod( jj_p, 27 )
                b_lc = (res1 / 9) + 1       !**!
                res2 = mod( res1, 9 )
                c_lc = (res2 / 3) + 1       !**!
                d_lc = mod( res2, 3 ) + 1   !**!

                dep_var3 = F(a_lc,fc_lc(1))%FArr7d(1, d_lc, c_lc, b_lc, fc_lc(4), fc_lc(3), fc_lc(2))

                dv3_chk: if ( dabs(dep_var3+9.0_dp) < EPS ) then

                    call find_ind9(Mat, ColLock, FC4, num_row, fc_lc(1), fc_lc(2), fc_lc(3), fc_lc(4), &
                                 & a_lc, b_lc, c_lc, d_lc, cf, dsp_col)

                else if ( dabs(dep_var3-1.0_dp) < EPS ) then dv3_chk

                    col_no = FC4%PosIndx(fc_lc(1))%Arr7d(d_lc, c_lc, b_lc, a_lc, fc_lc(4), fc_lc(3), fc_lc(2))

                    col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC4%ind_fc) ) then
                        write(*, 45)
                        45 FORMAT('Something going wrong, cannot be FILLVAL')
                        ERROR STOP

                    else col_no_chk

                        call OMP_SET_LOCK( ColLock(col_no) )
#ifdef _USE_DAXPY
                            call daxpy(num_row, cf, dsp_col, 1, Mat(:, col_no), 1)
#else
                            Mat(:, col_no) = Mat(:, col_no) + cf * dsp_col
#endif
                        call OMP_UNSET_LOCK( ColLock(col_no) )

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


    subroutine find_ind5(Mat, ColLock, FC4, num_row, mu, N2, N3, N4, &
                       & alpha, beta, gama, delta, Nbasis, F, dsp_col)

        implicit none
        integer, parameter                                  :: MAX_CF = 81

        real(dp), dimension(:, :), intent(inout)            :: Mat !Assumed-Shape Dummy Array
        integer(omp_lock_kind), dimension(:), intent(inout) :: ColLock !Assumed-Shape Dummy Array

        type(FC4_dat), intent(in)                           :: FC4
        integer, intent(in)                                 :: num_row, mu, N2, N3, N4, &
                                                             & alpha, beta, gama, delta, Nbasis

        type(FChunk7d), dimension(3,Nbasis), intent(in)     :: F

        real(dp), dimension(num_row), intent(in)            :: dsp_col !Explicit-Shape Dummy Array

        !.......................  Local variables   .........................!

        real(dp), dimension(MAX_CF)             :: coef_lc
        real(dp)                                :: cf, dep_var2
        integer                                 :: jj, jj_p, res1, res2, col_no
        integer                                 :: a_lc, b_lc, c_lc, d_lc

        !.......................  Local variables   .........................!

        coef_lc = F(alpha,mu)%FArr7d(7:, delta, gama, beta, N4, N3, N2)

        coef_lc_loop: do jj = 1, MAX_CF

            cf = coef_lc(jj)
            cf0_chk: if ( dabs(cf) > EPS ) then

                jj_p = jj - 1

                a_lc = (jj_p / 27) + 1      !**!
                res1 = mod( jj_p, 27 )
                b_lc = (res1 / 9) + 1       !**!
                res2 = mod( res1, 9 )
                c_lc = (res2 / 3) + 1       !**!
                d_lc = mod( res2, 3 ) + 1   !**!

                dep_var2 = F(a_lc,mu)%FArr7d(1, d_lc, c_lc, b_lc, N4, N3, N2)

                dv2_chk: if ( dabs(dep_var2+9.0_dp) < EPS ) then

                    call find_ind9(Mat, ColLock, FC4, num_row, mu, N2, N3, N4, &
                                 & a_lc, b_lc, c_lc, d_lc, cf, dsp_col)

                else if ( dabs(dep_var2-1.0_dp) < EPS ) then dv2_chk

                    col_no = FC4%PosIndx(mu)%Arr7d(d_lc, c_lc, b_lc, a_lc, N4, N3, N2)

                    col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC4%ind_fc) ) then
                        write(*, 45)
                        45 FORMAT('Something going wrong, cannot be FILLVAL')
                        ERROR STOP

                    else col_no_chk

                        call OMP_SET_LOCK( ColLock(col_no) )
#ifdef _USE_DAXPY
                            call daxpy(num_row, cf, dsp_col, 1, Mat(:, col_no), 1)
#else
                            Mat(:, col_no) = Mat(:, col_no) + cf * dsp_col
#endif
                        call OMP_UNSET_LOCK( ColLock(col_no) )

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


    subroutine find_ind9(Mat, ColLock, FC4, num_row, mu, N2, N3, N4, &
                       & a, b, c, d, cf, dsp_col)

        implicit none
        real(dp), dimension(:, :), intent(inout)            :: Mat !Assumed-Shape Dummy Array
        integer(omp_lock_kind), dimension(:), intent(inout) :: ColLock !Assumed-Shape Dummy Array

        type(FC4_dat), intent(in)                           :: FC4
        integer, intent(in)                                 :: num_row, mu, N2, N3, N4, &
                                                             & a, b, c, d
        real(dp), intent(in)                                :: cf
        real(dp), dimension(num_row), intent(in)            :: dsp_col !Explicit-Shape Dummy Array

        !.......................  Local variables     .........................!

        real(dp)                                :: mat_el, coeff
        integer                                 :: N, num_el, jj, N_ind, col_no
        integer, dimension(8)                   :: fc_ind
        integer                                 :: Indx_el

        integer                                 :: pivot
        integer, dimension(:), allocatable      :: afterI !Deferred-Shape Array
        real(dp), dimension(:), allocatable     :: afterC !Deferred-Shape Array

        !.......................  Local variables     .........................!

        N = FC4%PosIndx_old(mu)%Arr7d(d, c, b, a, N4, N3, N2)

        N_chk: if ( N == FILLVAL ) then
            write(*, 45)
            45 FORMAT('Something going wrong, cannot be FILLVAL')
            ERROR STOP
        end if N_chk

        pivot = FC4%DepIndxMap(N)

        pvt_chk: if ( (pivot == -2) .or. (pivot == -8) .or. (pivot == FILLVAL) ) then
            write(*, 786) pivot
            786 FORMAT('Something going wrong in FC4%DepIndxMap', I6)
            ERROR STOP
        end if pvt_chk

        num_el = size( FC4%SprMat(pivot)%Indx ) 

        ! ** This copy is Unnecessary ** !
        allocate( afterI(num_el) )
        allocate( afterC(num_el) )
        afterI = FC4%SprMat(pivot)%Indx
        afterC = FC4%SprMat(pivot)%Coeff

        !call dcopy(num_el, FC4%SprMat(pivot)%Indx, 1, afterI, 1)
        !call dcopy(num_el, FC4%SprMat(pivot)%Coeff, 1, afterC, 1)
        ! ** This copy is Unnecessary ** !

        mat_el_loop: do jj = 1, num_el

            Indx_el = afterI(jj)
            mat_el = afterC(jj)

            N_ind = N + Indx_el
            fc_ind = FC4%PosIndx_flat_old(:, N_ind)

            col_no = FC4%PosIndx( fc_ind(1) )%Arr7d(fc_ind(8), fc_ind(7), fc_ind(6), &
                                                  & fc_ind(5), fc_ind(4), fc_ind(3), fc_ind(2))

            col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC4%ind_fc) ) then
                write(*, 45)
                ERROR STOP

            else col_no_chk

                coeff = cf * mat_el

                call OMP_SET_LOCK( ColLock(col_no) )
#ifdef _USE_DAXPY
                    call daxpy(num_row, coeff, dsp_col, 1, Mat(:, col_no), 1)
#else
                    Mat(:, col_no) = Mat(:, col_no) + coeff * dsp_col
#endif
                call OMP_UNSET_LOCK( ColLock(col_no) )

            end if col_no_chk
                
        end do mat_el_loop

        deallocate( afterI )
        deallocate( afterC )

    end subroutine find_ind9


    subroutine find_disp_coeff( FC4, sup_dim, nstep, num_row, mu, N2, beta, N3, gama, &
                              & N4, delta, NumDef_sup, PointDef, cell_atm_def, disp, flags, dsp_col )

        implicit none

        type(FC4_dat), intent(in)                               :: FC4
        integer, intent(in)                                     :: nstep, num_row, mu, N2, beta, &
                                                                 & N3, gama, N4, delta
        integer, intent(in)                                     :: NumDef_sup !** Variable for point-defect calculation **!
        logical, intent(in)                                     :: PointDef !** Variable for point-defect calculation **!
        integer, dimension(4, NumDef_sup), intent(in)           :: cell_atm_def !** Variable for point-defect calculation **!
        integer, dimension(3), intent(in)                       :: sup_dim
        real(dp), dimension(:,:,:,:,:,:), intent(in)            :: disp !Assumed-Shape Dummy Array

        integer, dimension(num_row), intent(inout)              :: flags !** Variable for point-defect calculation **!
        real(dp), dimension(num_row), intent(out)               :: dsp_col !Explicit-Shape Dummy Array

        !.......................  Local variables     .........................!

        !-! real(dp)                                    :: u2, u3, u4
        integer, dimension(4)                       :: N2crt_nu, N3crt_eta, N4crt_ro
        integer, dimension(3)                       :: N2crt, N3crt, N4crt, &
                                                     & cell, cell2, cell3, cell4
        integer                                     :: nu, eta, ro, numd, &
                                                     & cx, cy, cz, ts
        integer, dimension(4, ORDER)                :: cell_atm_ifc !** Variable for point-defect calculation **!
        integer                                     :: flag !** Variable for point-defect calculation **!

        !.......................  Local variables     .........................!

        N2crt_nu = FC4%cb_Indx(:, N2, mu)
        N2crt = N2crt_nu(1:3)
        nu = N2crt_nu(4)

        N3crt_eta = FC4%cb_Indx(:, N3, mu)
        N3crt = N3crt_eta(1:3)
        eta = N3crt_eta(4)

        N4crt_ro = FC4%cb_Indx(:, N4, mu)
        N4crt = N4crt_ro(1:3)
        ro = N4crt_ro(4)

        !_! dsp_col = 0.0_dp

        numd = 0
        !-! numd = 1
        flag = 0

        cx_loop: do cx = 1, sup_dim(1)
            cy_loop: do cy = 1, sup_dim(2)
                cz_loop: do cz = 1, sup_dim(3)

                    cell = (/cx, cy, cz/)
                    cell2 = PeriodicMap(sup_dim, cell, N2crt)
                    cell3 = PeriodicMap(sup_dim, cell, N3crt)
                    cell4 = PeriodicMap(sup_dim, cell, N4crt)

                    !** Point defect realted calculation **!
                    if ( PointDef ) then
                        cell_atm_ifc(1:3, 1) = cell
                        cell_atm_ifc(4, 1) = mu
                        cell_atm_ifc(1:3, 2) = cell2
                        cell_atm_ifc(4, 2) = nu
                        cell_atm_ifc(1:3, 3) = cell3
                        cell_atm_ifc(4, 3) = eta
                        cell_atm_ifc(1:3, 4) = cell4
                        cell_atm_ifc(4, 4) = ro

                        flag = CheckDefCellAtom( NumDef_sup, cell_atm_def, cell_atm_ifc )
                    end if
                    !** Point defect realted calculation **!

                    ts_loop: do ts = 1, nstep
                        !-! u2 = disp(beta, nu, cell2(3), cell2(2), cell2(1), ts)
                        !-! u3 = disp(gama, eta, cell3(3), cell3(2), cell3(1), ts)
                        !-! u4 = disp(delta, ro, cell4(3), cell4(2), cell4(1), ts)

                        dsp_col(numd+ts) = ( disp(ts, beta, nu, cell2(3), cell2(2), cell2(1)) * &
                                           & disp(ts, gama, eta, cell3(3), cell3(2), cell3(1)) * &
                                           & disp(ts, delta, ro, cell4(3), cell4(2), cell4(1)) * m1_6val )

                        !-! dsp_col(numd) = ( u2 * u3 * u4 * m1_6val )
                        !-! numd = numd + 1

                      if ( PointDef .and. (flag /= 0) ) then
                          !$omp atomic
                              flags(numd+ts) = flags(numd+ts) + flag
                      end if

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
        real(dp), dimension(:,:,:,:,:,:), intent(in)            :: force !Assumed-Shape Dummy Array
        real(dp), dimension(:), allocatable, intent(out)        :: force_mua !Deferred-Shape Array

        !.......................  Local variables     .........................!
        integer                                     :: numf, cx, cy, cz, ts
        !.......................  Local variables     .........................!

        allocate( force_mua(num_row) )
        !_! force_mua = 0.0_dp

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

end module DispMatFC4

