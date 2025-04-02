
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
    use constants,      only : EPS, FILLVAL, m1_6val 
    use FC4_data_mod,   only : FC4_dat, RedInFCsp, Chunk7d
    use FCType,         only : FChunk7d
    use unit_cell,      only : cell

    use cudafor,        only : dim3, cudaGetLastError, cudaDeviceSynchronize, &
                             & cudaSuccess, cudaGetErrorString, cudaGetLastError, &
                             & cudaMemcpy, cudaMemcpyAsync

    use cublas_v2

    implicit none

    private

    public      :: IterOver_N2_N3_N4, find_force

contains

    subroutine  IterOver_N2_N3_N4(mu, alpha, num_row, nstep, &
                                & Nbasis, sup_dim, FC4, F, disp, Mat)

        implicit none

        integer, intent(in)                                 :: mu, alpha, num_row, &
                                                             & nstep, Nbasis
        integer, dimension(3), intent(in)                   :: sup_dim

        type(FC4_dat), intent(in)                           :: FC4
        type(FChunk7d), dimension(3,Nbasis), intent(in)     :: F

        real(dp), dimension(:,:,:,:,:,:), intent(in)        :: disp !Assumed-Shape Dummy Array

        real(dp), dimension(:,:), intent(inout), device     :: Mat !Assumed-Shape Dummy Array

        !.......................  Local variables     .........................!

        type(cublasHandle)                                      :: h

        real(dp), dimension(:,:,:,:,:,:), allocatable, device   :: disp_d
        real(dp), dimension(:), allocatable, device             :: dsp_col_d
        real(dp)                                                :: dv

        integer, dimension(:), allocatable, device              :: max_shape_d, N2crt_nu_d, N3crt_eta_d, N4crt_ro_d
        integer, dimension(4)                                   :: max_shape
        integer                                                 :: Natm
        integer                                                 :: N2, N3, N4, &
                                                                 & beta, gama, delta, &
                                                                 & istat, ierrSync, nElements
        !.......................  Local variables     .........................!

        istat = cublasCreate(h)
        if ( istat /= cudaSuccess ) write(*, *) "ERROR: cublasCreate(h) failed, istat = ", istat

        max_shape(1:3) = sup_dim
        max_shape(4) = nstep
        nElements = nstep * 3 * Nbasis * product( sup_dim )

        allocate( disp_d(nstep, 3, Nbasis, sup_dim(3), sup_dim(2), sup_dim(1)) )
        ierrSync = cudaGetLastError()
        if ( ierrSync /= cudaSuccess ) write(*, *) 'ERROR: ', cudaGetErrorString(ierrSync)

        istat = cudaMemcpy( disp_d, disp, nElements ) !**!
        if ( istat /= cudaSuccess ) write(*, *) "ERROR: cudaMemcpy( disp_d, disp, nElements ) failed, istat = ", istat

        allocate( dsp_col_d(num_row) )
        ierrSync = cudaGetLastError()
        if ( ierrSync /= cudaSuccess ) write(*, *) 'ERROR: ', cudaGetErrorString(ierrSync)

        allocate( max_shape_d(4) )
        ierrSync = cudaGetLastError()
        if ( ierrSync /= cudaSuccess ) write(*, *) 'ERROR: ', cudaGetErrorString(ierrSync)

        allocate( N2crt_nu_d(4) )
        ierrSync = cudaGetLastError()
        if ( ierrSync /= cudaSuccess ) write(*, *) 'ERROR: ', cudaGetErrorString(ierrSync)

        allocate( N3crt_eta_d(4) )
        ierrSync = cudaGetLastError()
        if ( ierrSync /= cudaSuccess ) write(*, *) 'ERROR: ', cudaGetErrorString(ierrSync)

        allocate( N4crt_ro_d(4) )
        ierrSync = cudaGetLastError()
        if ( ierrSync /= cudaSuccess ) write(*, *) 'ERROR: ', cudaGetErrorString(ierrSync)

        istat = cudaMemcpy( max_shape_d, max_shape, 4 ) !**!
        if ( istat /= cudaSuccess ) write(*, *) "ERROR: cudaMemcpy( max_shape_d, max_shape, 4 ) failed, istat = ", istat

        Natm = FC4%atmNum(mu)

        N2_loop: do N2 = 1, Natm
            N3_loop: do N3 = 1, Natm
                N4_loop: do N4 = 1, Natm
                    beta_loop: do beta = 1, 3
                        gama_loop: do gama = 1, 3
                            delta_loop: do delta = 1, 3

                                dv = F(alpha,mu)%FArr7d(1, delta, gama, beta, N4, N3, N2)

                                dv8_chk: if ( dabs(dv - 8.0_dp) > EPS ) then

                                    call FindLC(Mat, disp_d, dsp_col_d, Nbasis, num_row, nstep, &
                                              & max_shape_d, N2crt_nu_d, N3crt_eta_d, N4crt_ro_d, &
                                              & h, FC4, F, mu, N2, N3, N4, alpha, beta, gama, delta, dv)

                                end if dv8_chk

                            end do delta_loop
                        end do gama_loop
                    end do beta_loop
                end do N4_loop
            end do N3_loop
        end do N2_loop

        deallocate( dsp_col_d )
        deallocate( disp_d )
        deallocate( max_shape_d, N2crt_nu_d, N3crt_eta_d, N4crt_ro_d )

        ierrSync = cudaGetLastError()
        if ( ierrSync /= cudaSuccess ) write(*, *) 'ERROR: ', cudaGetErrorString(ierrSync)

        istat = cublasDestroy(h)
        if ( istat /= cudaSuccess ) write(*, *) "ERROR: cublasDestroy(h) failed, istat = ", istat

    end subroutine IterOver_N2_N3_N4


    subroutine FindLC(Mat, disp_d, dsp_col_d, Nbasis, num_row, nstep, &
                    & max_shape_d, N2crt_nu_d, N3crt_eta_d, N4crt_ro_d, &
                    & h, FC4, F, mu, N2, N3, N4, alpha, beta, gama, delta, dv)

        implicit none

        real(dp), dimension(:, :), intent(inout), device        :: Mat !Assumed-Shape Dummy Array

        real(dp), dimension(:,:,:,:,:,:), device, intent(in)    :: disp_d !Assumed-Shape Dummy Array

        real(dp), dimension(:), intent(inout), device           :: dsp_col_d

        integer, intent(in)                                     :: Nbasis, num_row, nstep, mu, N2, N3, N4, &
                                                                 & alpha, beta, gama, delta
        integer, dimension(4), device, intent(in)               :: max_shape_d
        integer, dimension(4), device, intent(inout)            :: N2crt_nu_d, N3crt_eta_d, N4crt_ro_d

        type(cublasHandle)                                      :: h
        type(FC4_dat), intent(in)                               :: FC4
        type(FChunk7d), dimension(3,Nbasis), intent(in)         :: F

        real(dp), intent(in)                                    :: dv

        !.......................  Local variables     .........................!

        integer                                     :: col_no, istat

        !.......................  Local variables     .........................!

        call find_disp_coeff(FC4, nstep, num_row, mu, N2, beta, &
                           & max_shape_d, N2crt_nu_d, N3crt_eta_d, N4crt_ro_d, &
                           & N3, gama, N4, delta, disp_d, dsp_col_d)

        dv_chk: if ( (dabs(dv+1.0_dp) < EPS) .or. (dabs(dv+7.0_dp) < EPS) ) then

            call find_ind17(Mat, h, FC4, num_row, mu, N2, N3, N4, &
                          & alpha, beta, gama, delta, Nbasis, F, dsp_col_d)

        else if ( dabs(dv+5.0_dp) < EPS ) then dv_chk

            call find_ind5(Mat, h, FC4, num_row, mu, N2, N3, N4, &
                         & alpha, beta, gama, delta, Nbasis, F, dsp_col_d)

        else if ( dabs(dv+9.0_dp) < EPS ) then dv_chk

            call find_ind9(Mat, h, FC4, num_row, mu, N2, N3, N4, &
                         & alpha, beta, gama, delta, 1.0_dp, dsp_col_d)

        else if ( dabs(dv-1.0_dp) < EPS ) then dv_chk

            col_no = FC4%PosIndx(mu)%Arr7d(delta, gama, beta, alpha, N4, N3, N2)

            col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC4%ind_fc) ) then
                write(*, 45) col_no
                45 FORMAT('Something going wrong, cannot be FILLVAL 1', I9)
                ERROR STOP

            else col_no_chk

                istat = cublasDaxpy_v2(h, num_row, 1.0_dp, dsp_col_d, 1, Mat(:, col_no), 1)
                if (istat /= cudaSuccess) write(*, *) "ERROR: cublasDaxpy_v2 failed at FindLC "

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


    subroutine find_ind17(Mat, h, FC4, num_row, mu, N2, N3, N4, &
                        & alpha, beta, gama, delta, Nbasis, F, dsp_col)

        implicit none
        integer, parameter                                  :: MAX_CF = 81

        real(dp), dimension(:, :), intent(inout), device    :: Mat !Assumed-Shape Dummy Array
        type(cublasHandle)                                  :: h

        type(FC4_dat), intent(in)                           :: FC4
        integer, intent(in)                                 :: num_row, mu, N2, N3, N4, &
                                                             & alpha, beta, gama, delta, Nbasis

        type(FChunk7d), dimension(3,Nbasis), intent(in)     :: F

        real(dp), dimension(num_row), intent(in), device    :: dsp_col !Explicit-Shape Dummy Array

        !.......................  Local variables     .........................!
        real(dp), dimension(MAX_CF)             :: coef_lc
        real(dp)                                :: cf, dep_var3
        integer, dimension(4)                   :: fc_lc
        integer                                 :: jj, jj_p, res1, res2, col_no
        integer                                 :: a_lc, b_lc, c_lc, d_lc, istat
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

                    call find_ind9(Mat, h, FC4, num_row, fc_lc(1), fc_lc(2), fc_lc(3), fc_lc(4), &
                                 & a_lc, b_lc, c_lc, d_lc, cf, dsp_col)

                else if ( dabs(dep_var3-1.0_dp) < EPS ) then dv3_chk

                    col_no = FC4%PosIndx(fc_lc(1))%Arr7d(d_lc, c_lc, b_lc, a_lc, fc_lc(4), fc_lc(3), fc_lc(2))

                    col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC4%ind_fc) ) then
                        write(*, 45) col_no
                        45 FORMAT('Something going wrong, cannot be FILLVAL 2', I9)
                        ERROR STOP

                    else col_no_chk

                        istat = cublasDaxpy_v2(h, num_row, cf, dsp_col, 1, Mat(:, col_no), 1)
                        if (istat /= cudaSuccess) write(*, *) "ERROR: cublasDaxpy_v2 failed at find_ind17 "

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


    subroutine find_ind5(Mat, h, FC4, num_row, mu, N2, N3, N4, &
                       & alpha, beta, gama, delta, Nbasis, F, dsp_col)

        implicit none
        integer, parameter                                  :: MAX_CF = 81

        real(dp), dimension(:, :), intent(inout), device    :: Mat !Assumed-Shape Dummy Array
        type(cublasHandle)                                  :: h

        type(FC4_dat), intent(in)                           :: FC4
        integer, intent(in)                                 :: num_row, mu, N2, N3, N4, &
                                                             & alpha, beta, gama, delta, Nbasis

        type(FChunk7d), dimension(3,Nbasis), intent(in)     :: F

        real(dp), dimension(num_row), intent(in), device    :: dsp_col !Explicit-Shape Dummy Array

        !.......................  Local variables   .........................!

        real(dp), dimension(MAX_CF)             :: coef_lc
        real(dp)                                :: cf, dep_var2
        integer                                 :: jj, jj_p, res1, res2, col_no
        integer                                 :: a_lc, b_lc, c_lc, d_lc, istat

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

                    call find_ind9(Mat, h, FC4, num_row, mu, N2, N3, N4, &
                                 & a_lc, b_lc, c_lc, d_lc, cf, dsp_col)

                else if ( dabs(dep_var2-1.0_dp) < EPS ) then dv2_chk

                    col_no = FC4%PosIndx(mu)%Arr7d(d_lc, c_lc, b_lc, a_lc, N4, N3, N2)

                    col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC4%ind_fc) ) then
                        write(*, 45) col_no
                        45 FORMAT('Something going wrong, cannot be FILLVAL 3', I9)
                        ERROR STOP

                    else col_no_chk

                        istat = cublasDaxpy_v2(h, num_row, cf, dsp_col, 1, Mat(:, col_no), 1)
                        if (istat /= cudaSuccess) write(*, *) "ERROR: cublasDaxpy_v2 failed at find_ind5 "

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


    subroutine find_ind9(Mat, h, FC4, num_row, mu, N2, N3, N4, &
                       & a, b, c, d, cf, dsp_col)

        implicit none
        real(dp), dimension(:, :), intent(inout), device    :: Mat !Assumed-Shape Dummy Array
        type(cublasHandle)                                  :: h

        type(FC4_dat), intent(in)                           :: FC4
        integer, intent(in)                                 :: num_row, mu, N2, N3, N4, &
                                                             & a, b, c, d
        real(dp), intent(in)                                :: cf
        real(dp), dimension(num_row), intent(in), device    :: dsp_col !Explicit-Shape Dummy Array

        !.......................  Local variables     .........................!

        real(dp)                                :: mat_el, coeff
        integer                                 :: N, num_el, jj, N_ind, col_no
        integer, dimension(8)                   :: fc_ind
        integer                                 :: Indx_el

        integer                                 :: pivot, istat

        !.......................  Local variables     .........................!

        N = FC4%PosIndx_old(mu)%Arr7d(d, c, b, a, N4, N3, N2)

        !?! write(*, *) N, FILLVAL
        N_chk: if ( N == FILLVAL ) then
            write(*, 45) N
            45 FORMAT('Something going wrong, cannot be FILLVAL 4', I6)
            ERROR STOP
        end if N_chk

        pivot = FC4%DepIndxMap(N)

        pvt_chk: if ( (pivot == -2) .or. (pivot == -8) .or. (pivot == FILLVAL) ) then
            write(*, 786) pivot
            786 FORMAT('Something going wrong in FC4%DepIndxMap', I6)
            ERROR STOP
        end if pvt_chk

        num_el = size( FC4%SprMat(pivot)%Indx ) 

        mat_el_loop: do jj = 1, num_el

            Indx_el = FC4%SprMat(pivot)%Indx(jj)
            mat_el = FC4%SprMat(pivot)%Coeff(jj)

            N_ind = N + Indx_el
            fc_ind = FC4%PosIndx_flat_old(:, N_ind)

            col_no = FC4%PosIndx( fc_ind(1) )%Arr7d(fc_ind(8), fc_ind(7), fc_ind(6), &
                                                  & fc_ind(5), fc_ind(4), fc_ind(3), fc_ind(2))

            !?! write(*, *) col_no
            col_no_chk: if ( (col_no == FILLVAL) .or. (col_no > FC4%ind_fc) ) then
                write(*, 46) col_no, FILLVAL, FC4%ind_fc
                46 FORMAT('Something going wrong, cannot be FILLVAL 5', I9, I9, I9)
                ERROR STOP

            else col_no_chk

                coeff = cf * mat_el

                istat = cublasDaxpy_v2(h, num_row, coeff, dsp_col, 1, Mat(:, col_no), 1)
                if (istat /= cudaSuccess) write(*, *) "ERROR: cublasDaxpy_v2 failed at find_ind9 "

            end if col_no_chk
                
        end do mat_el_loop

    end subroutine find_ind9


    subroutine find_disp_coeff(FC4, nstep, num_row, mu, N2, beta, &
                              & max_shape_d, N2crt_nu_d, N3crt_eta_d, N4crt_ro_d, &
                              & N3, gama, N4, delta, disp_d, dsp_col_d)


        implicit none

        type(FC4_dat), intent(in)                               :: FC4
        integer, intent(in)                                     :: nstep, num_row, mu, N2, beta, &
                                                                 & N3, gama, N4, delta
        integer, dimension(4), device, intent(in)               :: max_shape_d
        integer, dimension(4), device, intent(inout)            :: N2crt_nu_d, N3crt_eta_d, N4crt_ro_d

        real(dp), dimension(:,:,:,:,:,:), device, intent(in)    :: disp_d !Assumed-Shape Dummy Array
        real(dp), dimension(num_row), device, intent(inout)     :: dsp_col_d !Explicit-Shape Dummy Array

        !.......................  Local variables  .........................!

        type(dim3)                                  :: blocks, threads

        integer, dimension(4)                       :: N2crt_nu, N3crt_eta, N4crt_ro
        integer                                     :: tPB_x, Nblock_x
        integer                                     :: ierrSync, ierrAsync

        !.......................  Local variables  .........................!

        N2crt_nu = FC4%cb_Indx(:, N2, mu)
        N3crt_eta = FC4%cb_Indx(:, N3, mu)
        N4crt_ro = FC4%cb_Indx(:, N4, mu)

        ! ** Host to Device memory transfer ** !
        N2crt_nu_d = N2crt_nu
        N3crt_eta_d = N3crt_eta
        N4crt_ro_d = N4crt_ro
        ! ** Host to Device memory transfer ** !

        !tPB_x = nstep
        tPB_x = 2 * 32
        Nblock_x = CEILING( num_row / dble(tPB_x) )
        blocks = dim3( Nblock_x, 1, 1 )
        threads = dim3( tPB_x, 1, 1 )

        call dispColumnDevice <<<blocks, threads>>> ( max_shape_d, N2crt_nu_d, N3crt_eta_d, N4crt_ro_d, &
                                                    & beta, gama, delta, num_row, disp_d, dsp_col_d )

        ierrSync = cudaGetLastError()
        ierrAsync = cudaDeviceSynchronize()
        if ( ierrSync /= cudaSuccess ) write(*, *) 'Sync kernel ERROR: ', cudaGetErrorString(ierrSync)
        if ( ierrAsync /= cudaSuccess ) write(*, *) 'ASync kernel ERROR: ', &
                                                   & cudaGetErrorString( cudaGetLastError() )

    end subroutine find_disp_coeff


    attributes(global) subroutine dispColumnDevice( max_shape, N2crt_nu, N3crt_eta, N4crt_ro, &
                                                  & beta, gama, delta, nrow, disp, dsp_col )

        implicit none

        integer, dimension(4), intent(in)               :: max_shape
        integer, dimension(4), intent(in)               :: N2crt_nu, N3crt_eta, N4crt_ro
        integer, intent(in), value                      :: beta, gama, delta
        integer, intent(in), value                      :: nrow

        real(dp), dimension(:,:,:,:,:,:), intent(in)    :: disp !Assumed-Shape Dummy Array
        real(dp), dimension(nrow), intent(inout)        :: dsp_col !Explicit-Shape Dummy Array

        ! ======================================== Local Variables ========================================= !

        integer                                     :: ii, d, prd, my_thread_no
        integer, dimension(3)                       :: Center

        integer, dimension(3)                       :: N2crt, N3crt, N4crt !This can be shared memory
        integer                                     :: nu, eta, ro !This can be shared memory

        integer, dimension(3)                       :: CellEff, cell2, cell3, cell4
        integer                                     :: ts

        ! ======================================== Local Variables ========================================= !

        my_thread_no = ( (blockIdx%x - 1) * blockDim%x + threadIdx%x )

        check_bound: if (  my_thread_no <= nrow ) then

            ! ================================ Find basis and lattice index ================================ !
            N2crt = N2crt_nu(1:3)
            nu = N2crt_nu(4)

            N3crt = N3crt_eta(1:3)
            eta = N3crt_eta(4)

            N4crt = N4crt_ro(1:3)
            ro = N4crt_ro(4)
            ! ================================ Find basis and lattice index ================================ !

            ! ================ Find Lattice and timestep corresponding to the thread number ================ !
            ii = my_thread_no - 1
            do d = 1, 3

                prd = product( max_shape(d+1:) )

                Center(d) = (ii / prd)
                ii = mod(ii, prd)

            end do

            ts = ii + 1
            ! ================ Find Lattice and timestep corresponding to the thread number ================ !

            ! ======================================== Periodic map ======================================== !

            ! ------------------------- cell2 ------------------------- !
            CellEff = Center + N2crt

            do ii = 1, 3
                cell2(ii) = modulo( CellEff(ii), max_shape(ii) ) + 1
            end do
            ! ------------------------- cell2 ------------------------- !

            ! ------------------------- cell3 ------------------------- !
            CellEff = Center + N3crt

            do ii = 1, 3
                cell3(ii) = modulo( CellEff(ii), max_shape(ii) ) + 1
            end do
            ! ------------------------- cell3 ------------------------- !

            ! ------------------------- cell4 ------------------------- !
            CellEff = Center + N4crt

            do ii = 1, 3
                cell4(ii) = modulo( CellEff(ii), max_shape(ii) ) + 1
            end do
            ! ------------------------- cell4 ------------------------- !

            ! ======================================== Periodic map ======================================== !

            ! ================================ displacement column element ================================= !

            dsp_col(my_thread_no) = ( disp(ts, beta,  nu, cell2(3), cell2(2), cell2(1)) * &
                                    & disp(ts, gama, eta, cell3(3), cell3(2), cell3(1)) * &
                                    & disp(ts, delta, ro, cell4(3), cell4(2), cell4(1)) * m1_6val )

            ! ================================ displacement column element ================================= !

        end if check_bound

    end subroutine dispColumnDevice


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

end module DispMatFC4

