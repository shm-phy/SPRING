
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


program main

    use kinds,                      only : dp
    use parse_cmd_line,             only : get_command_line, ShowWelcomeBanner
    use unit_cell,                  only : cell
    use FC4_data_mod,               only : FC4_dat, RedInFCsp
    use FCtype,                     only : FChunk7d, set_FContg
    use hdf5_wrap,                  only : r_force_disp, wFD4_part, wFD_JoinAll
    use DispMatFC4,                 only : IterOver_N2_N3_N4, find_force

    use cudafor,                    only : cudaGetLastError, cudaGetErrorString, &
                                         & cudaSuccess, cudaMemcpy, cudaMemcpyAsync, &
                                         & cudaGetDeviceCount, cudaSetDevice

    use omp_lib,                    only : omp_get_thread_num

    implicit none

    type(cell)                                              :: sys
    character(len=256)                                      :: filename
    real(dp)                                                :: T

    type(FC4_dat)                                           :: FC4
    character(len=128)                                      :: fc4file

    type(FChunk7d), dimension(:,:), allocatable             :: NewFC

    character(len=24)                                       :: Temp_char, mu_char, alpha_char
    character(len=8)                                        :: frmt
    character(len=128)                                      :: fdfile
    real(dp), dimension(:,:,:,:,:,:), allocatable           :: disp
    real(dp), dimension(:,:,:,:,:,:), allocatable           :: force

    integer                                                 :: mu, alpha
    real(dp), allocatable, dimension(:,:)                   :: Mat
    real(dp), allocatable, dimension(:,:), device           :: Mat_d
    real(dp), allocatable, dimension(:)                     :: force_mua

    integer                                                 :: NIndFC, Nbasis, &
                                                             & nstep, num_row, &
                                                             & ierrSync, istat, &
                                                             & nElements, nDev_OMPthread

    integer                                                 :: omp_chunk, my_tid

    integer, dimension(3)                                   :: sup_dim

    logical                                                 :: dev_Mat_allocated

    character(len=64)                                       :: outfile, file_full


    call ShowWelcomeBanner()

    istat = cudaGetDeviceCount( nDev_OMPthread )
    if (istat /= cudaSuccess) write(*, *) "ERROR: cudaGetDeviceCount( nDev_OMPthread ) failed, istat = ", istat

    call get_command_line(filename, T, mu, alpha)
    call sys%init_data(filename)

    fc4file = 'FC_4th_common_'//trim(sys%prefix)//'_F.h5'
    call FC4%set_FC4dat(fc4file)

    call set_FContg(FC4, fc4file, NewFC)

    frmt = '(F6.1)'
    write(Temp_char, frmt) T
    fdfile = 'disp_forc_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
    call r_force_disp(fdfile, disp, force)


    nstep = size(disp, 1)
    sup_dim = sys%sup_cell
    num_row = nstep*product(sup_dim)
    NIndFC = FC4%ind_fc
    Nbasis = sys%natm
    nElements = num_row * NIndFC

    allocate( Mat(num_row, NIndFC))
    Mat = 0.0_dp

    write(*, *)
    write(*, "('No. of CUDA-capable devices detected: ', I3)") nDev_OMPthread
    write(*, *)

    frmt = '(I3)'
    dev_Mat_allocated = .false.
    omp_chunk = 2

    !$omp parallel do num_threads(nDev_OMPthread) default(shared) &
                !$omp schedule(dynamic, omp_chunk) &
                !$omp private(mu, alpha, istat, my_tid, ierrSync, force_mua, Mat_d, outfile) &
                !$omp private(mu_char, alpha_char) &
                !$omp firstprivate(sup_dim, dev_Mat_allocated, Mat, disp) &
                !$omp shared(num_row, nstep, Nbasis, NIndFC, nElements, force, frmt, Temp_char, FC4, NewFC, sys) &
                !$omp collapse(2)

    mu_loop: do mu = 1, Nbasis
        alpha_loop: do alpha = 1, 3

            my_tid = omp_get_thread_num()

            istat = cudaSetDevice( my_tid )
            if (istat /= cudaSuccess) write(*, *) "ERROR: cudaSetDevice( my_tid ) failed, istat = ", istat

            AllocateOnce: if ( .not. dev_Mat_allocated ) then

                allocate( Mat_d(num_row, NIndFC) )
                ierrSync = cudaGetLastError()
                if ( ierrSync /= cudaSuccess ) write(*, *) 'ERROR: ', cudaGetErrorString(ierrSync)

                dev_Mat_allocated = .true.

            end if AllocateOnce

            write(*, 100) mu, alpha, my_tid
            100 FORMAT( "Starting calculation of 4th-order displacement matrix with mu = ", I3, " and alpha = ", I2, &
                      & , ", OpenMP-thread / CUDA-device id = ", I4, " ..." )

            Mat_d = 0.0_dp
            call IterOver_N2_N3_N4(mu, alpha, num_row, nstep, &
                                 & Nbasis, sup_dim, FC4, NewFC, disp, Mat_d)

            istat = cudaMemcpy( Mat, Mat_d, nElements ) !**!
            if ( istat /= cudaSuccess ) write(*, *) "ERROR: cudaMemcpy( Mat, Mat_d, nElements ) failed, istat = ", istat

            call find_force(sup_dim, nstep, num_row, mu, alpha, force, force_mua)

            write(mu_char, frmt) mu
            write(alpha_char, frmt) alpha
            outfile = 'FD_4th_mu'//trim(adjustl(adjustr(mu_char)))//'_a'//trim(adjustl(adjustr(alpha_char)))&
                     &//'_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

            call wFD4_part(Mat, force_mua, outfile, mu, alpha)

        end do alpha_loop
    end do mu_loop

    !$omp end parallel do

    write(*, *)

    file_full = 'FD_4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'
    call wFD_JoinAll(Nbasis, num_row, NIndFC, Temp_char, sys%prefix, file_full)

end program main

