
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


module hdf5_wrap

    use hdf5,       only : hid_t, hsize_t, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER,&
                         & H5F_ACC_TRUNC_F, h5open_f, h5fcreate_f, &
                         & h5screate_simple_f, h5dcreate_f,  h5dwrite_f, &
                         & h5dclose_f, h5sclose_f, h5fclose_f, h5close_f, h5tclose_f, &
                         & HID_T, HSIZE_T, H5F_ACC_RDONLY_F, h5dget_type_f, &
                         & h5dget_space_f, h5sget_simple_extent_dims_f, h5dread_f, &
                         & h5dopen_f, h5fopen_f, h5gopen_f, h5gclose_f
    use kinds,      only : dp
    implicit none
    private
    public :: r2ndFC, r4thFC, r_force_disp, w2dh5_no_grp, wFC2_Renorm

contains

subroutine r2ndFC(filename, atmNum, atmIndx, FC2nd)

    implicit none

    character(len=*), intent(in)        :: filename
    integer, allocatable, intent(out)   :: atmNum(:)
    integer, allocatable, intent(out)   :: atmIndx(:, :, :)
    real(dp), allocatable, intent(out)  :: FC2nd(:, :, :, :, :)

    integer, parameter                  :: rank1 = 1, rank2 = 3, rank3 = 5
    character(len=50)                   :: dsetname1="NumAtoms", dsetname2="atm_Indx", dsetname3="FC2"

    integer(HID_T)                      :: file_id, dset_id1, type_id1, space_id1, &
                                                    dset_id2, type_id2, space_id2, &
                                                    dset_id3, type_id3, space_id3
    integer(HSIZE_T)                    :: dims1(rank1), dims2(rank2), dims3(rank3)
    integer(HSIZE_T)                    :: maxdims1(rank1), maxdims2(rank2), maxdims3(rank3)
    integer                             :: hdferr


    write(*, 155) filename
    155 FORMAT('Reading Second-Order Force constant file: ', A)

    call h5open_f(hdferr)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdferr) !output -> fileid

    call h5dopen_f(file_id, dsetname1, dset_id1, hdferr) !Opens an existing dataset.
                                                         !output -> dset_id1
    call h5dopen_f(file_id, dsetname2, dset_id2, hdferr) !Opens an existing dataset.
                                                         !output -> dset_id2
    call h5dopen_f(file_id, dsetname3, dset_id3, hdferr) !Opens an existing dataset.
                                                         !output -> dset_id3

    call h5dget_type_f(dset_id1, type_id1, hdferr) !output -> type_id1
    call h5dget_type_f(dset_id2, type_id2, hdferr) !output -> type_id2
    call h5dget_type_f(dset_id3, type_id3, hdferr) !output -> type_id3

    call h5dget_space_f(dset_id1, space_id1, hdferr) !output -> space_id1
    call h5sget_simple_extent_dims_f(space_id1, dims1, maxdims1, hdferr) !output -> dims1
                                                                                  ! maxdims1
    call h5dget_space_f(dset_id2, space_id2, hdferr) !output -> space_id2
    call h5sget_simple_extent_dims_f(space_id2, dims2, maxdims2, hdferr) !output -> dims2
                                                                                  ! maxdims2
    call h5dget_space_f(dset_id3, space_id3, hdferr) !output -> space_id3
    call h5sget_simple_extent_dims_f(space_id3, dims3, maxdims3, hdferr) !output -> dims3
                                                                                  ! maxdims3
    allocate(atmNum(1:dims1(1)))
    allocate(atmIndx(1:dims2(1), 1:dims2(2), 1:dims2(3)))
    allocate(FC2nd(1:dims3(1),1:dims3(2), 1:dims3(3), 1:dims3(4), 1:dims3(5)))

    call h5dread_f(dset_id1, type_id1, atmNum, dims1, hdferr) !Read in array 'atmNum'
    call h5dread_f(dset_id2, type_id2, atmIndx, dims2, hdferr) !Read in array 'atmIndx'
    call h5dread_f(dset_id3, type_id3, FC2nd, dims3, hdferr) !Read in array 'FC2'

    call h5sclose_f( space_id1, hdferr )
    call h5sclose_f( space_id2, hdferr )
    call h5sclose_f( space_id3, hdferr )
    
    call h5tclose_f( type_id1, hdferr )
    call h5tclose_f( type_id2, hdferr )
    call h5tclose_f( type_id3, hdferr )

    call h5dclose_f(dset_id1, hdferr)
    call h5dclose_f(dset_id2, hdferr)
    call h5dclose_f(dset_id3, hdferr)

    call h5fclose_f(file_id, hdferr)
    call h5close_f(hdferr)

end subroutine r2ndFC


subroutine r4thFC(filename, atmNum4, atmIndx4, FC_4th)

    implicit none

    character(len=*), intent(in)        :: filename
    integer, allocatable, intent(out)   :: atmNum4(:)
    integer, allocatable, intent(out)   :: atmIndx4(:, :, :)
    real(dp), allocatable, intent(out)  :: FC_4th(:,:,:,:,:,:,:,:)

    integer(HID_T), parameter                           :: rank1 = 1, rank2 = 3, rank3 = 7

    character(len=50)                                   :: fc_grp_name="FC4th", &
                                                         & dsetname1="NumAtoms", &
                                                         & dsetname2="atm_Indx", dsetname3

    integer(HID_T)                                      :: file_id, grp_id, &
                                                         & dset_id1, type_id1, space_id1, &
                                                         & dset_id2, type_id2, space_id2, &
                                                         & dset_id3, type_id3, space_id3

    integer(HSIZE_T)                                    :: dims1(rank1), dims2(rank2), dims3(rank3)
    integer(HSIZE_T)                                    :: maxdims1(rank1), maxdims2(rank2), maxdims3(rank3)

    real(dp), allocatable, dimension(:,:,:,:,:,:,:)     :: FC4_part
    integer                                             :: mu, Nbasis, Natm
    character(len=24)                                   :: mu_char
    character(len=8)                                    :: frmt

    integer                                             :: hdferr


    write(*, 155) filename
    155 FORMAT('Reading Fourth-Order Force constant file: ', A)

    call h5open_f(hdferr)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdferr) !output -> fileid


    call h5dopen_f(file_id, dsetname1, dset_id1, hdferr) !Opens an existing dataset.
                                                         !output -> dset_id1
    call h5dopen_f(file_id, dsetname2, dset_id2, hdferr) !Opens an existing dataset.
                                                         !output -> dset_id2

    call h5dget_type_f(dset_id1, type_id1, hdferr) !output -> type_id1
    call h5dget_type_f(dset_id2, type_id2, hdferr) !output -> type_id2

    call h5dget_space_f(dset_id1, space_id1, hdferr) !output -> space_id1
    call h5sget_simple_extent_dims_f(space_id1, dims1, maxdims1, hdferr) !output -> dims1
                                                                                  ! maxdims1
    call h5dget_space_f(dset_id2, space_id2, hdferr) !output -> space_id2
    call h5sget_simple_extent_dims_f(space_id2, dims2, maxdims2, hdferr) !output -> dims2
                                                                                  ! maxdims2
    allocate(atmNum4(1:dims1(1)))
    allocate(atmIndx4(1:dims2(1), 1:dims2(2), 1:dims2(3)))

    call h5dread_f(dset_id1, type_id1, atmNum4, dims1, hdferr) !Read in array 'atmNum'
    call h5dread_f(dset_id2, type_id2, atmIndx4, dims2, hdferr) !Read in array 'atmIndx'

    call h5sclose_f( space_id1, hdferr )
    call h5sclose_f( space_id2, hdferr )
    
    call h5tclose_f( type_id1, hdferr )
    call h5tclose_f( type_id2, hdferr )

    call h5dclose_f(dset_id1, hdferr)
    call h5dclose_f(dset_id2, hdferr)

    Nbasis = int( dims1(1) )
    Natm = maxval(atmNum4)
    allocate( FC_4th(3, 3, 3, 3, Natm, Natm, Natm, Nbasis) )
    FC_4th = 0.0_dp
    !..............................******************************................................!
    frmt = '(I2)'

    call h5gopen_f(file_id, fc_grp_name, grp_id, hdferr) !output -> grp_id

        mu_loop: do mu = 1, Nbasis

            write(mu_char, frmt) mu 
            dsetname3 = 'FC4_mu'//trim(adjustl(adjustr(mu_char)))

            call h5dopen_f(grp_id, dsetname3, dset_id3, hdferr) !Opens an existing dataset.
                                                                !output -> dset_id3
            call h5dget_type_f(dset_id3, type_id3, hdferr) !output -> type_id3
            call h5dget_space_f(dset_id3, space_id3, hdferr) !output -> space_id3
            call h5sget_simple_extent_dims_f(space_id3, dims3, maxdims3, hdferr) !output -> dims3

            allocate( FC4_part(1:dims3(1), 1:dims3(2), 1:dims3(3), 1:dims3(4), &
                              &1:dims3(5), 1:dims3(6), 1:dims3(7)) )

            call h5dread_f(dset_id3, type_id3, FC4_part, dims3, hdferr) 

            FC_4th(:,:,:,:,:,:,:, mu) = FC4_part(:,:,:,:,:,:,:)

            call h5sclose_f( space_id3, hdferr )
            call h5tclose_f( type_id3, hdferr )

            call h5dclose_f(dset_id3, hdferr)

            deallocate( FC4_part )

        end do mu_loop

    call h5gclose_f(grp_id, hdferr)
    !..............................******************************................................!

    call h5fclose_f(file_id, hdferr)
    call h5close_f(hdferr)

end subroutine r4thFC


subroutine w2dh5_no_grp(arr, filename)

    implicit none

    integer, parameter                      :: rank = 2 !for 2D matrix rank is always 2
    real(dp), dimension(:,:), intent(in)    :: arr
    character(len=*) , intent(in)           :: filename
    character(len=100)                      :: dsetname="freq_dat" !the default dataset name

    integer                                 :: h5err
    integer(hid_t)                          :: file_id, dspace_id, dset_id
    integer(hsize_t), dimension(1:2)        :: dims 

    dims = (/size(arr, 1), size(arr, 2)/)

    write(*, 255) filename
    255 FORMAT('Phonon dispersion is written in file: ', A)

    call h5open_f(h5err)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)  !output -> file_id

    call h5screate_simple_f(rank, dims, dspace_id, h5err) !output -> dspace_id

    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, &
                     dset_id, h5err) !output -> dset_id

    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, arr, dims, h5err)

    call h5dclose_f(dset_id, h5err)
    call h5sclose_f(dspace_id, h5err)
    call h5fclose_f(file_id, h5err)
    call h5close_f(h5err)

end subroutine w2dh5_no_grp


subroutine r_force_disp(filename, disp, force)

    implicit none

    character(len=*), intent(in)                                    :: filename
    real(dp), dimension(:,:,:,:,:,:), allocatable, intent(out)      :: disp
    real(dp), dimension(:,:,:,:,:,:), allocatable, intent(out)      :: force

    integer, parameter                                              :: rank=6
    character(len=50)                                               :: disp_dset='disp_dataset', &
                                                                       force_dset='force_dataset'

    real(dp), dimension(:,:,:,:,:,:), allocatable                   :: force_tmp, &
                                                                     & disp_tmp

    integer(HID_T)                                                  :: file_id, dset_id1, &
                                                                      & type_id1, space_id1, &
                                                                      & dset_id2, type_id2, space_id2!, hdferr
    integer(HSIZE_T)                                                :: dims1(rank), dims2(rank)
    integer(HSIZE_T)                                                :: maxdims1(rank), maxdims2(rank)
    integer                                                         :: hdferr, ts


    write(*, 125) filename
    125 FORMAT('Reading force displacemt dataset file: ', A64)

    call h5open_f(hdferr)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdferr) !output -> fileid

    call h5dopen_f(file_id, disp_dset, dset_id1, hdferr)  !Opens an existing dataset.
                                                          !output -> dset_id1
    call h5dopen_f(file_id, force_dset, dset_id2, hdferr) !Opens an existing dataset.
                                                          !output -> dset_id2
    call h5dget_type_f(dset_id1, type_id1, hdferr) !output -> type_id1
    call h5dget_type_f(dset_id2, type_id2, hdferr) !output -> type_id2

    call h5dget_space_f(dset_id1, space_id1, hdferr) !output -> space_id1
    call h5sget_simple_extent_dims_f(space_id1, dims1, maxdims1, hdferr) !output -> dims1
                                                                                  ! maxdims1
    call h5dget_space_f(dset_id2, space_id2, hdferr) !output -> space_id2
    call h5sget_simple_extent_dims_f(space_id2, dims2, maxdims2, hdferr) !output -> dims2
                                                                                  ! maxdims2
    allocate(disp_tmp(1:dims1(1), 1:dims1(2), 1:dims1(3), 1:dims1(4), 1:dims1(5), 1:dims1(6)))
    allocate(force_tmp(1:dims2(1), 1:dims2(2), 1:dims2(3), 1:dims2(4), 1:dims2(5), 1:dims2(6)))

    call h5dread_f(dset_id1, type_id1, disp_tmp, dims1, hdferr)  !Read in array 'disp_tmp'
    call h5dread_f(dset_id2, type_id2, force_tmp, dims2, hdferr) !Read in array 'force_tmp'

    call h5sclose_f( space_id1, hdferr )
    call h5sclose_f( space_id2, hdferr )

    call h5tclose_f( type_id1, hdferr )
    call h5tclose_f( type_id2, hdferr )

    call h5dclose_f(dset_id1, hdferr)
    call h5dclose_f(dset_id2, hdferr)
    call h5fclose_f(file_id, hdferr)
    call h5close_f(hdferr)

    allocate( disp(1:dims1(6), 1:dims1(1), 1:dims1(2), 1:dims1(3), 1:dims1(4), 1:dims1(5)) )
    allocate( force(1:dims2(6), 1:dims2(1), 1:dims2(2), 1:dims2(3), 1:dims2(4), 1:dims2(5)) )

    ts_loop: do ts = 1, int(dims1(6))
        disp(ts, :, :, :, :, :) = disp_tmp(:, :, :, :, :, ts)
        force(ts, :, :, :, :, :) = force_tmp(:, :, :, :, :, ts)
    end do ts_loop

    deallocate( disp_tmp, force_tmp )

end subroutine r_force_disp


subroutine wFC2_Renorm(filename, FC2_Renorm, atmNum, atm_Indx, &
                     & force_renorm, force_md)

    implicit none

    integer, parameter                                      :: rank1 = 1, rank3 = 3, &
                                                             & rank5 = 5

    character(len=*) , intent(in)                           :: filename
    real(dp), dimension(:,:,:,:,:), intent(in)              :: FC2_Renorm
    integer, dimension(:), intent(in)                       :: atmNum
    integer, dimension(:,:,:), intent(in)                   :: atm_Indx
    real(dp), dimension(:), intent(in)                      :: force_renorm, force_md

    !................................... Local variable .....................................!

    character(len=100)                                      :: dsetname1='NumAtoms', &
                                                             & dsetname2='FC2', &
                                                             & dsetname3='atm_Indx', &
                                                             & dsetname4="force_renorm", &
                                                             & dsetname5="force_md"

                                                             !& fc_grp_name='FC2nd'
    
    integer(hid_t)                                          :: file_id , &!grp_id, &
                                                             & dspace_id1, dset_id1, &
                                                             & dspace_id2, dset_id2, &
                                                             & dspace_id3, dset_id3, &
                                                             & dspace_id4, dset_id4, &
                                                             & dspace_id5, dset_id5


    integer(hsize_t), dimension(1:rank1)                    :: dims1
    integer(hsize_t), dimension(1:rank3)                    :: dims3
    integer(hsize_t), dimension(1:rank5)                    :: dims5
    integer(hsize_t), dimension(1:rank1)                    :: dims_fr
    integer(hsize_t), dimension(1:rank1)                    :: dims_fmd

    !integer                                                 :: mu, Nbasis
    !character(len=24)                                       :: mu_char
    !character(len=8)                                        :: frmt

    integer                                                 :: h5err

    !................................... Local variable .....................................!

    !Nbasis = size( atmNum )

    dims1 = size( atmNum )
    dims_fr = size(force_renorm)
    dims_fmd = size(force_md)
    dims3 = shape( atm_Indx )

    dims5 = shape( FC2_Renorm )

    write(*, 225) filename
    225 FORMAT('Writing Renormalize FC2 in file: ', A64)

    !frmt = '(I2)'

    call h5open_f(h5err)

        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)  !output -> file_id


            call h5screate_simple_f(rank1, dims1, dspace_id1, h5err) !output -> dspace_id1
            call h5screate_simple_f(rank3, dims3, dspace_id3, h5err) !output -> dspace_id3
            call h5screate_simple_f(rank1, dims_fr, dspace_id4, h5err) !output -> dspace_id4
            call h5screate_simple_f(rank1, dims_fmd, dspace_id5, h5err) !output -> dspace_id5


            call h5dcreate_f(file_id, dsetname1, H5T_NATIVE_INTEGER, dspace_id1, &
                             dset_id1, h5err) !output -> dset_id1
            call h5dcreate_f(file_id, dsetname3, H5T_NATIVE_INTEGER, dspace_id3, &
                             dset_id3, h5err) !output -> dset_id3
            call h5dcreate_f(file_id, dsetname4, H5T_NATIVE_DOUBLE, dspace_id4, &
                             dset_id4, h5err) !output -> dset_id4
            call h5dcreate_f(file_id, dsetname5, H5T_NATIVE_DOUBLE, dspace_id5, &
                             dset_id5, h5err) !output -> dset_id5


            call h5dwrite_f(dset_id1, H5T_NATIVE_INTEGER, atmNum, dims1, h5err)
            call h5dwrite_f(dset_id3, H5T_NATIVE_INTEGER, atm_Indx, dims3, h5err)
            call h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, force_renorm, dims_fr, h5err)
            call h5dwrite_f(dset_id5, H5T_NATIVE_DOUBLE, force_md, dims_fmd, h5err)

            call h5screate_simple_f(rank5, dims5, dspace_id2, h5err) !output -> dspace_id2
            call h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id2, &
                             dset_id2, h5err) !output -> dset_id2

            call h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, FC2_Renorm, dims5, h5err)

            call h5dclose_f(dset_id2, h5err)
            call h5sclose_f(dspace_id2, h5err)

            !!================================= FC2 write ===================================!
            !call h5gcreate_f(file_id, fc_grp_name, grp_id, h5err) !output -> grp_id

            !    mu_loop: do mu = 1, Nbasis

            !        write(mu_char, frmt) mu
            !        dsetname2 = 'FC2_mu'//trim(adjustl(adjustr(mu_char)))

            !        call h5screate_simple_f(rank4, dims4, dspace_id2, h5err) !output -> dspace_id2
            !        call h5dcreate_f(grp_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id2, &
            !                         dset_id2, h5err) !output -> dset_id2

            !        call h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, &
            !                      & FC3_Reconst(:,:,:,mu) , dims4, h5err)

            !        call h5dclose_f(dset_id2, h5err)
            !        call h5sclose_f(dspace_id2, h5err)

            !    end do mu_loop

            !call h5gclose_f(grp_id, h5err)
            !!================================= FC2 write ===================================!

            call h5dclose_f(dset_id5, h5err)
            call h5dclose_f(dset_id4, h5err)
            call h5dclose_f(dset_id3, h5err)
            call h5dclose_f(dset_id1, h5err)

            call h5sclose_f(dspace_id5, h5err)
            call h5sclose_f(dspace_id4, h5err)
            call h5sclose_f(dspace_id3, h5err)
            call h5sclose_f(dspace_id1, h5err)

        call h5fclose_f(file_id, h5err)

    call h5close_f(h5err)

end subroutine wFC2_Renorm

end module hdf5_wrap

