
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

    use hdf5,       only : H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER,&
                         & H5F_ACC_TRUNC_F, H5F_ACC_RDWR_F, h5open_f, h5fcreate_f, &
                         & h5screate_simple_f, h5dcreate_f,  h5dwrite_f, &
                         & h5dclose_f, h5sclose_f, h5fclose_f, h5close_f,&
                         & HID_T, HSIZE_T, H5F_ACC_RDONLY_F, h5dget_type_f, &
                         & h5dget_space_f, h5sget_simple_extent_dims_f, h5dread_f, &
                         & h5dopen_f, h5fopen_f, h5gopen_f, h5gclose_f, h5tclose_f
    use kinds,      only : dp

    implicit none
    private
    public :: r2ndFC, r3rdFC, r4thFC, w_RestartLW, w_Restart, r_Restart, r_RestartLW, &
            & w1dh5_no_grp, w2dh5_no_grp, w2dh5_no_grp_type2, w_ScatterMatEl4th, r_ScatterMatEl4th

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


    img1_chk: if ( this_image() == 1 ) then
        write(*, 155) filename
        155 FORMAT('Reading Second-Order Force constant file: ', A)
    end if img1_chk

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


subroutine r3rdFC(filename, atmNum3, atmIndx3, FC_3rd)

    implicit none

    character(len=*), intent(in)        :: filename
    integer, allocatable, intent(out)   :: atmNum3(:)
    integer, allocatable, intent(out)   :: atmIndx3(:, :, :)
    real(dp), allocatable, intent(out)  :: FC_3rd(:,:,:,:,:,:)

    integer(HID_T), parameter                           :: rank1 = 1, rank2 = 3, rank3 = 5

    character(len=50)                                   :: fc_grp_name="FC3rd", &
                                                         & dsetname1="NumAtoms", &
                                                         & dsetname2="atm_Indx", dsetname3

    integer(HID_T)                                      :: file_id, grp_id, &
                                                         & dset_id1, type_id1, space_id1, &
                                                         & dset_id2, type_id2, space_id2, &
                                                         & dset_id3, type_id3, space_id3

    integer(HSIZE_T)                                    :: dims1(rank1), dims2(rank2), dims3(rank3)
    integer(HSIZE_T)                                    :: maxdims1(rank1), maxdims2(rank2), maxdims3(rank3)

    real(dp), allocatable, dimension(:,:,:,:,:)         :: FC3_part
    integer                                             :: mu, Nbasis, Natm
    character(len=24)                                   :: mu_char
    character(len=8)                                    :: frmt

    integer                                             :: hdferr


    img1_chk: if ( this_image() == 1 ) then
        write(*, 155) filename
        155 FORMAT('Reading Third-Order Force constant file: ', A)
    end if img1_chk

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
        allocate(atmNum3(1:dims1(1)))
        allocate(atmIndx3(1:dims2(1), 1:dims2(2), 1:dims2(3)))

        call h5dread_f(dset_id1, type_id1, atmNum3, dims1, hdferr) !Read in array 'atmNum'
        call h5dread_f(dset_id2, type_id2, atmIndx3, dims2, hdferr) !Read in array 'atmIndx'

        call h5sclose_f( space_id1, hdferr )
        call h5sclose_f( space_id2, hdferr )

        call h5tclose_f( type_id1, hdferr )
        call h5tclose_f( type_id2, hdferr )

        call h5dclose_f(dset_id1, hdferr)
        call h5dclose_f(dset_id2, hdferr)

        Nbasis = int( dims1(1) )
        Natm = maxval( atmNum3 )
        allocate( FC_3rd(3, 3, 3, Natm, Natm, Nbasis) )
        FC_3rd = 0.0_dp
        !..............................******************************................................!
        frmt = '(I2)'

        call h5gopen_f(file_id, fc_grp_name, grp_id, hdferr) !output -> grp_id

            mu_loop: do mu = 1, Nbasis

                write(mu_char, frmt) mu 
                dsetname3 = 'FC3_mu'//trim(adjustl(adjustr(mu_char)))

                call h5dopen_f(grp_id, dsetname3, dset_id3, hdferr) !Opens an existing dataset.
                                                                    !output -> dset_id3
                call h5dget_type_f(dset_id3, type_id3, hdferr) !output -> type_id3
                call h5dget_space_f(dset_id3, space_id3, hdferr) !output -> space_id3
                call h5sget_simple_extent_dims_f(space_id3, dims3, maxdims3, hdferr) !output -> dims3

                allocate( FC3_part(1:dims3(1), 1:dims3(2), 1:dims3(3), 1:dims3(4), 1:dims3(5)) )

                call h5dread_f(dset_id3, type_id3, FC3_part, dims3, hdferr) 

                FC_3rd(:,:,:,:,:, mu) = FC3_part(:,:,:,:,:)

                call h5sclose_f( space_id3, hdferr )
                call h5tclose_f( type_id3, hdferr )

                call h5dclose_f(dset_id3, hdferr)

                deallocate( FC3_part )

            end do mu_loop

        call h5gclose_f(grp_id, hdferr)
        !..............................******************************................................!

    call h5fclose_f(file_id, hdferr)
    call h5close_f(hdferr)

end subroutine r3rdFC


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


    img1_chk: if ( this_image() == 1 ) then
        write(*, 155) filename
        155 FORMAT('Reading Fourth-Order Force constant file: ', A)
    end if img1_chk

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


subroutine w_RestartLW(filename, Ndof, Nq0points, i0, complete, arr)

    implicit none

    integer, parameter                                  :: rank1 = 1, rank2 = 2 

    character(len=*), intent(in)                        :: filename

    integer, intent(in)                                 :: Ndof, Nq0points
    integer, intent(in)                                 :: i0, complete

    real(dp), dimension(Ndof, Nq0points), intent(in)    :: arr

    ! ================================== Local Variables ================================== !

    character(len=100)                      :: dsetname_m = "LW_data", &
                                             & prgrs_dset = "progress"

    integer, dimension(2)                   :: progress

    integer                                 :: h5err
    integer(hid_t)                          :: file_id, dspace_id, dset_id, &
                                             & dspace_idm, dset_idm
    integer(hsize_t), dimension(1:rank1)    :: dims1 
    integer(hsize_t), dimension(1:rank2)    :: dims2 

    ! ================================== Local Variables ================================== !

    progress(:) = (/i0, complete/)

    dims1 = 2
    dims2 = (/Ndof, Nq0points/)

    call h5open_f(h5err)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)  !output -> file_id
    
        ! --------------------------------------- progress --------------------------------------- !
        call h5screate_simple_f( rank1, dims1, dspace_id, h5err ) !output -> dspace_id

        call h5dcreate_f( file_id, prgrs_dset, H5T_NATIVE_INTEGER, dspace_id, &
                        & dset_id, h5err ) !output -> dset_id
        call h5dwrite_f( dset_id, H5T_NATIVE_INTEGER, progress, dims1, h5err )

        call h5dclose_f( dset_id, h5err )
        call h5sclose_f( dspace_id, h5err )
        ! --------------------------------------- progress --------------------------------------- !

        ! -------------------------------------- Linewidth --------------------------------------- !
        call h5screate_simple_f(rank2, dims2, dspace_idm, h5err) !output -> dspace_idm

        call h5dcreate_f(file_id, dsetname_m, H5T_NATIVE_DOUBLE, dspace_idm, &
                       & dset_idm, h5err) !output -> dset_idm
        call h5dwrite_f(dset_idm, H5T_NATIVE_DOUBLE, arr, dims2, h5err)

        call h5dclose_f(dset_idm, h5err)
        call h5sclose_f(dspace_idm, h5err)
        ! -------------------------------------- Linewidth --------------------------------------- !

    call h5fclose_f(file_id, h5err)
    call h5close_f(h5err)

end subroutine w_RestartLW


subroutine w_Restart( Restrt_file, Ndof, QPerProcsMax, ii, i0, complete, Sctrq1q2 )

    implicit none
    
    integer, parameter                                                  :: rank1 = 1, rank5 = 5

    character(len=*), intent(in)                                        :: Restrt_file
    integer, intent(in)                                                 :: Ndof, QPerProcsMax
    integer, intent(in)                                                 :: ii, i0, complete
    real(dp), dimension(Ndof,Ndof,Ndof,Ndof, QPerProcsMax), intent(in)  :: Sctrq1q2

    ! ===================================== Local Variables ===================================== !

    character(len=128)                                  :: prgrs_dset='progress', &
                                                         & Sctr_dset = 'ScatterMat'

    integer, dimension(3)                               :: progress

    integer                                             :: hdferr
    integer(hid_t)                                      :: file_id, dspace_id, dset_id, &
                                                         & dspace_idm, dset_idm
    integer(hsize_t), dimension(1:rank1)                :: dims1
    integer(hsize_t), dimension(1:rank5)                :: dims5

    ! ===================================== Local Variables ===================================== !

    dims1 = 3
    dims5(:) = (/Ndof,Ndof,Ndof,Ndof, QPerProcsMax/)

    progress(:) = (/ii, i0, complete/)

    call h5open_f( hdferr )
    call h5fcreate_f( Restrt_file, H5F_ACC_TRUNC_F, file_id, hdferr )   !output -> file_id

        ! --------------------------------------- progress --------------------------------------- !
        call h5screate_simple_f( rank1, dims1, dspace_id, hdferr ) !output -> dspace_id

        call h5dcreate_f( file_id, prgrs_dset, H5T_NATIVE_INTEGER, dspace_id, &
                        & dset_id, hdferr ) !output -> dset_id

        call h5dwrite_f( dset_id, H5T_NATIVE_INTEGER, progress, dims1, hdferr )

        call h5dclose_f( dset_id, hdferr )
        call h5sclose_f( dspace_id, hdferr )
        ! --------------------------------------- progress --------------------------------------- !

        ! ------------------------------ Scattering Matrix Elements ------------------------------ !
        call h5screate_simple_f( rank5, dims5, dspace_idm, hdferr ) !output -> dspace_idm

        call h5dcreate_f( file_id, Sctr_dset, H5T_NATIVE_DOUBLE, dspace_idm, &
                        & dset_idm, hdferr ) !output -> dset_idm

        call h5dwrite_f( dset_idm, H5T_NATIVE_DOUBLE, Sctrq1q2, dims5, hdferr )

        call h5dclose_f( dset_idm, hdferr )
        call h5sclose_f( dspace_idm, hdferr )
        ! ------------------------------ Scattering Matrix Elements ------------------------------ !

    call h5fclose_f( file_id, hdferr )
    call h5close_f( hdferr )

end subroutine w_Restart

subroutine r_Restart( Restrt_file, Ndof, QPerProcsMax, ii, i0, complete, Sctrq1q2 )

    implicit none
    
    integer, parameter                                                      :: rank1 = 1, rank5 = 5

    character(len=*), intent(in)                                            :: Restrt_file
    integer, intent(in)                                                     :: Ndof, QPerProcsMax

    integer, intent(out)                                                    :: ii, i0, complete
    real(dp), dimension(Ndof,Ndof,Ndof,Ndof, QPerProcsMax), intent(inout)   :: Sctrq1q2

    ! ===================================== Local Variables ===================================== !

    character(len=128)                                  :: prgrs_dset='progress', &
                                                         & Sctr_dset = 'ScatterMat'

    integer, dimension(3)                               :: progress

    integer                                             :: hdferr
    integer(hid_t)                                      :: file_id, &
                                                         & dset_id, dset_idm, &
                                                         & type_id, type_idm, &
                                                         & space_id, space_idm
    integer(hsize_t), dimension(1:rank1)                :: dims1, maxdims1
    integer(hsize_t), dimension(1:rank5)                :: dims5, maxdims5, &
                                                         & dims5_c, maxdims5_c


    ! ===================================== Local Variables ===================================== !

    dims1 = 3
    dims5(:) = (/Ndof,Ndof,Ndof,Ndof, QPerProcsMax/)

    ! ** Debug ** !
    dims5_c = shape( Sctrq1q2 )
    maxdims5_c = dims5_c
    ! ** Debug ** !

    call h5open_f( hdferr )
    call h5fopen_f( Restrt_file, H5F_ACC_RDONLY_F, file_id, hdferr )   !output -> file_id

        ! --------------------------------------- progress --------------------------------------- !
        call h5dopen_f(file_id, prgrs_dset, dset_id, hdferr) !Opens an existing dataset.
                                                             !output -> dset_id

        call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
        call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
        call h5sget_simple_extent_dims_f(space_id, dims1, maxdims1, hdferr) !output -> dims1
                                                                                     ! maxdims1
        call h5dread_f(dset_id, type_id, progress, dims1, hdferr) !Read in array 'progress'

        call h5sclose_f( space_id, hdferr )
        call h5tclose_f( type_id, hdferr )
        call h5dclose_f(dset_id, hdferr)
        ! --------------------------------------- progress --------------------------------------- !

        ! ------------------------------ Scattering Matrix Elements ------------------------------ !
        call h5dopen_f(file_id, Sctr_dset, dset_idm, hdferr) !Opens an existing dataset.
                                                             !output -> dset_idm

        call h5dget_type_f(dset_idm, type_idm, hdferr) !output -> type_idm
        call h5dget_space_f(dset_idm, space_idm, hdferr) !output -> space_idm
        call h5sget_simple_extent_dims_f(space_idm, dims5, maxdims5, hdferr) !output -> dims5
                                                                                      ! maxdims5
        Debug: if ( any((dims5_c - dims5) /= 0) .or. any((maxdims5_c - maxdims5) /= 0) ) then

            write(*, 100)
            write(*, "( '[ (', 5I6, ') /= (', 5I6, ') ] .or. ', '[ (', 5I6, ') /= (', 5I6, ') ]' )") &
                   & dims5_c, dims5, maxdims5_c, maxdims5

            100 FORMAT( "ERROR in restart read! ")

            ERROR STOP

        end if Debug

        call h5dread_f(dset_idm, type_idm, Sctrq1q2, dims5, hdferr) !Read in array 'Sctrq1q2'

        call h5sclose_f( space_idm, hdferr )
        call h5tclose_f( type_idm, hdferr )

        call h5dclose_f(dset_idm, hdferr)
        ! ------------------------------ Scattering Matrix Elements ------------------------------ !

    call h5fclose_f( file_id, hdferr )
    call h5close_f( hdferr )

    ii = progress(1)
    i0 = progress(2)
    complete = progress(3)

end subroutine r_Restart


subroutine r_RestartLW(filename, Ndof, Nq0points, progress, arr)

    implicit none

    integer, parameter                                  :: rank1 = 1, rank2 = 2 

    character(len=*) , intent(in)                       :: filename

    integer, intent(in)                                 :: Ndof, Nq0points

    integer, dimension(2), intent(out)                  :: progress
    real(dp), dimension(Ndof, Nq0points), intent(inout) :: arr

    ! ================================== Local Variables ================================== !

    character(len=100)                      :: dsetname_m = "LW_data", &
                                             & prgrs_dset = "progress"

    integer                                 :: hdferr
    integer(hid_t)                          :: file_id, &
                                             & dset_id, dset_idm, &
                                             & type_id, type_idm, &
                                             & space_id, space_idm
    integer(hsize_t), dimension(1:rank1)    :: dims1, maxdims1
    integer(hsize_t), dimension(1:rank2)    :: dims2, maxdims2, &
                                             & dims2_c, maxdims2_c

    ! ================================== Local Variables ================================== !

    dims1 = 2
    dims2 = (/Ndof, Nq0points/)

    !** Debug **!
    dims2_c = shape( arr )
    maxdims2_c = dims2_c
    !** Debug **!

    call h5open_f(hdferr)
    call h5fopen_f( filename, H5F_ACC_RDONLY_F, file_id, hdferr )   !output -> file_id
    
        ! --------------------------------------- progress --------------------------------------- !
        call h5dopen_f(file_id, prgrs_dset, dset_id, hdferr) !Opens an existing dataset.
                                                             !output -> dset_id

        call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
        call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
        call h5sget_simple_extent_dims_f(space_id, dims1, maxdims1, hdferr) !output -> dims1
                                                                                     ! maxdims1
        call h5dread_f(dset_id, type_id, progress, dims1, hdferr) !Read in array 'progress'

        call h5sclose_f( space_id, hdferr )
        call h5tclose_f( type_id, hdferr )
        call h5dclose_f(dset_id, hdferr)
        ! --------------------------------------- progress --------------------------------------- !

        ! -------------------------------------- Linewidth --------------------------------------- !
        call h5dopen_f(file_id, dsetname_m, dset_idm, hdferr) !Opens an existing dataset.
                                                              !output -> dset_idm

        call h5dget_type_f(dset_idm, type_idm, hdferr) !output -> type_idm
        call h5dget_space_f(dset_idm, space_idm, hdferr) !output -> space_idm
        call h5sget_simple_extent_dims_f(space_idm, dims2, maxdims2, hdferr) !output -> dims2
                                                                                      ! maxdims2
        Debug: if ( any((dims2_c - dims2) /= 0) .or. any((maxdims2_c - maxdims2) /= 0) ) then

            write(*, 100)
            write(*, "( '[ (', 2I6, ') /= (', 2I6, ') ] .or. ', '[ (', 2I6, ') /= (', 2I6, ') ]' )") &
                   & dims2_c, dims2, maxdims2_c, maxdims2

            100 FORMAT( "ERROR in restart read! ")

            ERROR STOP

        end if Debug

        call h5dread_f(dset_idm, type_idm, arr, dims2, hdferr) !Read in array 'arr'

        call h5sclose_f( space_idm, hdferr )
        call h5tclose_f( type_idm, hdferr )
        call h5dclose_f(dset_idm, hdferr)
        ! -------------------------------------- Linewidth --------------------------------------- !

    call h5fclose_f(file_id, hdferr)
    call h5close_f(hdferr)

end subroutine r_RestartLW


subroutine w1dh5_no_grp(arr, filename)

    implicit none

    integer, parameter                      :: rank = 1 

    real(dp), dimension(9), intent(in)      :: arr
    character(len=*) , intent(in)           :: filename

    ! ============================== Local Variables ============================== !
    character(len=128)                      :: dsetname="ThermalConductivityTensor" 

    integer                                 :: h5err
    integer(hid_t)                          :: file_id, dspace_id, dset_id
    integer(hsize_t), dimension(1:rank)     :: dims
    ! ============================== Local Variables ============================== !

    dims = (/9/)

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

end subroutine w1dh5_no_grp


subroutine w2dh5_no_grp(arr, filename, dsetname)

    implicit none

    integer, parameter                      :: rank = 2 

    real(dp), dimension(:, :), intent(in)   :: arr
    character(len=*), intent(in)            :: filename
    character(len=*), intent(in)            :: dsetname

    ! ============================== Local Variables ============================== !
    integer                                 :: h5err
    integer(hid_t)                          :: file_id, dspace_id, dset_id
    integer(hsize_t), dimension(1:rank)     :: dims
    ! ============================== Local Variables ============================== !

    dims = (/size(arr, 1), size(arr, 2)/)

    call h5open_f(h5err)
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, h5err)  !output -> file_id

    call h5screate_simple_f(rank, dims, dspace_id, h5err) !output -> dspace_id

    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, &
                     dset_id, h5err) !output -> dset_id

    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, arr, dims, h5err)

    call h5dclose_f(dset_id, h5err)
    call h5sclose_f(dspace_id, h5err)
    call h5fclose_f(file_id, h5err)
    call h5close_f(h5err)

end subroutine w2dh5_no_grp

subroutine w2dh5_no_grp_type2( Ndof, Nq0, Omega, Inarr, filename, dsetname )

    implicit none

    integer, parameter                          :: rank = 2 

    integer, intent(in)                         :: Ndof, Nq0
    real(dp), dimension(Ndof, Nq0), intent(in)  :: Omega
    real(dp), dimension(Ndof, Nq0), intent(in)  :: Inarr
    character(len=*), intent(in)                :: filename
    character(len=*), intent(in)                :: dsetname

    ! ============================== Local Variables ============================== !
    real(dp), dimension(:, :), allocatable  :: arr

    integer                                 :: h5err
    integer(hid_t)                          :: file_id, dspace_id, dset_id
    integer(hsize_t), dimension(1:rank)     :: dims
    ! ============================== Local Variables ============================== !

    allocate( arr(2*Ndof, Nq0) )
    arr(1:Ndof, :) = Omega !/ (2.0_dp * PI)
    arr((Ndof+1):, :) = Inarr

    dims = (/size(arr, 1), size(arr, 2)/)

    call h5open_f(h5err)
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, h5err)  !output -> file_id

    call h5screate_simple_f(rank, dims, dspace_id, h5err) !output -> dspace_id

    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, &
                     dset_id, h5err) !output -> dset_id

    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, arr, dims, h5err)

    call h5dclose_f(dset_id, h5err)
    call h5sclose_f(dspace_id, h5err)
    call h5fclose_f(file_id, h5err)
    call h5close_f(h5err)

    deallocate( arr )

end subroutine w2dh5_no_grp_type2


subroutine w_ScatterMatEl4th( SctrFile4th, Ndof, W_s0s1s2s3 )

    implicit none
    
    integer, parameter                                                  :: rank4 = 4

    character(len=*), intent(in)                                        :: SctrFile4th
    integer, intent(in)                                                 :: Ndof
    real(dp), dimension(Ndof,Ndof,Ndof,Ndof), intent(in)                :: W_s0s1s2s3

    ! ===================================== Local Variables ===================================== !

    character(len=128)                                  :: Sctr_dset = 'ScatterMatEl'

    integer                                             :: hdferr
    integer(hid_t)                                      :: file_id, dspace_idm, dset_idm
    integer(hsize_t), dimension(1:rank4)                :: dims4

    ! ===================================== Local Variables ===================================== !

    dims4(:) = (/Ndof, Ndof, Ndof, Ndof/)

    call h5open_f( hdferr )
    call h5fcreate_f( SctrFile4th, H5F_ACC_TRUNC_F, file_id, hdferr ) !output -> file_id

        ! ------------------------------ Scattering Matrix Elements ------------------------------ !
        call h5screate_simple_f( rank4, dims4, dspace_idm, hdferr ) !output -> dspace_idm

        call h5dcreate_f( file_id, Sctr_dset, H5T_NATIVE_DOUBLE, dspace_idm, &
                        & dset_idm, hdferr ) !output -> dset_idm

        call h5dwrite_f( dset_idm, H5T_NATIVE_DOUBLE, W_s0s1s2s3, dims4, hdferr )

        call h5dclose_f( dset_idm, hdferr )
        call h5sclose_f( dspace_idm, hdferr )
        ! ------------------------------ Scattering Matrix Elements ------------------------------ !

    call h5fclose_f( file_id, hdferr )
    call h5close_f( hdferr )

end subroutine w_ScatterMatEl4th


subroutine r_ScatterMatEl4th( SctrFile4th, Ndof, W_s0s1s2s3 )

    implicit none
    
    integer, parameter                                                      :: rank4 = 4

    character(len=*), intent(in)                                            :: SctrFile4th
    integer, intent(in)                                                     :: Ndof

    real(dp), dimension(Ndof,Ndof,Ndof,Ndof), intent(out)                   :: W_s0s1s2s3

    ! ===================================== Local Variables ===================================== !

    character(len=128)                                  :: Sctr_dset = 'ScatterMatEl'

    integer                                             :: hdferr
    integer(hid_t)                                      :: file_id, dset_idm, &
                                                         & type_idm, space_idm
    integer(hsize_t), dimension(1:rank4)                :: dims4, maxdims4, &
                                                         & dims4_c, maxdims4_c

    ! ===================================== Local Variables ===================================== !

    dims4(:) = (/Ndof,Ndof,Ndof,Ndof/)

    ! ** Debug ** !
    dims4_c = shape( W_s0s1s2s3 )
    maxdims4_c = dims4_c
    ! ** Debug ** !

    call h5open_f( hdferr )
    call h5fopen_f( SctrFile4th, H5F_ACC_RDONLY_F, file_id, hdferr )   !output -> file_id

        ! ------------------------------ Scattering Matrix Elements ------------------------------ !
        call h5dopen_f(file_id, Sctr_dset, dset_idm, hdferr) !Opens an existing dataset.
                                                             !output -> dset_idm

        call h5dget_type_f(dset_idm, type_idm, hdferr) !output -> type_idm
        call h5dget_space_f(dset_idm, space_idm, hdferr) !output -> space_idm
        call h5sget_simple_extent_dims_f(space_idm, dims4, maxdims4, hdferr) !output -> dims4
                                                                                      ! maxdims4
        Debug: if ( any((dims4_c - dims4) /= 0) .or. any((maxdims4_c - maxdims4) /= 0) ) then

            write(*, 100)
            write(*, "( '[ (', 5I6, ') /= (', 5I6, ') ] .or. ', '[ (', 5I6, ') /= (', 5I6, ') ]' )") &
                   & dims4_c, dims4, maxdims4_c, maxdims4

            100 FORMAT( "ERROR in r_ScatterMatEl4th! ")

            ERROR STOP

        end if Debug

        call h5dread_f(dset_idm, type_idm, W_s0s1s2s3, dims4, hdferr) !Read in array 'W_s0s1s2s3'

        call h5sclose_f( space_idm, hdferr )
        call h5tclose_f( type_idm, hdferr )

        call h5dclose_f( dset_idm, hdferr )
        ! ------------------------------ Scattering Matrix Elements ------------------------------ !

    call h5fclose_f( file_id, hdferr )
    call h5close_f( hdferr )

end subroutine r_ScatterMatEl4th

end module hdf5_wrap

