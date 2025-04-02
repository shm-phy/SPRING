
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
                         & h5dopen_f, h5fopen_f
    use kinds,      only : dp
    implicit none
    private
    public :: wFC2_dd, r_force_disp

contains

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


    subroutine wFC2_dd(FC2_dd, Rcb, force_dd, force_md, filename)
    
        implicit none
    
        integer, parameter                                      :: rank1 = 4, rank2 = 3, &
                                                                 & rank3 = 1
    
        real(dp), contiguous, dimension(:,:,:,:), intent(in)    :: FC2_dd
        integer, contiguous, dimension(:,:,:), intent(in)       :: RCb
        real(dp), contiguous, dimension(:), intent(in)          :: force_dd, force_md
        character(len=*) , intent(in)                           :: filename
    
        character(len=100)                                      :: dsetname1="FC2_dd", & 
                                                                 & dsetname2="N2_indx", & 
                                                                 & dsetname3="force_dd", &
                                                                 & dsetname4="force_md"
    
        integer                                                 :: h5err
    
        integer(hid_t)                                          :: file_id, &
                                                                 & dspace_id1, dset_id1, &
                                                                 & dspace_id2, dset_id2, &
                                                                 & dspace_id3, dset_id3, &
                                                                 & dspace_id4, dset_id4
    
        integer(hsize_t), dimension(1:rank1)                        :: dims1
        integer(hsize_t), dimension(1:rank2)                        :: dims2
        integer(hsize_t), dimension(1:rank3)                        :: dims3, dims4
    
        dims1 = shape(FC2_dd)
        dims2 = shape(Rcb)
        dims3 = size(force_dd)
        dims4 = size(force_md)
    
        write(*, 225) filename
        225 FORMAT('Writing d-d FC and longrange force in file: ', A64)
    
        call h5open_f(h5err)
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)  !output -> file_id
    
        call h5screate_simple_f(rank1, dims1, dspace_id1, h5err) !output -> dspace_id1
        call h5screate_simple_f(rank2, dims2, dspace_id2, h5err) !output -> dspace_id2
        call h5screate_simple_f(rank3, dims3, dspace_id3, h5err) !output -> dspace_id3
        call h5screate_simple_f(rank3, dims4, dspace_id4, h5err) !output -> dspace_id4
    
        call h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id1, &
                         dset_id1, h5err) !output -> dset_id1
        call h5dcreate_f(file_id, dsetname2, H5T_NATIVE_INTEGER, dspace_id2, &
                         dset_id2, h5err) !output -> dset_id2
        call h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE, dspace_id3, &
                         dset_id3, h5err) !output -> dset_id3
        call h5dcreate_f(file_id, dsetname4, H5T_NATIVE_DOUBLE, dspace_id4, &
                         dset_id4, h5err) !output -> dset_id4
    
    
        call h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, FC2_dd, dims1, h5err)
        call h5dwrite_f(dset_id2, H5T_NATIVE_INTEGER, RCb, dims2, h5err)
        call h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, force_dd, dims3, h5err)
        call h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, force_md, dims4, h5err)
    
        call h5dclose_f(dset_id4, h5err)
        call h5dclose_f(dset_id3, h5err)
        call h5dclose_f(dset_id2, h5err)
        call h5dclose_f(dset_id1, h5err)

        call h5sclose_f(dspace_id4, h5err)
        call h5sclose_f(dspace_id3, h5err)
        call h5sclose_f(dspace_id2, h5err)
        call h5sclose_f(dspace_id1, h5err)

        call h5fclose_f(file_id, h5err)
        call h5close_f(h5err)
    
    end subroutine wFC2_dd

end module hdf5_wrap

