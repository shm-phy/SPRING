
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
    public :: r2ndFC, w2dh5_no_grp

contains

    subroutine r2ndFC(filename, atmNum, atmIndx, FC2)
    
        implicit none
    
        character(len=*), intent(in)        :: filename
        integer, allocatable, intent(out)   :: atmNum(:)
        integer, allocatable, intent(out)   :: atmIndx(:, :, :)
        real(dp), allocatable, intent(out)  :: FC2(:, :, :, :, :)
    
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
        allocate(FC2(1:dims3(1),1:dims3(2), 1:dims3(3), 1:dims3(4), 1:dims3(5)))
    
        call h5dread_f(dset_id1, type_id1, atmNum, dims1, hdferr) !Read in array 'atmNum'
        call h5dread_f(dset_id2, type_id2, atmIndx, dims2, hdferr) !Read in array 'atmIndx'
        call h5dread_f(dset_id3, type_id3, FC2, dims3, hdferr) !Read in array 'FC2'

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

end module hdf5_wrap

