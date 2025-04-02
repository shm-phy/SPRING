
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
    public :: r2ndFC, wCellBasisDisp

contains

    subroutine r2ndFC(filename, atmNum, atmIndx, FC2nd)
    
        implicit none

        integer, parameter                  :: rank1 = 1, rank2 = 3, rank3 = 5
    
        character(len=*), intent(in)        :: filename
        integer, allocatable, intent(out)   :: atmNum(:)
        integer, allocatable, intent(out)   :: atmIndx(:, :, :)
        real(dp), allocatable, intent(out)  :: FC2nd(:, :, :, :, :)
    
        ! ====================================== Local Variables ====================================== !

        character(len=50)                   :: dsetname1="NumAtoms", dsetname2="atm_Indx", dsetname3="FC2"
    
        integer(HID_T)                      :: file_id, dset_id1, type_id1, space_id1, &
                                                        dset_id2, type_id2, space_id2, &
                                                        dset_id3, type_id3, space_id3
        integer(HSIZE_T)                    :: dims1(rank1), dims2(rank2), dims3(rank3)
        integer(HSIZE_T)                    :: maxdims1(rank1), maxdims2(rank2), maxdims3(rank3)
        integer                             :: hdferr

        ! ====================================== Local Variables ====================================== !
    
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
    
    
    subroutine wCellBasisDisp(cell_basis_record, InitPos, dU, filename)
    
        implicit none
    
        integer, parameter                              :: rank2=2, rank6=6

        integer, dimension(:,:), intent(in)             :: cell_basis_record
        real(dp), dimension(:,:), intent(in)            :: InitPos
        real(dp), dimension(:,:,:,:,:,:), intent(in)    :: dU
        character(len=*), intent(in)                    :: filename
    
        ! ====================================== Local Variables ====================================== !
    
        integer, dimension(:,:), allocatable            :: cb_tmp
    
        character(len=128)                              :: disp_dset='disp_dataset', &
                                                         & cb_dset='cell_basis', cart_dset='cart_coord'
    
        integer(hid_t)                                  :: file_id, &
                                                        & dspace_id1, dset_id1, &
                                                        & dspace_id2, dset_id2, &
                                                        & dspace_id3, dset_id3
        
        integer(hsize_t), dimension(1:rank2)            :: dims2
        integer(hsize_t), dimension(1:rank2)            :: dims2_2nd
        integer(hsize_t), dimension(1:rank6)            :: dims6

        integer                                         :: h5err
    
        ! ====================================== Local Variables ====================================== !
    
        dims2 = shape( cell_basis_record )
        dims2_2nd = shape( InitPos )
        dims6 = shape( dU )
    
        allocate( cb_tmp(dims2(1), dims2(2)) )
        cb_tmp = cell_basis_record
        cb_tmp(4, :) = cb_tmp(4, :) - 1
    
        write(*, 225) filename
        225 FORMAT('Writing displacement data in file: ', A64)
    
        call h5open_f(h5err)
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)  !output -> file_id
    
            call h5screate_simple_f(rank2, dims2, dspace_id1, h5err) !output -> dspace_id1
            call h5screate_simple_f(rank2, dims2_2nd, dspace_id3, h5err) !output -> dspace_id3
            call h5screate_simple_f(rank6, dims6, dspace_id2, h5err) !output -> dspace_id2
    
            call h5dcreate_f(file_id, cb_dset, H5T_NATIVE_INTEGER, dspace_id1, &
                           & dset_id1, h5err) !output -> dset_id1
            call h5dcreate_f(file_id, cart_dset, H5T_NATIVE_DOUBLE, dspace_id3, &
                           & dset_id3, h5err) !output -> dset_id3
            call h5dcreate_f(file_id, disp_dset, H5T_NATIVE_DOUBLE, dspace_id2, &
                           & dset_id2, h5err) !output -> dset_id2
    
            call h5dwrite_f(dset_id1, H5T_NATIVE_INTEGER, cb_tmp, dims2, h5err)
            call h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, InitPos, dims2_2nd, h5err)
            call h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, dU, dims6, h5err)
    
            call h5dclose_f(dset_id2, h5err)
            call h5dclose_f(dset_id3, h5err)
            call h5dclose_f(dset_id1, h5err)
    
            call h5sclose_f(dspace_id2, h5err)
            call h5sclose_f(dspace_id3, h5err)
            call h5sclose_f(dspace_id1, h5err)
    
        call h5fclose_f(file_id, h5err)
        call h5close_f(h5err)
    
        deallocate( cb_tmp )
    
    end subroutine wCellBasisDisp

end module hdf5_wrap

