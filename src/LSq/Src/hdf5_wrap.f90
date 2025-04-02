
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

    use hdf5,               only : H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER,&
                                 & H5F_ACC_TRUNC_F, H5F_ACC_RDWR_F, h5open_f, h5fcreate_f, &
                                 & h5screate_simple_f, h5dcreate_f,  h5dwrite_f, h5gget_info_by_name_f, &
                                 & h5dclose_f, h5sclose_f, h5fclose_f, h5close_f, h5gunlink_f, &
                                 & HID_T, HSIZE_T, H5F_ACC_RDONLY_F, h5dget_type_f, h5tclose_f, &
                                 & h5dget_space_f, h5sget_simple_extent_dims_f, h5dread_f, &
                                 & h5dopen_f, h5fopen_f, h5gcreate_f, h5gopen_f, h5gclose_f
    use kinds,              only : dp

    implicit none
    private

    public  :: ReadDispMat, ReadForcedd, ReadForceRenorm, WriteIndFC, ReadFlags

contains

    subroutine ReadDispMat(filename, Mat, force, MatShape)

        implicit none

        integer, parameter                                      :: rank1=1, rank2=2

        character(len=*), intent(in)                            :: filename
        real(dp), dimension(:,:), allocatable, intent(out)      :: Mat
        real(dp), dimension(:), allocatable, intent(out)        :: force
        integer, dimension(2), intent(out)                      :: MatShape

        ! ===================================== Local Variables ===================================== !

        real(dp), dimension(:,:), allocatable                   :: MatTmp

        character(len=128)                                      :: grp_name='DispForceDset', &
                                                                 & dsetname1='disp_mat', &
                                                                 & dsetname2='force_row'

        integer(hid_t)                                          :: file_id, grp_id, &
                                                                 & dset_id1, dset_id2, &
                                                                 & type_id1, type_id2, &
                                                                 & space_id1, space_id2

        integer(hsize_t), dimension(1:rank1)                    :: dims1, maxdims1
        integer(hsize_t), dimension(1:rank2)                    :: dims2, maxdims2

        integer                                                 :: h5err

        ! ===================================== Local Variables ===================================== !

        write(*, 100) filename
        100 FORMAT("Reading file: ", A64)

        call h5open_f(h5err)
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, h5err) !output -> file_id
        call h5gopen_f(file_id, grp_name, grp_id, h5err) !output -> grp_id
            
            call h5dopen_f(grp_id, dsetname1, dset_id1, h5err) !Opens an existing dataset.
                                                               !output -> dset_id1
            call h5dopen_f(grp_id, dsetname2, dset_id2, h5err) !Opens an existing dataset.
                                                               !output -> dset_id2
            
                call h5dget_type_f(dset_id1, type_id1, h5err) !output -> type_id1
                call h5dget_type_f(dset_id2, type_id2, h5err) !output -> type_id2
            
                call h5dget_space_f(dset_id1, space_id1, h5err) !output -> space_id1
                call h5dget_space_f(dset_id2, space_id2, h5err) !output -> space_id2
            
                call h5sget_simple_extent_dims_f(space_id1, dims2, maxdims2, h5err) !output -> dims2
                                                                                    ! maxdims2
                call h5sget_simple_extent_dims_f(space_id2, dims1, maxdims1, h5err) !output -> dims1
                                                                                    ! maxdims1
            
                allocate( MatTmp(1:dims2(1), 1:dims2(2)) )
                call h5dread_f(dset_id1, type_id1, MatTmp, dims2, h5err) !Read in array 'MatTmp'
            
                allocate( force(1:dims1(1)) )
                call h5dread_f(dset_id2, type_id2, force, dims1, h5err) !Read in array 'force'

                call h5sclose_f( space_id1, h5err )
                call h5sclose_f( space_id2, h5err )

                call h5tclose_f( type_id1, h5err )
                call h5tclose_f( type_id2, h5err )
            
            call h5dclose_f(dset_id2, h5err)
            call h5dclose_f(dset_id1, h5err)

        call h5gclose_f(grp_id, h5err)
        call h5fclose_f(file_id, h5err)
        call h5close_f(h5err)

        allocate( Mat( size(MatTmp, 2), size(MatTmp, 1) ) )
        Mat = transpose( MatTmp )
        deallocate( MatTmp )

        MatShape = shape( Mat )

    end subroutine ReadDispMat


    subroutine ReadForcedd(filename, Forcedd, Forcemd)

        implicit none

        integer, parameter                                      :: rank1=1

        character(len=*), intent(in)                            :: filename
        real(dp), dimension(:), allocatable, intent(out)        :: Forcedd, Forcemd

        ! ===================================== Local Variables ===================================== !

        character(len=128)                                      :: dsetname1='force_dd', &
                                                                 & dsetname2='force_md'

        integer(hid_t)                                          :: file_id, &
                                                                 & dset_id1, dset_id2, &
                                                                 & type_id1, type_id2, &
                                                                 & space_id1, space_id2

        integer(hsize_t), dimension(1:rank1)                    :: dims1, maxdims1
        integer(hsize_t), dimension(1:rank1)                    :: dims2, maxdims2

        integer                                                 :: h5err

        ! ===================================== Local Variables ===================================== !

        write(*, 100) filename
        100 FORMAT("Reading file: ", A64)

        call h5open_f(h5err)
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, h5err) !output -> file_id
            
            call h5dopen_f(file_id, dsetname1, dset_id1, h5err) !Opens an existing dataset.
                                                                !output -> dset_id1
            call h5dopen_f(file_id, dsetname2, dset_id2, h5err) !Opens an existing dataset.
                                                                !output -> dset_id2
                call h5dget_type_f(dset_id1, type_id1, h5err) !output -> type_id1
                call h5dget_type_f(dset_id2, type_id2, h5err) !output -> type_id2
            
                call h5dget_space_f(dset_id1, space_id1, h5err) !output -> space_id1
                call h5dget_space_f(dset_id2, space_id2, h5err) !output -> space_id2
            
                call h5sget_simple_extent_dims_f(space_id1, dims1, maxdims1, h5err) !output -> dims1
                                                                                    ! maxdims1
                call h5sget_simple_extent_dims_f(space_id2, dims2, maxdims2, h5err) !output -> dims2
                                                                                    ! maxdims2
            
                allocate( Forcedd(1:dims1(1)) )
                call h5dread_f(dset_id1, type_id1, Forcedd, dims2, h5err) !Read in array 'Forcedd'
            
                allocate( Forcemd(1:dims2(1)) )
                call h5dread_f(dset_id2, type_id2, Forcemd, dims1, h5err) !Read in array 'Forcemd'

                call h5sclose_f( space_id1, h5err )
                call h5sclose_f( space_id2, h5err )

                call h5tclose_f( type_id1, h5err )
                call h5tclose_f( type_id2, h5err )
            
            call h5dclose_f(dset_id2, h5err)
            call h5dclose_f(dset_id1, h5err)

        call h5fclose_f(file_id, h5err)
        call h5close_f(h5err)

    end subroutine ReadForcedd


    subroutine WriteIndFC(filename, grp_name, ind_fc) 

        implicit none
        
        integer, parameter                          :: rank1=1

        character(len=*), intent(in)                :: filename, grp_name
        real(dp), dimension(:), intent(in)          :: ind_fc

        ! ===================================== Local Variables ===================================== !

        character(len=28)                           :: dsetname = "indFC"

        integer(hid_t)                              :: file_id, grp_id, &
                                                     & dspace_id, dset_id
        
        integer(hsize_t), dimension(1:rank1)        :: dims1

        !-! integer                                     :: storage_type, nlinks, max_corder
        integer                                     :: h5err

        ! ===================================== Local Variables ===================================== !

        write(*, *)
        write(*, 120) filename

        dims1 = size(ind_fc)

        call h5open_f(h5err)
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, h5err) !output -> fileid

        call h5gcreate_f(file_id, grp_name, grp_id, h5err) !output -> grp_id

        Exists: if ( h5err /= 0 ) then

            write(*, *)
            write(*, 100)
            write(*, 150) filename
            write(*, 160) h5err
            write(*, 200)
            write(*, 100)

            call h5gunlink_f(file_id, grp_name, h5err)
            call h5gcreate_f(file_id, grp_name, grp_id, h5err) !output -> grp_id
           
        end if Exists

            call h5screate_simple_f(rank1, dims1, dspace_id, h5err) !output -> dspace_id
            
            call h5dcreate_f(grp_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, &
                           & dset_id, h5err) !output -> dset_id
            
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ind_fc, dims1, h5err)
            
            call h5dclose_f(dset_id, h5err)
            call h5sclose_f(dspace_id, h5err)

        call h5gclose_f(grp_id, h5err)
        call h5fclose_f(file_id, h5err)
        call h5close_f(h5err)

        100 FORMAT("============----------------------------------------------------------============")
        120 FORMAT("Writing independent FC in file: ", A64)
        150 FORMAT("Unable to create independent-FC group in file ", A64)
        160 FORMAT("Probably the group already exists. Error code: ", I3)
        200 FORMAT("Unlinking the existing group and overwriting (IGNORE THE ABOVE MESSAGES FROM HDF5)")


    end subroutine WriteIndFC


    subroutine ReadForceRenorm(filename, ForceRenorm, Forcemd)

        implicit none

        integer, parameter                                      :: rank1=1

        character(len=*), intent(in)                            :: filename
        real(dp), dimension(:), allocatable, intent(out)        :: ForceRenorm, Forcemd

        ! ===================================== Local Variables ===================================== !

        character(len=128)                                      :: dsetname1='force_renorm', &
                                                                 & dsetname2='force_md'

        integer(hid_t)                                          :: file_id, &
                                                                 & dset_id1, dset_id2, &
                                                                 & type_id1, type_id2, &
                                                                 & space_id1, space_id2

        integer(hsize_t), dimension(1:rank1)                    :: dims1, maxdims1
        integer(hsize_t), dimension(1:rank1)                    :: dims2, maxdims2

        integer                                                 :: h5err

        ! ===================================== Local Variables ===================================== !

        write(*, 100) filename
        100 FORMAT("Reading force from renormalized 2nd-IFCs. Filename: ", A64)

        call h5open_f(h5err)
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, h5err) !output -> file_id
            
            call h5dopen_f(file_id, dsetname1, dset_id1, h5err) !Opens an existing dataset.
                                                                !output -> dset_id1
            call h5dopen_f(file_id, dsetname2, dset_id2, h5err) !Opens an existing dataset.
                                                                !output -> dset_id2
                call h5dget_type_f(dset_id1, type_id1, h5err) !output -> type_id1
                call h5dget_type_f(dset_id2, type_id2, h5err) !output -> type_id2
            
                call h5dget_space_f(dset_id1, space_id1, h5err) !output -> space_id1
                call h5dget_space_f(dset_id2, space_id2, h5err) !output -> space_id2
            
                call h5sget_simple_extent_dims_f(space_id1, dims1, maxdims1, h5err) !output -> dims1
                                                                                    ! maxdims1
                call h5sget_simple_extent_dims_f(space_id2, dims2, maxdims2, h5err) !output -> dims2
                                                                                    ! maxdims2
            
                allocate( ForceRenorm(1:dims1(1)) )
                call h5dread_f(dset_id1, type_id1, ForceRenorm, dims2, h5err) !Read in array 'ForceRenorm'
            
                allocate( Forcemd(1:dims2(1)) )
                call h5dread_f(dset_id2, type_id2, Forcemd, dims1, h5err) !Read in array 'Forcemd'

                call h5sclose_f( space_id1, h5err )
                call h5sclose_f( space_id2, h5err )

                call h5tclose_f( type_id1, h5err )
                call h5tclose_f( type_id2, h5err )
            
            call h5dclose_f(dset_id2, h5err)
            call h5dclose_f(dset_id1, h5err)

        call h5fclose_f(file_id, h5err)
        call h5close_f(h5err)

    end subroutine ReadForceRenorm


    subroutine ReadFlags( filename, flags )

        implicit none

        integer, parameter                                      :: rank1=1

        character(len=*), intent(in)                            :: filename
        integer, dimension(:), allocatable, intent(out)         :: flags

        ! ===================================== Local Variables ===================================== !

        character(len=28)                                       :: dsetname1='flags_row'

        integer(hid_t)                                          :: file_id, dset_id1, &
                                                                 & type_id1, space_id1

        integer(hsize_t), dimension(1:rank1)                    :: dims1, maxdims1

        integer                                                 :: h5err

        ! ===================================== Local Variables ===================================== !

        write(*, 100) filename
        100 FORMAT("Reading file: ", A64)

        call h5open_f(h5err)
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, h5err) !output -> file_id
            
            call h5dopen_f(file_id, dsetname1, dset_id1, h5err) !Opens an existing dataset.
                                                                !output -> dset_id1
                call h5dget_type_f(dset_id1, type_id1, h5err) !output -> type_id1
                call h5dget_space_f(dset_id1, space_id1, h5err) !output -> space_id1
            
                call h5sget_simple_extent_dims_f(space_id1, dims1, maxdims1, h5err) !output -> dims1
                                                                                    ! maxdims1
                allocate( flags(1:dims1(1)) )
                call h5dread_f(dset_id1, type_id1, flags, dims1, h5err) !Read in array 'flags'

                call h5sclose_f( space_id1, h5err )
                call h5tclose_f( type_id1, h5err )
            
            call h5dclose_f(dset_id1, h5err)

        call h5fclose_f(file_id, h5err)
        call h5close_f(h5err)

    end subroutine ReadFlags

end module hdf5_wrap

