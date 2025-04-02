
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

module FCType

    use kinds,          only : dp
    use hdf5,           only : H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER,&
                             & H5F_ACC_TRUNC_F, h5open_f, h5fcreate_f, &
                             & h5screate_simple_f, h5dcreate_f,  h5dwrite_f, &
                             & h5dclose_f, h5sclose_f, h5fclose_f, h5close_f, h5tclose_f, &
                             & HID_T, HSIZE_T, H5F_ACC_RDONLY_F, h5dget_type_f, &
                             & h5dget_space_f, h5sget_simple_extent_dims_f, h5dread_f, &
                             & h5dopen_f, h5fopen_f, h5gopen_f, h5gclose_f

    use FC4_data_mod,   only : FC4_dat

    implicit none
    private

    type, public    :: FChunk7d

        real(dp), dimension(:,:,:,:,:,:,:), allocatable                 :: FArr7d

    end type FChunk7d

    public          :: set_FContg


contains
    
    subroutine set_FContg(FC4, filename, Farr)
    
        implicit none
    
        integer, parameter                                              :: LAST_DIM = 87

        type(FC4_dat), intent(in)                                       :: FC4
        character(len=*), intent(in)                                    :: filename
        type(FChunk7d), dimension(:,:), allocatable, intent(out)        :: Farr
    
        ! ============================= Local Variables ============================= !

        integer(HID_T), parameter                           :: rank7=7

        character(len=64)                                   :: grp_name1 = 'Final', &
                                                             & grp_name3 = 'FC_fin_list', &
                                                             & grp_name_var, dset_name_var

        integer(HID_T)                                      :: file_id, grp_id1, &
                                                             & dset_id, type_id, space_id, &
                                                             & grp_id_in, grp_id_inin

        integer(HSIZE_T)                                    :: dims7(rank7), maxdims7(rank7)
    
        integer                                             :: Nbasis, Natm_mu
        integer                                             :: mu, alpha, chunk

        real(dp), dimension(:,:,:,:,:,:,:), allocatable     :: F_part7d, Arr7tmp
        real(dp), dimension(:,:,:,:,:,:,:,:), allocatable   :: Fmupart8d

        character(len=24)                                   :: mu_char, chunk_char
        character(len=8)                                    :: frmt

        integer                                             :: hdferr
    
        ! ============================= Local Variables ============================= !

        frmt = '(I6)'
    
        Nbasis = size( FC4%atmNum )

        allocate( Farr(3, Nbasis) )

        !.......................................... Farr ............................................!

        call h5open_f(hdferr)
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdferr) !output -> file_id

            call h5gopen_f(file_id, grp_name1, grp_id1, hdferr) !Open and existing group, output -> grp_id1

                call h5gopen_f(grp_id1, grp_name3, grp_id_in, hdferr) !Open and existing group, output -> grp_id_in

                    mu_loop2: do mu = 1, Nbasis

                        Natm_mu = FC4%atmNum(mu)

                        allocate( Fmupart8d(LAST_DIM, 3, 3, 3, 3, Natm_mu, Natm_mu, Natm_mu) )

                        write(mu_char, frmt) mu
                        grp_name_var = 'FC4_mu'//trim(adjustl(adjustr(mu_char)))

                        call h5gopen_f(grp_id_in, grp_name_var, grp_id_inin, hdferr) !Open and existing group, output -> grp_id_inin

                            chunk_loop: do chunk = 1, LAST_DIM

                                write(chunk_char, frmt) chunk
                                dset_name_var = 'FC4_chunk'//trim(adjustl(adjustr(chunk_char)))

                                call h5dopen_f(grp_id_inin, dset_name_var, dset_id, hdferr) !Opens an existing dataset.
                                                                                            !output -> dset_id
                                call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
                                call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
                                call h5sget_simple_extent_dims_f(space_id, dims7, maxdims7, hdferr) !output -> dims1
                                                                                                     ! maxdims1
                                allocate( F_part7d(1:dims7(1), 1:dims7(2), 1:dims7(3), 1:dims7(4), &
                                                 & 1:dims7(5), 1:dims7(6), 1:dims7(7)) )

                                call h5dread_f(dset_id, type_id, F_part7d, dims7, hdferr) !Read in array 'F_part'
                                Fmupart8d(chunk,:,:,:,:,:,:,:) = F_part7d(:,:,:,:,:,:,:)

                                call h5sclose_f( space_id, hdferr )
                                call h5tclose_f( type_id, hdferr )

                                call h5dclose_f(dset_id, hdferr)
                                deallocate(F_part7d )

                            end do chunk_loop

                        call h5gclose_f(grp_id_inin, hdferr)

                        allocate( Arr7tmp(LAST_DIM, 3, 3, 3, Natm_mu, Natm_mu, Natm_mu) )

                        alpha_loop: do alpha = 1, 3

                            allocate( Farr(alpha, mu)%FArr7d(LAST_DIM, 3, 3, 3, Natm_mu, Natm_mu, Natm_mu) )

                            Arr7tmp = Fmupart8d(:,:,:,:, alpha,:,:,:)

                            Farr(alpha, mu)%FArr7d(:,:,:,:, :,:,:) = Arr7tmp

                        end do alpha_loop

                        deallocate( Arr7tmp )
                        deallocate( Fmupart8d )

                    end do mu_loop2
                call h5gclose_f(grp_id_in, hdferr)

            call h5gclose_f(grp_id1, hdferr)
        call h5fclose_f(file_id, hdferr)
        call h5close_f(hdferr)
        !.......................................... Farr ............................................!

    end subroutine set_FContg

end module FCType 

