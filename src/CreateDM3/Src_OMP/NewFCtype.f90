
module FCType

    use kinds,          only : dp
    use hdf5,           only : H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER,&
                             & H5F_ACC_TRUNC_F, h5open_f, h5fcreate_f, &
                             & h5screate_simple_f, h5dcreate_f,  h5dwrite_f, &
                             & h5dclose_f, h5sclose_f, h5fclose_f, h5close_f, h5tclose_f, &
                             & HID_T, HSIZE_T, H5F_ACC_RDONLY_F, h5dget_type_f, &
                             & h5dget_space_f, h5sget_simple_extent_dims_f, h5dread_f, &
                             & h5dopen_f, h5fopen_f, h5gopen_f, h5gclose_f

    use FC3_data_mod,   only : FC3_dat

    implicit none
    private

    public          :: set_FContg


contains
    
    subroutine set_FContg(FC3, filename, Farr)
    
        implicit none
    
        integer, parameter                                  :: LAST_DIM = 32

        type(FC3_dat), intent(in)                                       :: FC3
        character(len=*), intent(in)                                    :: filename
        real(dp), dimension(:,:,:,:,:,:,:), allocatable, intent(out)    :: Farr
    
        ! ============================= Local Variables ============================= !

        integer(HID_T), parameter                           :: rank6=6

        character(len=64)                                   :: grp_name1 = 'Final', &
                                                             & grp_name3 = 'FC_fin_list', &
                                                             & dset_name_var

        integer(HID_T)                                      :: file_id, grp_id1, &
                                                             & dset_id, type_id, space_id, &
                                                             & grp_id_in

        integer(HSIZE_T)                                    :: dims6(rank6), maxdims6(rank6)
    
        integer                                             :: Nbasis, NatmMax, Natm_mu
        integer                                             :: mu, alpha

        real(dp), dimension(:,:,:,:,:), allocatable         :: F_alpha_part
        real(dp), dimension(:,:,:,:,:,:), allocatable       :: F_part

        character(len=24)                                   :: mu_char
        character(len=8)                                    :: frmt

        integer                                             :: hdferr
    
        ! ============================= Local Variables ============================= !

        frmt = '(I6)'
    
        Nbasis = size( FC3%atmNum )
        NatmMax = maxval( FC3%atmNum )

        allocate( Farr(LAST_DIM,3,3,NatmMax,NatmMax,3,Nbasis) )
        Farr = 0.0_dp

        !.......................................... Farr ............................................!

        call h5open_f(hdferr)
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdferr) !output -> file_id

            call h5gopen_f(file_id, grp_name1, grp_id1, hdferr) !Open and existing group, output -> grp_id1

                call h5gopen_f(grp_id1, grp_name3, grp_id_in, hdferr) !Open and existing group, output -> grp_id_in

                    mu_loop2: do mu = 1, Nbasis

                        Natm_mu = FC3%atmNum(mu)

                        write(mu_char, frmt) mu
                        dset_name_var = 'FC3_mu'//trim(adjustl(adjustr(mu_char)))

                        call h5dopen_f(grp_id_in, dset_name_var, dset_id, hdferr) !Opens an existing dataset.
                                                                                  !output -> dset_id
                        call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
                        call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
                        call h5sget_simple_extent_dims_f(space_id, dims6, maxdims6, hdferr) !output -> dims6
                                                                                            ! maxdims6
                        allocate( F_part(1:dims6(1), 1:dims6(2), 1:dims6(3), 1:dims6(4), &
                                       & 1:dims6(5), 1:dims6(6)) )

                        call h5dread_f(dset_id, type_id, F_part, dims6, hdferr) !Read in array 'F_part'

                        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
                        allocate( F_alpha_part(LAST_DIM, 3, 3, Natm_mu, Natm_mu) )
                        alpha_loop: do alpha = 1, 3

                            F_alpha_part = F_part(:,:,:, alpha, :,:)

                            Farr(:,:,:,1:Natm_mu,1:Natm_mu, alpha, mu) = F_alpha_part

                        end do alpha_loop
                        deallocate( F_alpha_part )
                        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

                        call h5sclose_f( space_id, hdferr )
                        call h5tclose_f( type_id, hdferr )

                        call h5dclose_f(dset_id, hdferr)
                        deallocate(F_part )

                    end do mu_loop2

                call h5gclose_f(grp_id_in, hdferr)

            call h5gclose_f(grp_id1, hdferr)
        call h5fclose_f(file_id, hdferr)
        call h5close_f(hdferr)
        !.......................................... Farr ............................................!

    end subroutine set_FContg

end module FCType 

