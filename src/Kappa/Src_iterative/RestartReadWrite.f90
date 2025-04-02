
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

subroutine WriteRestartFile( this, Nbasis, NumQ1, progress, RestartDir )

    implicit none

    integer, parameter                                          :: rank1 = 1, rank3 = 3

    class(AllQ), intent(in)                                     :: this

    integer, intent(in)                                         :: Nbasis, NumQ1
    integer, dimension(2), intent(in)                           :: progress
    character(len=*), intent(in)                                :: RestartDir

    ! ================================= Local Variables ================================= !

    character(len=512)                                  :: filename
    character(len=12)                                   :: frmt, img_num_char
    character(len=12)                                   :: qchar, q1char

    character(len=128)                                  :: prgrs_dset='progress', &
                                                         & GrpNameQ_dyn, DsetNameQ1_dyn, &
                                                         & DsetNameQ1m_dyn

    real(dp), dimension(3*Nbasis, 3*Nbasis, 3*Nbasis)   :: W_ss1s2_p, W_ss1s2_m ! Automatic Array

    integer                                             :: qi, q1i, q1end
    integer                                             :: img_num
    integer                                             :: hdferr

    integer(hid_t)                                      :: file_id, dspace_id, dset_id, &
                                                         & grp_id, dspace_idm, dset_idm
    integer(hsize_t), dimension(1:rank1)                :: dims1
    integer(hsize_t), dimension(1:rank3)                :: dims3

    ! ================================= Local Variables ================================= !

    frmt = '(I6)'

    img_num = this_image()
    write(img_num_char, frmt) img_num
    filename = trim(adjustl(adjustr(RestartDir)))//'/Img'//trim(adjustl(adjustr(img_num_char)))//'_Restrt.h5'

    write(*, 132) img_num, filename
    132 FORMAT('Writing restart file for image = ', I4, ', in file: ', A64, ' ...')

    dims3 = (/3*Nbasis, 3*Nbasis, 3*Nbasis/)

    call h5open_f( hdferr )
    call h5fcreate_f( filename, H5F_ACC_TRUNC_F, file_id, hdferr )   !output -> file_id

        ! --------------------------------------- progress --------------------------------------- !
        dims1 = 2
        call h5screate_simple_f( rank1, dims1, dspace_id, hdferr ) !output -> dspace_id

        call h5dcreate_f( file_id, prgrs_dset, H5T_NATIVE_INTEGER, dspace_id, &
                        & dset_id, hdferr ) !output -> dset_id

        call h5dwrite_f( dset_id, H5T_NATIVE_INTEGER, progress, dims1, hdferr )

        call h5dclose_f( dset_id, hdferr )
        call h5sclose_f( dspace_id, hdferr )
        ! --------------------------------------- progress --------------------------------------- !

        q_loop: do qi = 1, progress(1)

            if ( qi == progress(1) ) then
                q1end = progress(2)

            else
                q1end = NumQ1

            end if

            write(qchar, frmt) qi
            GrpNameQ_dyn = 'ScatterProbQ_'//trim(adjustl(adjustr(qchar)))

            call h5gcreate_f( file_id, GrpNameQ_dyn, grp_id, hdferr ) !output -> grp_id

                q1_loop: do q1i = 1, q1end

                    W_ss1s2_p = this%WPlus(qi)%q1q2(q1i)%ss1s2
                    W_ss1s2_m = this%WMinus(qi)%q1q2(q1i)%ss1s2

                    write(q1char, frmt) q1i
                    DsetNameQ1_dyn = 'ScatterProbQ1p_'//trim(adjustl(adjustr(q1char)))
                    DsetNameQ1m_dyn = 'ScatterProbQ1m_'//trim(adjustl(adjustr(q1char)))

                    call h5screate_simple_f( rank3, dims3, dspace_id, hdferr ) !output -> dspace_id
                    call h5screate_simple_f( rank3, dims3, dspace_idm, hdferr ) !output -> dspace_idm

                    call h5dcreate_f( grp_id, DsetNameQ1_dyn, H5T_NATIVE_DOUBLE, dspace_id, &
                                    & dset_id, hdferr ) !output -> dset_id
                    call h5dcreate_f( grp_id, DsetNameQ1m_dyn, H5T_NATIVE_DOUBLE, dspace_idm, &
                                    & dset_idm, hdferr ) !output -> dset_idm

                    call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, W_ss1s2_p, dims3, hdferr )
                    call h5dwrite_f( dset_idm, H5T_NATIVE_DOUBLE, W_ss1s2_m, dims3, hdferr )

                    call h5dclose_f( dset_id, hdferr )
                    call h5sclose_f( dspace_id, hdferr )
                    call h5dclose_f( dset_idm, hdferr )
                    call h5sclose_f( dspace_idm, hdferr )

                end do q1_loop

            call h5gclose_f( grp_id, hdferr )

        end do q_loop

    call h5fclose_f( file_id, hdferr )
    call h5close_f( hdferr )

    write(*, 142) img_num, filename
    142 FORMAT('Done writing restart file for image = ', I4, ', in file: ', A64)

end subroutine WriteRestartFile


subroutine ReadRestartFile( this, Nbasis, NumQ1, RestartDir, progress )

    implicit none

    integer, parameter                                          :: rank1 = 1, rank3 = 3

    class(AllQ)                                                 :: this

    integer, intent(in)                                         :: Nbasis, NumQ1
    character(len=*), intent(in)                                :: RestartDir
    integer, dimension(2), intent(out)                          :: progress

    ! ================================= Local Variables ================================= !

    character(len=512)                                  :: filename
    character(len=12)                                   :: frmt, img_num_char
    character(len=12)                                   :: qchar, q1char

    character(len=128)                                  :: prgrs_dset='progress', &
                                                         & GrpNameQ_dyn, DsetNameQ1_dyn, &
                                                         & DsetNameQ1m_dyn

    real(dp), dimension(3*Nbasis, 3*Nbasis, 3*Nbasis)   :: W_ss1s2_p, W_ss1s2_m ! Automatic Array

    integer                                             :: qi, q1i, q1end
    integer                                             :: img_num
    integer                                             :: hdferr

    integer(hid_t)                                      :: file_id, dset_id, &
                                                         & grp_id, dset_idm, &
                                                         & type_id, type_idm, &
                                                         & space_id, space_idm
    integer(hsize_t), dimension(1:rank1)                :: dims1, maxdims1
    integer(hsize_t), dimension(1:rank3)                :: dims3, maxdims3, dims3m, maxdims3m

    ! ================================= Local Variables ================================= !

    frmt = '(I6)'

    img_num = this_image()
    write(img_num_char, frmt) img_num
    filename = trim(adjustl(adjustr(RestartDir)))//'/Img'//trim(adjustl(adjustr(img_num_char)))//'_Restrt.h5'

    write(*, 132) img_num, filename
    132 FORMAT('Reading restart file for image = ', I4, ', from file: ', A64, ' ...')

    dims3 = (/3*Nbasis, 3*Nbasis, 3*Nbasis/)

    call h5open_f( hdferr )
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

        q_loop: do qi = 1, progress(1)

            if ( qi == progress(1) ) then
                q1end = progress(2)

            else
                q1end = NumQ1

            end if

            allocate( this%WPlus(qi)%q1q2( NumQ1 ) )
            allocate( this%WMinus(qi)%q1q2( NumQ1 ) )

            write(qchar, frmt) qi
            GrpNameQ_dyn = 'ScatterProbQ_'//trim(adjustl(adjustr(qchar)))

            call h5gopen_f( file_id, GrpNameQ_dyn, grp_id, hdferr ) !output -> grp_id

                q1_loop: do q1i = 1, q1end

                    write(q1char, frmt) q1i
                    DsetNameQ1_dyn = 'ScatterProbQ1p_'//trim(adjustl(adjustr(q1char)))
                    DsetNameQ1m_dyn = 'ScatterProbQ1m_'//trim(adjustl(adjustr(q1char)))

                    call h5dopen_f(grp_id, DsetNameQ1_dyn, dset_id, hdferr) !Opens an existing dataset.
                                                                            !output -> dset_id
                    call h5dopen_f(grp_id, DsetNameQ1m_dyn, dset_idm, hdferr) !Opens an existing dataset.
                                                                              !output -> dset_idm

                    call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
                    call h5dget_type_f(dset_idm, type_idm, hdferr) !output -> type_idm

                    call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
                    call h5dget_space_f(dset_idm, space_idm, hdferr) !output -> space_idm

                    call h5sget_simple_extent_dims_f(space_id, dims3, maxdims3, hdferr) !output -> dims3
                                                                                                 ! maxdims3
                    call h5sget_simple_extent_dims_f(space_idm, dims3m, maxdims3m, hdferr) !output -> dims3m
                                                                                                    ! maxdims3m

                    call h5dread_f(dset_id, type_id, W_ss1s2_p, dims3, hdferr) !Read in array 'W_ss1s2_p'
                    call h5dread_f(dset_idm, type_idm, W_ss1s2_m, dims3m, hdferr) !Read in array 'W_ss1s2_m'

                    allocate( this%WPlus(qi)%q1q2(q1i)%ss1s2(1:dims3(1), 1:dims3(2), 1:dims3(3)) )
                    this%WPlus(qi)%q1q2(q1i)%ss1s2 = W_ss1s2_p

                    allocate( this%WMinus(qi)%q1q2(q1i)%ss1s2(1:dims3m(1), 1:dims3m(2), 1:dims3m(3)) )
                    this%WMinus(qi)%q1q2(q1i)%ss1s2 = W_ss1s2_m

                    call h5sclose_f( space_id, hdferr )
                    call h5sclose_f( space_idm, hdferr )

                    call h5tclose_f( type_id, hdferr )
                    call h5tclose_f( type_idm, hdferr )

                    call h5dclose_f(dset_id, hdferr)
                    call h5dclose_f(dset_idm, hdferr)

                end do q1_loop

            call h5gclose_f( grp_id, hdferr )

        end do q_loop

    call h5fclose_f( file_id, hdferr )
    call h5close_f( hdferr )

    write(*, 142) img_num, filename
    142 FORMAT('Done reading restart file for image = ', I4, ', from file: ', A64)

end subroutine ReadRestartFile

