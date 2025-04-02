
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


module FC3_data_mod

    use kinds,          only : dp
    use hdf5,           only : hid_t, hsize_t, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER,&
                             & H5F_ACC_TRUNC_F, h5open_f, h5fcreate_f, &
                             & h5screate_simple_f, h5dcreate_f,  h5dwrite_f, &
                             & h5dclose_f, h5sclose_f, h5fclose_f, h5close_f, h5tclose_f, &
                             & HID_T, HSIZE_T, H5F_ACC_RDONLY_F, h5dget_type_f, &
                             & h5dget_space_f, h5sget_simple_extent_dims_f, h5dread_f, &
                             & h5dopen_f, h5fopen_f, h5gopen_f, h5gclose_f
    implicit none
    private

    ! ============================== Sparse ============================== !
    type, public    :: RedInFCsp

        integer, dimension(:), allocatable                          :: Indx
        real(dp), dimension(:), allocatable                         :: Coeff

    end type RedInFCsp
    ! ============================== Sparse ============================== !

    type, public    :: FC3_dat

        integer                                                     :: ind_fc
        integer, dimension(:), allocatable                          :: atmNum
        integer, dimension(:,:,:), allocatable                      :: cb_Indx
        integer, dimension(:,:,:,:,:,:), allocatable                :: PosIndx

        integer, dimension(:, :), allocatable                       :: PosIndx_flat_old
        integer, dimension(:,:,:,:,:,:), allocatable                :: PosIndx_old

        !real(dp), dimension(:, :), allocatable                      :: RedInFC
        ! ============================== Sparse ============================== !
        integer                                                     :: dep_fc
        integer, dimension(:), allocatable                          :: DepIndxMap
        type(RedInFCsp), dimension(:), allocatable                  :: SprMat
        ! ============================== Sparse ============================== !

        real(dp), dimension(:,:,:,:,:,:,:), allocatable             :: F

    contains
        procedure, public, pass         :: set_FC3dat

    end type FC3_dat

contains

subroutine set_FC3dat(this, filename)

    implicit none

    class(FC3_dat)                                                          :: this
    character(len=*), intent(in)                                            :: filename

    !..................................... Local variable .....................................!
    integer, parameter                                      :: LAST_DIM = 32
    integer(HID_T), parameter                               :: rank1=1, rank3=3, rank2=2, &
                                                             & rank5=5, rank6=6

    character(len=64)                                       :: grp_name1 = 'Final', &
                                                             & grp_name2 = 'PosIndx_fin_list', &
                                                             & grp_name3 = 'FC_fin_list', &
                                                             & grp_spMat = 'DepFC_sp', &                ! Sparse
                                                             & grp_Indx = 'IndxSet', &                  ! Sparse
                                                             & grp_Coeff = 'CoefSet', &                 ! Sparse
                                                             & grp_name4 = 'Old', &
                                                             & grp_name5 = 'PosIndx_init_list'

    character(len=64)                                       :: dset_name1 = 'NumFC', &
                                                             & dset_name2 = 'NumAtoms', &
                                                             & dset_name3 = 'atm_Indx', &
                                                             & dset_depfc = 'dep_fc', &                 ! Sparse
                                                             & dset_DepMap = 'DepIndxMap', &            ! Sparse
                                                             & dset_name5 = 'PosIndx_flat_init', &
                                                             & dset_name_C, dset_name_I, &              ! Sparse
                                                             & dset_name_var

    integer(HID_T)                                          :: file_id, grp_id1, &
                                                             & dset_id, type_id, space_id, &
                                                             & grp_id_in, &
                                                             & grp_id_sp, grp_id_Indx, grp_id_Coeff, &  ! Sparse
                                                             & dset_id_I, dset_id_C, &                  ! Sparse
                                                             & type_id_I, type_id_C, &                  ! Sparse
                                                             & space_id_I, space_id_C                   ! Sparse

    integer(HSIZE_T)                                        :: dims1(rank1), dims3(rank3), &
                                                             & dims2(rank2), dims5(rank5), &
                                                             & dims6(rank6), &
                                                             & dimsC(rank1), dimsI(rank1)               ! Sparse
                                                             
    integer(HSIZE_T)                                        :: maxdims1(rank1), maxdims3(rank3), &
                                                             & maxdims2(rank2), maxdims5(rank5), &
                                                             & maxdims6(rank6), &
                                                             & maxdimsC(rank1), maxdimsI(rank1)         ! Sparse

    integer, dimension(:,:,:,:,:), allocatable              :: PosIndx_part
    real(dp), dimension(:,:,:,:,:,:), allocatable           :: F_part

    integer                                                 :: Nbasis, Natm, mu, &
                                                             & Ndep, ii                                 ! Sparse

    character(len=24)                                       :: mu_char, &
                                                             & Indx_char                                ! Sparse

    character(len=8)                                        :: frmt

    integer                                                 :: hdferr
    !..................................... Local variable .....................................!

    img1_chk: if ( this_image() == 1 ) then
        write(*, 132) filename
        132 FORMAT('Reading Third order FC data from file: ', A)
    end if img1_chk

    call h5open_f(hdferr)

    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdferr) !output -> file_id
        call h5gopen_f(file_id, grp_name1, grp_id1, hdferr) !Open and existing group, output -> grp_id1
            
            !....................................... ind_fc ..........................................!
            call h5dopen_f(grp_id1, dset_name1, dset_id, hdferr) !Opens an existing dataset.
                                                                 !output -> dset_id
            call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
            call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
            dims1 = 1

            call h5dread_f(dset_id, type_id, this%ind_fc, dims1, hdferr) !Read in array 'read_data'

            call h5sclose_f( space_id, hdferr )
            call h5tclose_f( type_id, hdferr )

            call h5dclose_f(dset_id, hdferr)
            !....................................... ind_fc ..........................................!

            !....................................... atmNum ..........................................!
            call h5dopen_f(grp_id1, dset_name2, dset_id, hdferr) !Opens an existing dataset.
                                                                 !output -> dset_id
            call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
            call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
            call h5sget_simple_extent_dims_f(space_id, dims1, maxdims1, hdferr) !output -> dims1
                                                                                ! maxdims1
            allocate( this%atmNum(1:dims1(1)) )
            call h5dread_f(dset_id, type_id, this%atmNum, dims1, hdferr) !Read in array 'atmNum'

            call h5sclose_f( space_id, hdferr )
            call h5tclose_f( type_id, hdferr )

            call h5dclose_f(dset_id, hdferr)
            !....................................... atmNum ..........................................!

            !...................................... cb_Indx ..........................................!
            call h5dopen_f(grp_id1, dset_name3, dset_id, hdferr) !Opens an existing dataset.
                                                                 !output -> dset_id
            call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
            call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
            call h5sget_simple_extent_dims_f(space_id, dims3, maxdims3, hdferr) !output -> dims3
                                                                                ! maxdims3
            allocate( this%cb_Indx(1:dims3(1), 1:dims3(2), 1:dims3(3)) )
            call h5dread_f(dset_id, type_id, this%cb_Indx, dims3, hdferr) !Read in array 'cb_indx'

            call h5sclose_f( space_id, hdferr )
            call h5tclose_f( type_id, hdferr )

            call h5dclose_f(dset_id, hdferr)
            !...................................... cb_Indx ..........................................!

            !...................................... RedInFC ..........................................!
            !*! call h5dopen_f(grp_id1, dset_name4, dset_id, hdferr) !Opens an existing dataset.
            !*!                                                      !output -> dset_id
            !*! call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
            !*! call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
            !*! call h5sget_simple_extent_dims_f(space_id, dims2, maxdims2, hdferr) !output -> dims1
            !*!                                                                      ! maxdims1
            !*! allocate( this%RedInFC(1:dims2(1), 1:dims2(2)) )
            !*! call h5dread_f(dset_id, type_id, this%RedInFC, dims2, hdferr) !Read in array 'RedInFC'
            !*! call h5sclose_f( space_id, hdferr )
            !*! call h5tclose_f( type_id, hdferr )
            !*! call h5dclose_f(dset_id, hdferr)
            !...................................... RedInFC ..........................................!

            !................................... RedInFC_sparse ......................................!

            call h5gopen_f(grp_id1, grp_spMat, grp_id_sp, hdferr) !Open and existing group, output -> grp_id_sp
            
                !!                          ********** dep_fc **********                            !!
                call h5dopen_f(grp_id_sp, dset_depfc, dset_id, hdferr) !Opens an existing dataset.
                                                                       !output -> dset_id
                call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
                call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
                dims1 = 1

                call h5dread_f(dset_id, type_id, Ndep, dims1, hdferr) !Read in array 'Ndep'
                this%dep_fc = Ndep

                call h5sclose_f( space_id, hdferr )
                call h5tclose_f( type_id, hdferr )

                call h5dclose_f(dset_id, hdferr)
                !!                          ********** dep_fc **********                            !!

                !!                        ********** DepIndxMap **********                          !!
                call h5dopen_f(grp_id_sp, dset_DepMap, dset_id, hdferr) !Opens an existing dataset.
                                                                        !output -> dset_id
                call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
                call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
                call h5sget_simple_extent_dims_f(space_id, dims1, maxdims1, hdferr) !output -> dims1
                                                                                    ! maxdims1
                allocate( this%DepIndxMap(1:dims1(1)) )
                call h5dread_f(dset_id, type_id, this%DepIndxMap, dims1, hdferr) !Read in array 'DepIndxMap'

                call h5sclose_f( space_id, hdferr )
                call h5tclose_f( type_id, hdferr )

                call h5dclose_f(dset_id, hdferr)
                !!                        ********** DepIndxMap **********                          !!

                !!                    ********** CoefSet & IndxSet **********                       !!

                allocate( this%SprMat(Ndep) )
                frmt = '(I6)'

                call h5gopen_f(grp_id_sp, grp_Indx, grp_id_Indx, hdferr) !Open and existing group, output -> grp_id_Indx
                call h5gopen_f(grp_id_sp, grp_Coeff, grp_id_Coeff, hdferr) !Open and existing group, output -> grp_id_Coeff

                    depFCloop: do ii = 1, Ndep

                        write(Indx_char, frmt) ii
                        dset_name_C = 'Coeff'//trim(adjustl(adjustr(Indx_char)))
                        dset_name_I = 'Indx'//trim(adjustl(adjustr(Indx_char)))

                        call h5dopen_f(grp_id_Indx, dset_name_I, dset_id_I, hdferr) !Opens an existing dataset.
                                                                                    !output -> dset_id_I
                        call h5dopen_f(grp_id_Coeff, dset_name_C, dset_id_C, hdferr) !Opens an existing dataset.
                                                                                     !output -> dset_id_C
                        call h5dget_type_f(dset_id_I, type_id_I, hdferr) !output -> type_id_I
                        call h5dget_space_f(dset_id_I, space_id_I, hdferr) !output -> space_id_I

                        call h5dget_type_f(dset_id_C, type_id_C, hdferr) !output -> type_id_C
                        call h5dget_space_f(dset_id_C, space_id_C, hdferr) !output -> space_id_C

                        call h5sget_simple_extent_dims_f(space_id_I, dimsI, maxdimsI, hdferr) !output -> dimsI
                                                                                                       ! maxdimsI
                        call h5sget_simple_extent_dims_f(space_id_C, dimsC, maxdimsC, hdferr) !output -> dimsC
                                                                                                       ! maxdimsC
                        allocate( this%SprMat(ii)%Indx(1:dimsI(1)) )
                        call h5dread_f(dset_id_I, type_id_I, this%SprMat(ii)%Indx, dimsI, hdferr) !Read in array 'Indx'

                        allocate( this%SprMat(ii)%Coeff(1:dimsC(1)) )
                        call h5dread_f(dset_id_C, type_id_C, this%SprMat(ii)%Coeff, dimsC, hdferr) !Read in array 'Coeff'

                        call h5sclose_f( space_id_I, hdferr )
                        call h5tclose_f( type_id_I, hdferr )

                        call h5sclose_f( space_id_C, hdferr )
                        call h5tclose_f( type_id_C, hdferr )

                        call h5dclose_f(dset_id_I, hdferr)
                        call h5dclose_f(dset_id_C, hdferr)

                    end do depFCloop

                call h5gclose_f(grp_id_Indx, hdferr)
                call h5gclose_f(grp_id_Coeff, hdferr)

                !!                    ********** CoefSet & IndxSet **********                       !!

            call h5gclose_f(grp_id_sp, hdferr)

            !................................... RedInFC_sparse ......................................!

            Nbasis = size( this%atmNum )
            Natm = maxval( this%atmNum )
            frmt = '(I3)'

            !...................................... PosIndx ..........................................!
            allocate( this%PosIndx(3, 3, 3, Natm, Natm, Nbasis) )
            this%PosIndx = 0

            call h5gopen_f(grp_id1, grp_name2, grp_id_in, hdferr) !Open and existing group, output -> grp_id_in
                                                                     
                mu_loop1: do mu = 1, Nbasis

                    write(mu_char, frmt) mu
                    dset_name_var = 'PosIndx_fin_mu'//trim(adjustl(adjustr(mu_char)))

                    call h5dopen_f(grp_id_in, dset_name_var, dset_id, hdferr) !Opens an existing dataset.
                                                                              !output -> dset_id
                    call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
                    call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
                    call h5sget_simple_extent_dims_f(space_id, dims5, maxdims5, hdferr) !output -> dims1
                                                                                        ! maxdims1
                    allocate( PosIndx_part(1:dims5(1), 1:dims5(2), 1:dims5(3), 1:dims5(4), 1:dims5(5)) )

                    call h5dread_f(dset_id, type_id, PosIndx_part, dims5, hdferr) !Read in array 'PosIndx_part'
                    this%PosIndx(:,:,:,:,:,mu) = PosIndx_part(:,:,:,:,:)

                    call h5sclose_f( space_id, hdferr )
                    call h5tclose_f( type_id, hdferr )

                    call h5dclose_f(dset_id, hdferr)
                    deallocate( PosIndx_part )

                end do mu_loop1
            call h5gclose_f(grp_id_in, hdferr)
            !...................................... PosIndx ..........................................!

            !.......................................... F ............................................!
            allocate( this%F(LAST_DIM, 3, 3, 3, Natm, Natm, Nbasis) )
            this%F = 0.0_dp

            call h5gopen_f(grp_id1, grp_name3, grp_id_in, hdferr) !Open and existing group, output -> grp_id_in

                mu_loop2: do mu = 1, Nbasis

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
                    this%F(:,:,:,:,:,:,mu) = F_part(:,:,:,:,:,:)

                    call h5sclose_f( space_id, hdferr )
                    call h5tclose_f( type_id, hdferr )

                    call h5dclose_f(dset_id, hdferr)
                    deallocate(F_part )

                end do mu_loop2

            call h5gclose_f(grp_id_in, hdferr)
            !.......................................... F ............................................!

        call h5gclose_f(grp_id1, hdferr)

        call h5gopen_f(file_id, grp_name4, grp_id1, hdferr) !Open and existing group, output -> grp_id
            
            !................................. PosIndx_flat_old ......................................!
            call h5dopen_f(grp_id1, dset_name5, dset_id, hdferr) !Opens an existing dataset.
                                                                 !output -> dset_id
            call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
            call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
            call h5sget_simple_extent_dims_f(space_id, dims2, maxdims2, hdferr) !output -> dims1
                                                                                 ! maxdims1
            allocate( this%PosIndx_flat_old(1:dims2(1), 1:dims2(2)) )
            call h5dread_f(dset_id, type_id, this%PosIndx_flat_old, dims2, hdferr) !Read in array 'PosIndx_flat_old'

            call h5sclose_f( space_id, hdferr )
            call h5tclose_f( type_id, hdferr )

            call h5dclose_f(dset_id, hdferr)
            !................................. PosIndx_flat_old ......................................!

            !................................... PosIndx_old .........................................!
            allocate( this%PosIndx_old(3, 3, 3, Natm, Natm, Nbasis) )
            this%PosIndx_old = 0

            call h5gopen_f(grp_id1, grp_name5, grp_id_in, hdferr) !Open and existing group, output -> grp_id_in
                                                                     
                mu_loop3: do mu = 1, Nbasis

                    write(mu_char, frmt) mu
                    dset_name_var = 'PosIndx_init_mu'//trim(adjustl(adjustr(mu_char)))

                    call h5dopen_f(grp_id_in, dset_name_var, dset_id, hdferr) !Opens an existing dataset.
                                                                              !output -> dset_id
                    call h5dget_type_f(dset_id, type_id, hdferr) !output -> type_id
                    call h5dget_space_f(dset_id, space_id, hdferr) !output -> space_id
                    call h5sget_simple_extent_dims_f(space_id, dims5, maxdims5, hdferr) !output -> dims5
                                                                                        ! maxdims5
                    allocate( PosIndx_part(1:dims5(1), 1:dims5(2), 1:dims5(3), 1:dims5(4), 1:dims5(5)) )

                    call h5dread_f(dset_id, type_id, PosIndx_part, dims5, hdferr) !Read in array 'PosIndx_part'
                    this%PosIndx_old(:,:,:,:,:,mu) = PosIndx_part(:,:,:,:,:)

                    call h5sclose_f( space_id, hdferr )
                    call h5tclose_f( type_id, hdferr )

                    call h5dclose_f(dset_id, hdferr)
                    deallocate( PosIndx_part )

                end do mu_loop3

            call h5gclose_f(grp_id_in, hdferr)
            !................................... PosIndx_old .........................................!

        call h5gclose_f(grp_id1, hdferr)

    call h5fclose_f(file_id, hdferr)
        

    call h5close_f(hdferr)

end subroutine set_FC3dat

end module FC3_data_mod

