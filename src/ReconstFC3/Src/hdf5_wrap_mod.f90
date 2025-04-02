
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
                         & H5F_ACC_TRUNC_F, H5F_ACC_RDWR_F, h5open_f, h5fcreate_f, &
                         & h5screate_simple_f, h5dcreate_f,  h5dwrite_f, &
                         & h5dclose_f, h5sclose_f, h5fclose_f, h5close_f, h5tclose_f, &
                         & HID_T, HSIZE_T, H5F_ACC_RDONLY_F, h5dget_type_f, &
                         & h5dget_space_f, h5sget_simple_extent_dims_f, h5dread_f, &
                         & h5dopen_f, h5fopen_f, h5gcreate_f, h5gopen_f, h5gclose_f

    use kinds,      only : dp
    use unit_cell,  only : cell

    implicit none
    private
    public :: r_force_disp, wFD_part, wFC3_Reconst, ShengBTE_FC3rd

contains

subroutine r_force_disp(filename, disp, force)

    implicit none

    character(len=*), intent(in)                                    :: filename
    real(dp), dimension(:,:,:,:,:,:), allocatable, intent(out)      :: disp
    real(dp), dimension(:,:,:,:,:,:), allocatable, intent(out)      :: force

    integer, parameter          :: rank=6
    character(len=50)           :: disp_dset='disp_dataset', &
                                   force_dset='force_dataset'
    integer(HID_T)              :: file_id, dset_id1, type_id1, space_id1, &
                                   dset_id2, type_id2, space_id2!, hdferr
    integer(HSIZE_T)            :: dims1(rank), dims2(rank)
    integer(HSIZE_T)            :: maxdims1(rank), maxdims2(rank)
    integer                     :: hdferr


    img1_chk: if ( this_image() == 1 ) then
        write(*, 125) filename
        125 FORMAT('Reading force displacemt dataset file: ', A64)
    end if img1_chk

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
    allocate(disp(1:dims1(1), 1:dims1(2), 1:dims1(3), 1:dims1(4), 1:dims1(5), 1:dims1(6)))
    allocate(force(1:dims2(1), 1:dims2(2), 1:dims2(3), 1:dims2(4), 1:dims2(5), 1:dims2(6)))

    call h5dread_f(dset_id1, type_id1, disp, dims1, hdferr)  !Read in array 'disp'
    call h5dread_f(dset_id2, type_id2, force, dims2, hdferr) !Read in array 'force'

    call h5sclose_f( space_id1, hdferr )
    call h5sclose_f( space_id2, hdferr )
    
    call h5tclose_f( type_id1, hdferr )
    call h5tclose_f( type_id2, hdferr )

    call h5dclose_f(dset_id1, hdferr)
    call h5dclose_f(dset_id2, hdferr)
    call h5fclose_f(file_id, hdferr)
    call h5close_f(hdferr)

end subroutine r_force_disp


subroutine wFD_part(Mat, force_mua, filename, mu, alpha)

    implicit none

    integer, parameter                                      :: rank1 = 2, rank2 = 1

    real(dp), contiguous, dimension(:,:), intent(in)        :: Mat
    real(dp), contiguous, dimension(:), intent(in)          :: force_mua
    character(len=*) , intent(in)                           :: filename
    integer, intent(in)                                     :: mu, alpha

    !................................... Local variable .....................................!
    real(dp), dimension(:, :), allocatable                  :: Mat_trans

    character(len=100)                                      :: grp_name='DispForceDset_parts', &
                                                             & dsetname1, dsetname2

    integer(hid_t)                                          :: file_id, grp_id, &
                                                             & dspace_id1, dset_id1, &
                                                             & dspace_id2, dset_id2

    integer(hsize_t), dimension(1:rank1)                    :: dims1
    integer(hsize_t), dimension(1:rank2)                    :: dims2

    character(len=8)                                        :: frmt
    character(len=24)                                       :: mu_char, alpha_char
    integer                                                 :: h5err, mua
    !................................... Local variable .....................................!

    frmt = '(I3)'
    write(mu_char, frmt) mu
    write(alpha_char, frmt) alpha

    dsetname1 = 'disp_mat_mu'//trim(adjustl(adjustr(mu_char)))//'_alpha'//trim(adjustl(adjustr(alpha_char)))
    dsetname2 = 'force_row_mu'//trim(adjustl(adjustr(mu_char)))//'_alpha'//trim(adjustl(adjustr(alpha_char)))

    allocate( Mat_trans(size(Mat, 2), size(Mat, 1)) )
    Mat_trans = transpose( Mat )

    dims1 = shape(Mat_trans)
    dims2 = size(force_mua)

    img1_chk: if ( this_image() == 1 ) then
        write(*, 225) filename
        225 FORMAT('Writing FD data in file: ', A64)
    end if img1_chk

    mua = 3*(mu-1) + (alpha-1)  + 1

    call h5open_f(h5err)

    mua_chk : if (mua == 1) then
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)  !output -> file_id
        call h5gcreate_f(file_id, grp_name, grp_id, h5err) !output -> grp_id 

    else mua_chk
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, h5err) !output -> fileid
        call h5gopen_f(file_id, grp_name, grp_id, h5err) !output -> grp_id

    end if mua_chk

            call h5screate_simple_f(rank1, dims1, dspace_id1, h5err) !output -> dspace_id1
            call h5screate_simple_f(rank2, dims2, dspace_id2, h5err) !output -> dspace_id2

            call h5dcreate_f(grp_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id1, &
                             dset_id1, h5err) !output -> dset_id1
            call h5dcreate_f(grp_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id2, &
                             dset_id2, h5err) !output -> dset_id2

            call h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, Mat_trans, dims1, h5err)
            call h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, force_mua, dims2, h5err)

            call h5dclose_f(dset_id2, h5err)
            call h5dclose_f(dset_id1, h5err)
            call h5sclose_f(dspace_id2, h5err)
            call h5sclose_f(dspace_id1, h5err)

        call h5gclose_f(grp_id, h5err)
    call h5fclose_f(file_id, h5err)
    call h5close_f(h5err)

end subroutine wFD_part


subroutine wFC3_Reconst(filename, FC3_Reconst, atmNum, atm_Indx)

    implicit none

    integer, parameter                                      :: rank1 = 1, rank3 = 3, &
                                                             & rank5 = 5

    character(len=*) , intent(in)                           :: filename
    real(dp), dimension(:,:,:,:,:,:), intent(in)            :: FC3_Reconst
    integer, dimension(:), intent(in)                       :: atmNum
    integer, dimension(:,:,:), intent(in)                   :: atm_Indx

    !................................... Local variable .....................................!

    character(len=100)                                      :: dsetname1='NumAtoms', &
                                                             & dsetname2, &
                                                             & dsetname3='atm_Indx', &
                                                             & fc_grp_name='FC3rd'
    
    integer(hid_t)                                          :: file_id, grp_id, &
                                                             & dspace_id1, dset_id1, &
                                                             & dspace_id2, dset_id2, &
                                                             & dspace_id3, dset_id3

    integer(hsize_t), dimension(1:rank1)                    :: dims1
    integer(hsize_t), dimension(1:rank3)                    :: dims3
    integer(hsize_t), dimension(1:rank5)                    :: dims5

    integer                                                 :: mu, Nbasis
    character(len=24)                                       :: mu_char
    character(len=8)                                        :: frmt

    integer(hsize_t), dimension(1:6)                        :: dims6
    integer                                                 :: h5err

    !................................... Local variable .....................................!

    Nbasis = size( atmNum )

    dims1 = size( atmNum )
    dims3 = shape( atm_Indx )
    dims6 = shape( FC3_Reconst )
    dims5 = dims6(1:5)

    img1_chk: if ( this_image() == 1 ) then
        write(*, 225) filename
        225 FORMAT('Writing FC3 in file: ', A64)
    end if img1_chk

    frmt = '(I2)'

    call h5open_f(h5err)

        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)  !output -> file_id


            call h5screate_simple_f(rank1, dims1, dspace_id1, h5err) !output -> dspace_id1
            call h5screate_simple_f(rank3, dims3, dspace_id3, h5err) !output -> dspace_id3

            call h5dcreate_f(file_id, dsetname1, H5T_NATIVE_INTEGER, dspace_id1, &
                             dset_id1, h5err) !output -> dset_id1
            call h5dcreate_f(file_id, dsetname3, H5T_NATIVE_INTEGER, dspace_id3, &
                             dset_id3, h5err) !output -> dset_id3

            call h5dwrite_f(dset_id1, H5T_NATIVE_INTEGER, atmNum, dims1, h5err)
            call h5dwrite_f(dset_id3, H5T_NATIVE_INTEGER, atm_Indx, dims3, h5err)

            !================================= FC3 write ===================================!
            call h5gcreate_f(file_id, fc_grp_name, grp_id, h5err) !output -> grp_id

                mu_loop: do mu = 1, Nbasis

                    write(mu_char, frmt) mu
                    dsetname2 = 'FC3_mu'//trim(adjustl(adjustr(mu_char)))

                    call h5screate_simple_f(rank5, dims5, dspace_id2, h5err) !output -> dspace_id2
                    call h5dcreate_f(grp_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id2, &
                                     dset_id2, h5err) !output -> dset_id2

                    call h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, &
                                  & FC3_Reconst(:,:,:,:,:,mu) , dims5, h5err)

                    call h5dclose_f(dset_id2, h5err)
                    call h5sclose_f(dspace_id2, h5err)

                end do mu_loop

            call h5gclose_f(grp_id, h5err)
            !================================= FC3 write ===================================!

            call h5dclose_f(dset_id3, h5err)
            call h5dclose_f(dset_id1, h5err)

            call h5sclose_f(dspace_id3, h5err)
            call h5sclose_f(dspace_id1, h5err)

        call h5fclose_f(file_id, h5err)

    call h5close_f(h5err)

end subroutine wFC3_Reconst


    subroutine ShengBTE_FC3rd(sys, Temp, FC3_Reconst, atmNum, atm_Indx)
    
        implicit none
    
        integer, parameter                                      :: order = 3

        type(cell), intent(in)                                  :: sys
        real(dp), intent(in)                                    :: Temp
        real(dp), dimension(:,:,:,:,:,:), intent(in)            :: FC3_Reconst
        integer, dimension(:), intent(in)                       :: atmNum
        integer, dimension(:,:,:), intent(in)                   :: atm_Indx
    
        !................................... Local variable .....................................!

        character(len=19)                                   :: filename='FORCE_CONSTANTS_3RD'

        real(dp), dimension(3, order-1)                     :: Lat_xyz
        real(dp)                                            :: fc

        integer                                             :: Nbasis, Natm_tot, Natm_mu, &
                                                             & atom_count
        integer, dimension(3)                               :: Lattice2, Lattice3
        integer                                             :: basis_indx1, basis_indx2, basis_indx3
        integer                                             :: alpha, beta, gama
        integer                                             :: bb, atm2, atm3, ii

        integer                                             :: ui, err
        character(len=512)                                  :: err_msg
        !................................... Local variable .....................................!

        Nbasis = size( atmNum )

        Natm_tot = 0
        basis_loop1: do bb = 1, Nbasis
            Natm_tot = Natm_tot + ( atmNum( bb ) ** (order-1) )
        end do basis_loop1
    
        write(*, 5) filename, Temp
        5 FORMAT("Writing Third-Order IFCs in ShengBTE format: ", A20, "  .Temperature: ", F8.2)

        ui = 5
        open(unit=ui, file=filename, status='REPLACE', action='WRITE', &
           & iostat=err, iomsg=err_msg)

            open_chk: if ( err /= 0 ) then
                 write(*, *) 'File OPEN failed: iostat = ', err
                 write(*, *) 'Error message = ', err_msg

            else open_chk

                !*! write(unit=ui, fmt=10) sys%prefix, Temp
                !*! 10 FORMAT("Fourth-Order IFCs of ", A12, ". Temperature : ", F6.1)

                write(unit=ui, fmt=20) Natm_tot
                20 FORMAT(I8)

                atom_count = 0
                basis_loop2: do bb = 1, Nbasis

                    basis_indx1 = bb
                    Natm_mu = atmNum( bb )

                    atom2: do atm2 = 1, Natm_mu

                        Lattice2 = atm_Indx(1:3, atm2, bb)
                        basis_indx2 = atm_Indx(4, atm2, bb)

                        Lat_xyz(:, 1) = matmul(sys%latvec, Lattice2)

                        atom3: do atm3 = 1, Natm_mu

                            Lattice3 = atm_Indx(1:3, atm3, bb)
                            basis_indx3 = atm_Indx(4, atm3, bb)

                            Lat_xyz(:, 2) = matmul(sys%latvec, Lattice3)

                            write(unit=ui, fmt=33)
                            33 FORMAT(" ")

                            atom_count = atom_count + 1
                            write(unit=ui, fmt=30) atom_count
                            30 FORMAT(I8)

                            write_latt: do ii = 1, (order-1)

                                write(unit=ui, fmt=40) Lat_xyz(:, ii)
                                !40 FORMAT(2(F17.11, ' '), F17.11)
                                40 FORMAT(2(E20.12, ' '), E20.12)

                            end do write_latt

                            write(unit=ui, fmt=50) basis_indx1, basis_indx2, basis_indx3
                            50 FORMAT( 3I6 )

                            cart1: do alpha = 1, 3
                                cart2: do beta = 1, 3
                                    cart3: do gama = 1, 3

                                        fc = FC3_Reconst(gama, beta, alpha, atm3, atm2, bb)

                                        write(unit=ui, fmt=60) alpha, beta, gama, fc
                                        60 FORMAT(2(I2, ' '), I2, '     ', E20.12)

                                    end do cart3
                                end do cart2
                            end do cart1
                                
                        end do atom3
                    end do atom2

                end do basis_loop2

                check: if ( atom_count /= Natm_tot ) then

                    write(*, 14)
                    14 FORMAT("ERROR: atom_count /= Natm_tot ")

                end if check

            end if open_chk

        close(unit=ui, status='KEEP', iostat=err, iomsg=err_msg)

        close_chk: if ( err /= 0 ) then
             write(*, *) 'File close failed: iostat = ', err
             write(*, *) 'Error message = ', err_msg
        end if close_chk
    
    end subroutine ShengBTE_FC3rd

end module hdf5_wrap

