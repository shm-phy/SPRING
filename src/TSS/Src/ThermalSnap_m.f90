
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

module ThermalStochasticSnap

    use kinds,              only : dp
    use constants,          only : PI, AngConv
    use unit_cell,          only : cell
    use Irr_q_point,        only : q_points_data
    use phonon_m,           only : Phon
    use mklWrap,            only : dinverse

#ifndef GNU
    use ifport,             only : SYSTEM !, SRAND, RAND ! *Old-version-compatibility*
#endif

    implicit none

    private

    public      :: GenerateInitialPos, CreateSnapshotType1, &
                 & CreateSnapshotType2

contains

    subroutine GenerateInitialPos(sys, InitPos, R_Ncell, cell_basis_record)

        implicit none

        type(cell), intent(in)                              :: sys

        real(dp), dimension(:,:), allocatable, intent(out)  :: InitPos
        real(dp), dimension(:,:), allocatable, intent(out)  :: R_Ncell
        integer, dimension(:,:), allocatable, intent(out)   :: cell_basis_record

        ! ================================= Local variables ================================= !

        integer                                         :: Nbasis, Ncell, Natm, atm_count
        integer, dimension(3)                           :: Ntrans
        real(dp), dimension(3)                          :: atm_pos_frac

        integer                                         :: nx, ny, nz, nb

        ! ================================= Local variables ================================= !

        Nbasis = sys%natm
        Ncell = product( sys%sup_cell )
        Natm = Ncell * Nbasis

        atm_count = 1

        allocate( InitPos(3, Natm) )
        allocate( cell_basis_record(4, Natm) )
        allocate( R_Ncell(3, Natm) )

        nx_loop: do nx = 0, (sys%sup_cell(1) - 1)
            ny_loop: do ny = 0, (sys%sup_cell(2) - 1)
                nz_loop: do nz = 0, (sys%sup_cell(3) - 1)

                    Ntrans = (/nx, ny, nz/)

                    bs_loop: do nb = 1, Nbasis

                        cell_basis_record(1:3, atm_count) = Ntrans
                        cell_basis_record(4, atm_count) = nb

                        atm_pos_frac = dble( Ntrans ) + sys%basis_frac(:, nb)

                        InitPos(:, atm_count) = matmul( sys%latvec, atm_pos_frac )
                        R_Ncell(:, atm_count) = matmul( sys%latvec, dble(Ntrans) )

                        atm_count = atm_count + 1

                    end do bs_loop

                end do nz_loop
            end do ny_loop
        end do nx_loop

        !*! debug !*!
        !*! write(*, *)
        !*! write(*, *) atm_count, shape(InitPos), shape(cell_basis_record)
        !*! write(*, *) "Cell-basis record: "
        !*! do nx = 1, ( atm_count - 1 )
        !*!     write(*, 45) cell_basis_record(:, nx)
        !*! end do

        !*! write(*, *)
        !*! write(*, *) "Atom positions: "
        !*! do nx = 1, ( atm_count - 1 )
        !*!     write(*, 55) InitPos(:, nx)
        !*! end do

        !*! 45 FORMAT("[", 3(I5, ', '), I5, "]")
        !*! 55 FORMAT("[", 2(F12.5, ', '), F12.5, "]")
        !*! debug !*!

    end subroutine GenerateInitialPos


    subroutine CreateSnapshotType1(sys, T, Nsnap, seed, qmesh, Qpoints, phonon_dat, &
                                 & InitPos, R_Ncell, cell_basis_record, calculator, abinit_xred, dU)

        implicit none

        complex(dp), parameter                              :: iu = dcmplx(0.0_dp, 1.0_dp)
        !-! integer, parameter                                  :: seed = 88

        type(cell), intent(in)                                      :: sys
        real(dp), intent(in)                                        :: T
        integer, intent(in)                                         :: Nsnap, seed
        integer, dimension(3)                                       :: qmesh
        type(q_points_data), intent(in)                             :: Qpoints
        type(Phon), intent(in)                                      :: phonon_dat

        real(dp), dimension(:,:), intent(in)                        :: InitPos
        real(dp), dimension(:,:), intent(in)                        :: R_Ncell
        integer, dimension(:,:), intent(in)                         :: cell_basis_record
        character(len=8), intent(in)                                :: calculator
        logical, intent(in)                                         :: abinit_xred

        real(dp), dimension(:,:,:,:,:,:), allocatable, intent(out)  :: dU

        ! ================================= Local variables ================================= !

        real(dp), dimension(:,:), allocatable       :: Udisp, SnapPos

        integer                                     :: N_irrq, degeneracy, Ndof, Natm, &
                                                     & rowno_mu_alpha, Nq, mu, err

        integer, dimension(3)                       :: q_int, CellN
        logical                                     :: q0chk
        real(dp), dimension(3)                      :: q
        real(dp), dimension(3, 3)                   :: Lat_vec_inv
        real(dp), dimension(:), allocatable         :: Mq
        real(dp), dimension(:), allocatable         :: nBEq
        complex(dp), dimension(:,:), allocatable    :: Evecq

        real(dp)                                    :: qdotR, rand1_qs, rand2_qs, &
                                                     & ev_exp, mass_fac, qs_fact, &
                                                     & NqsqrtInv, tmp, randqs_tmp

        complex(dp)                                 :: expn, wqs_mu_alpha, evExpn

        integer                                     :: ns, qi, s, strt, na, alpha, ii

        integer                                     :: seed_size
        integer, dimension(:), allocatable          :: put_seed

        ! ================================= Local variables ================================= !

        Lat_vec_inv = sys%latvec
        do ii = 1, 3
            Lat_vec_inv(:, ii) = Lat_vec_inv(:, ii) * dble(sys%sup_cell(ii))
        end do
        call dinverse( Lat_vec_inv )

        call random_seed( SIZE = seed_size )
        allocate( put_seed(seed_size) )

        N_irrq = size( Qpoints%irr_q_int, 2 )
        Nq = product( qmesh )
        Ndof = size( phonon_dat%omega, 1 )
        Natm = size( InitPos, 2 )

        NqsqrtInv = 1.0_dp / dsqrt( dble(Nq) )

        !               *** Create Directory, only unix based os works here ***               !
#ifdef GNU
        call SYSTEM('mkdir -p -v ./SnapShots/', status=err)
#else
        err = SYSTEM('mkdir -p -v ./SnapShots/')
#endif

        dir_create: if ( err /= 0 ) then
             write(*, *) 'Directory creating failed: iostat = ', err
        end if dir_create
        !               *** Create Directory, only unix based os works here ***               !

        calculator_chk0: if ( calculator == 'qe' ) then
            call WriteSnapshotQE(sys, InitPos, cell_basis_record, T, Natm, 0)

        else if ( calculator == 'siesta' ) then calculator_chk0
            call WriteSnapshotSIESTA(sys, InitPos, cell_basis_record, T, Natm, 0)

        else if ( calculator == 'abinit' ) then calculator_chk0
            call WriteSnapshotAbinit( sys, InitPos, cell_basis_record, T, Lat_vec_inv, &
                                    & abinit_xred, Natm, 0 )

        else calculator_chk0
            write(*, 97) calculator
            97 FORMAT("Unsupported Force-calculator specified (qe and siesta available) : ", A8)
            STOP

        end if calculator_chk0

        allocate( dU( 3, sys%natm, sys%sup_cell(3), sys%sup_cell(2), sys%sup_cell(1), Nsnap ) )
        dU = 0.0_dp

        allocate( Mq(Ndof) )
        allocate( nBEq(Ndof) )
        allocate( Evecq(Ndof, Ndof) )

        allocate( Udisp(3, Natm) )

        SnapShotLoop: do ns = 1, Nsnap

            !*Old-version-compatible*! call SRAND( (seed+ns) )
            !*Old-version-compatible*! randqs_tmp = RAND()

            put_seed(:) = ( seed + ns )
            call random_seed( PUT=put_seed(1:seed_size) )

            call random_number( HARVEST=randqs_tmp )
            call random_number( HARVEST=randqs_tmp )

            Udisp = 0.0_dp

            irr_q_loop: do qi = 1, N_irrq

                q_int = Qpoints%irr_q_int(1:3, qi)
                q0chk = all( q_int == 0 )
                degeneracy = Qpoints%irr_q_int(4, qi)
                q = Qpoints%irr_q(:, qi)

                Mq = phonon_dat%omega(:, qi)
                nBEq = phonon_dat%nBE(:, qi)
                Evecq = phonon_dat%Evec(:, :, qi)

                if ( q0chk ) then
                    strt = 4
                else
                    strt = 1
                end if 

                band_loop: do s = strt, Ndof

                    !*Old-version-campatible*! rand1_qs = RAND()
                    !*Old-version-campatible*! rand2_qs = RAND()

                    call random_number( HARVEST=rand1_qs )
                    call random_number( HARVEST=rand2_qs )

                    qs_fact = dsqrt ( (2.0_dp * nBEq(s) + 1.0_dp) / Mq(s) ) * &
                            & dcos( 2.0_dp * PI * rand1_qs ) * & 
                            & dsqrt( -1.0_dp * dlog( 1.0_dp - rand2_qs ) )

                    atm_loop: do na = 1, Natm

                        CellN = (cell_basis_record(1:3, na) + 1)
                        mu = cell_basis_record(4, na)

                        mass_fac = 1.0_dp / dsqrt( sys%mass(mu) )

                        qdotR = dot_product( q, R_Ncell(:, na) )
                        expn = cdexp( iu * qdotR )

                        cart_comp: do alpha = 1, 3

                            rowno_mu_alpha = 3*(mu-1) + (alpha-1) + 1
                            wqs_mu_alpha = Evecq(rowno_mu_alpha, s)

                            evExpn = wqs_mu_alpha * expn

                            ev_exp = dble(degeneracy) * dble( evExpn )

                            tmp = (AngConv * qs_fact * ev_exp * mass_fac * NqsqrtInv)

                            Udisp(alpha, na)  =  Udisp(alpha, na) + tmp

                            dU(alpha, mu, CellN(3), CellN(2), CellN(1), ns) = &
                          & dU(alpha, mu, CellN(3), CellN(2), CellN(1), ns) + tmp

                        end do cart_comp

                    end do atm_loop

                end do band_loop

            end do irr_q_loop

            SnapPos = ( InitPos + Udisp )

            calculator_chk: if ( calculator == 'qe' ) then
                call WriteSnapshotQE(sys, SnapPos, cell_basis_record, T, Natm, ns)

            else if ( calculator == 'siesta' ) then calculator_chk
                call WriteSnapshotSIESTA(sys, SnapPos, cell_basis_record, T, Natm, ns)

            else if ( calculator == 'abinit' ) then calculator_chk
                call WriteSnapshotAbinit( sys, SnapPos, cell_basis_record, T, Lat_vec_inv, &
                                        & abinit_xred, Natm, ns )

            else calculator_chk
                write(*, 87) calculator
                87 FORMAT("Unsupported Force-calculator specified (qe and siesta available) : ", A8)
                STOP

            end if calculator_chk

        end do SnapShotLoop

        !Udisp = Udisp / dsqrt( dble(Nq) )
        !*! debug !*!
        !*! write(*, *)
        !*! write(*, *) "Displacements: "
        !*! do na = 1, Natm
        !*!     write(*, 56) Udisp(:, na)
        !*! end do
        !*! 56 FORMAT("[", 2(F12.6, ', '), F12.6, "]")
        !*! debug !*!

        deallocate( Evecq )
        deallocate( Mq, nBEq )
        deallocate( Udisp, SnapPos )
        deallocate( put_seed )

    end subroutine CreateSnapshotType1


    subroutine CreateSnapshotType2(sys, T, Nsnap, seed, qmesh, Qpoints, phonon_dat, &
                                 & InitPos, R_Ncell, cell_basis_record, calculator, abinit_xred, dU)

        implicit none

        complex(dp), parameter                              :: iu = dcmplx(0.0_dp, 1.0_dp)
        !-! integer, parameter                                  :: seed = 88

        type(cell), intent(in)                                      :: sys
        real(dp), intent(in)                                        :: T
        integer, intent(in)                                         :: Nsnap, seed
        integer, dimension(3)                                       :: qmesh
        type(q_points_data), intent(in)                             :: Qpoints
        type(Phon), intent(in)                                      :: phonon_dat

        real(dp), dimension(:,:), intent(in)                        :: InitPos
        real(dp), dimension(:,:), intent(in)                        :: R_Ncell
        integer, dimension(:,:), intent(in)                         :: cell_basis_record
        character(len=8), intent(in)                                :: calculator
        logical, intent(in)                                         :: abinit_xred

        real(dp), dimension(:,:,:,:,:,:), allocatable, intent(out)  :: dU

        ! ================================= Local variables ================================= !

        real(dp), dimension(:,:), allocatable       :: Udisp, SnapPos

        integer                                     :: N_irrq, degeneracy, Ndof, Natm, &
                                                     & rowno_mu_alpha, Nq, mu, &
                                                     & new_seed, err

        integer, dimension(3)                       :: q_int, CellN
        logical                                     :: q0chk
        real(dp), dimension(3)                      :: q
        real(dp), dimension(3, 3)                   :: Lat_vec_inv
        real(dp), dimension(:), allocatable         :: Mq
        real(dp), dimension(:), allocatable         :: nBEq
        complex(dp), dimension(:,:), allocatable    :: Evecq

        real(dp)                                    :: qdotR, rand1_qs, rand2_qs, &
                                                     & ev_exp, mass_fac, qs_fact, &
                                                     & randqs_tmp, NqsqrtInv
        
        complex(dp)                                 :: expn, wqs_mu_alpha, evExpn

        integer                                     :: ns, qi, s, strt, na, alpha, ii

        integer                                     :: seed_size
        integer, dimension(:), allocatable          :: put_seed

        ! ================================= Local variables ================================= !

        Lat_vec_inv = sys%latvec
        do ii = 1, 3
            Lat_vec_inv(:, ii) = Lat_vec_inv(:, ii) * dble(sys%sup_cell(ii))
        end do
        call dinverse( Lat_vec_inv )

        call random_seed( SIZE = seed_size )
        allocate( put_seed(seed_size) )

        N_irrq = size( Qpoints%irr_q_int, 2 )
        Nq = product( qmesh )
        Ndof = size( phonon_dat%omega, 1 )
        Natm = size( InitPos, 2 )

        NqsqrtInv = 1.0_dp / dsqrt( dble(Nq) )

        allocate( dU( 3, sys%natm, sys%sup_cell(3), sys%sup_cell(2), sys%sup_cell(1), Nsnap ) )
        dU = 0.0_dp

        !               *** Create Directory, only unix based os works here ***               !
#ifdef GNU
        call SYSTEM('mkdir -p -v ./SnapShots/', status=err)
#else    
        err = SYSTEM('mkdir -p -v ./SnapShots/')
#endif

        dir_create: if ( err /= 0 ) then
             write(*, *) 'Directory creation failed: iostat = ', err
        end if dir_create
        !               *** Create Directory, only unix based os works here ***               !

        calculator_chk0: if ( calculator == 'qe' ) then
            call WriteSnapshotQE(sys, InitPos, cell_basis_record, T, Natm, 0)

        else if ( calculator == 'siesta' ) then calculator_chk0
            call WriteSnapshotSIESTA(sys, InitPos, cell_basis_record, T, Natm, 0)

        else if ( calculator == 'abinit' ) then calculator_chk0
            call WriteSnapshotAbinit( sys, InitPos, cell_basis_record, T, Lat_vec_inv, &
                                    & abinit_xred, Natm, 0 )

        else calculator_chk0
            write(*, 97) calculator
            97 FORMAT("Unsupported Force-calculator specified (qe and siesta available) : ", A8)
            STOP

        end if calculator_chk0

        allocate( Mq(Ndof) )
        allocate( nBEq(Ndof) )
        allocate( Evecq(Ndof, Ndof) )

        allocate( Udisp(3, Natm) )
        allocate( SnapPos(3, Natm) )

        SnapShotLoop: do ns = 1, Nsnap

            Udisp = 0.0_dp
            atm_loop: do na = 1, Natm

                CellN = (cell_basis_record(1:3, na) + 1)
                mu = cell_basis_record(4, na)

                mass_fac = 1.0_dp / dsqrt( sys%mass(mu) )

                cart_comp: do alpha = 1, 3

                    rowno_mu_alpha = 3*(mu-1) + (alpha-1) + 1

                    new_seed = seed + ( (Natm*3)*(ns-1) + 3*(na-1) + alpha ) !-1 + 1
                    
                    put_seed(:) = new_seed
                    call random_seed( PUT=put_seed(1:seed_size) )

                    call random_number( HARVEST=randqs_tmp )
                    call random_number( HARVEST=randqs_tmp )

                    !*Old-version-compatible*! call SRAND( new_seed )
                    !*Old-version-compatible*! randqs_tmp = RAND()
                    !-! write(*, *) "=======---------------------------------------------======="

                    irr_q_loop: do qi = 1, N_irrq

                        q_int = Qpoints%irr_q_int(1:3, qi)
                        q0chk = all( q_int == 0 )
                        degeneracy = Qpoints%irr_q_int(4, qi)
                        q = Qpoints%irr_q(:, qi)

                        Mq = phonon_dat%omega(:, qi)
                        nBEq = phonon_dat%nBE(:, qi)
                        Evecq = phonon_dat%Evec(:, :, qi)

                        qdotR = dot_product( q, R_Ncell(:, na) )
                        expn = cdexp( iu * qdotR )

                        avoid_0div: if ( q0chk ) then
                            strt = 4
                        else avoid_0div
                            strt = 1
                        end if avoid_0div

                        band_loop: do s = strt, Ndof

                            !-! randqs_tmp = RAND() ! ?
                            !*Old-version-compatible*! rand1_qs = RAND()
                            !*Old-version-compatible*! rand2_qs = RAND()
                            !-! write(*, *) rand1_qs, rand2_qs

                            call random_number( HARVEST=rand1_qs )
                            call random_number( HARVEST=rand2_qs )

                            qs_fact = dsqrt ( (2.0_dp * nBEq(s) + 1.0_dp) / Mq(s) ) * &
                                    & dcos( 2.0_dp * PI * rand1_qs ) * & 
                                    & dsqrt( -1.0_dp * dlog( 1.0_dp - rand2_qs ) )

                            wqs_mu_alpha = Evecq(rowno_mu_alpha, s)

                            evExpn = wqs_mu_alpha * expn

                            ev_exp = dble(degeneracy) * dble( evExpn )

                            Udisp(alpha, na)  =  Udisp(alpha, na) + (AngConv * &
                          & qs_fact * ev_exp * mass_fac * NqsqrtInv)

                          !_!   dU(alpha, mu, CellN(3), CellN(2), CellN(1), ns) = &
                          !_! & dU(alpha, mu, CellN(3), CellN(2), CellN(1), ns) + tmp

                        end do band_loop

                    end do irr_q_loop

                    dU(alpha, mu, CellN(3), CellN(2), CellN(1), ns) = Udisp(alpha, na)

                end do cart_comp
            end do atm_loop

            SnapPos = ( InitPos + Udisp )

            calculator_chk: if ( calculator == 'qe' ) then
                call WriteSnapshotQE(sys, SnapPos, cell_basis_record, T, Natm, ns)

            else if ( calculator == 'siesta' ) then calculator_chk
                call WriteSnapshotSIESTA(sys, SnapPos, cell_basis_record, T, Natm, ns)

            else if ( calculator == 'abinit' ) then calculator_chk
                call WriteSnapshotAbinit( sys, SnapPos, cell_basis_record, T, Lat_vec_inv, &
                                        & abinit_xred, Natm, ns )

            else calculator_chk
                write(*, 87) calculator
                87 FORMAT("Unsupported Force-calculator specified (qe and siesta available) : ", A8)
                STOP

            end if calculator_chk

        end do SnapShotLoop

        !Udisp = Udisp / dsqrt( dble(Nq) )
        !*! debug !*!
        !*! write(*, *)
        !*! write(*, *) "Displacements: "
        !*! do na = 1, Natm
        !*!     write(*, 56) Udisp(:, na)
        !*! end do
        !*! 56 FORMAT("[", 2(F12.6, ', '), F12.6, "]")
        !*! debug !*!

        deallocate( Evecq )
        deallocate( Mq, nBEq )
        deallocate( Udisp, SnapPos )
        deallocate( put_seed )

    end subroutine CreateSnapshotType2


    subroutine WriteSnapshotQE(sys, snap_shot, cell_basis_record, T, Natm, num)

        implicit none

        type(cell), intent(in)                          :: sys
        real(dp), dimension(:,:), intent(in)            :: snap_shot
        integer, dimension(:,:), intent(in)             :: cell_basis_record
        real(dp), intent(in)                            :: T
        integer, intent(in)                             :: Natm, num

        ! ================================= Local variables ================================= !

        integer                                         :: na, ii, ui, err
        character(len=3)                                :: atm_lable
        character(len=8)                                :: frmt1, frmt2
        character(len=24)                               :: Temp_char, num_char, pseudo_file
        character(len=512)                              :: out_filename, err_msg

        ! ================================= Local variables ================================= !

        frmt1 = '(F6.1)'
        write(Temp_char, frmt1) T

        frmt2 = '(I2)'
        write(num_char, frmt2) num

        out_filename = './SnapShots/'//trim(sys%prefix)//'_T'//&
                      &trim(adjustl(adjustr(Temp_char)))//'K_'//&
                      &trim(adjustl(adjustr(num_char)))//'.in'

        ui = 5
        open(unit=ui, file=out_filename, status='REPLACE', action='WRITE', &
           & iostat=err, iomsg=err_msg)
            
            open_chk: if ( err /= 0 ) then
                 write(*, *) 'File OPEN failed: iostat = ', err
                 write(*, *) 'Error message = ', err_msg

            else open_chk

                write(unit=ui, fmt=10)

                write(unit=ui, fmt=200)
                species_lbl2: do ii = 1, sys%ntyp

                    pseudo_file = trim(adjustl(adjustr(sys%typ_lbl_unq(ii))))//'.upf'
                    write(unit=ui, fmt=210) sys%typ_lbl_unq(ii), sys%mass_unq(ii), pseudo_file

                end do species_lbl2
                write(unit=ui, fmt=10)

                write(unit=ui, fmt=25)
                write(unit=ui, fmt=10)

                write(unit=ui, fmt=35)
                do ii = 1, 3

                    write(unit=ui, fmt=50) ( dble(sys%sup_cell(ii)) * &
                                           & sys%latvec(:, ii) )
                end do

                write(unit=ui, fmt=10)

                write(unit=ui, fmt=75)
                atm_loop: do na = 1, Natm

                    atm_lable = sys%typ_lbl_unq( sys%typ_indx( cell_basis_record(4, na) ) )

                    write(unit=ui, fmt=100) atm_lable, snap_shot(:, na)

                end do atm_loop

                write(unit=ui, fmt=10)

            end if open_chk

        close(unit=ui, status='KEEP', iostat=err, iomsg=err_msg)

        close_chk: if ( err /= 0 ) then
             write(*, *) 'File close failed: iostat = ', err
             write(*, *) 'Error message = ', err_msg
        end if close_chk

        !10 FORMAT('(/)')
        10 FORMAT(" ")

        200 FORMAT("ATOMIC_SPECIES")
        210 FORMAT( A3, '   ', F15.9, '   ', A12) 
        25 FORMAT("K_POINTS gamma")
        35 FORMAT("CELL_PARAMETERS angstrom")
        50 FORMAT(4X, 2(F17.11, ' '), F17.11)
        75 FORMAT("ATOMIC_POSITIONS angstrom")
        100 FORMAT( A3, ' ', 2(F17.11, ' '), F17.11 )

    end subroutine WriteSnapshotQE


    subroutine WriteSnapshotSIESTA(sys, snap_shot, cell_basis_record, T, Natm, num)

        implicit none

        type(cell), intent(in)                          :: sys
        real(dp), dimension(:,:), intent(in)            :: snap_shot
        integer, dimension(:,:), intent(in)             :: cell_basis_record
        real(dp), intent(in)                            :: T
        integer, intent(in)                             :: Natm, num

        ! ================================= Local variables ================================= !

        integer                                         :: na, nspecies, unq_indx, ii, ui, err
        character(len=4)                                :: atm_lable
        character(len=8)                                :: frmt1, frmt2
        character(len=24)                               :: Temp_char, num_char, pseudo_file
        character(len=512)                              :: out_filename, err_msg

        ! ================================= Local variables ================================= !

        frmt1 = '(F6.1)'
        write(Temp_char, frmt1) T

        frmt2 = '(I2)'
        write(num_char, frmt2) num

        out_filename = './SnapShots/'//trim(sys%prefix)//'_T'//&
                      &trim(adjustl(adjustr(Temp_char)))//'K_'//&
                      &trim(adjustl(adjustr(num_char)))//'.fdf'

        ui = 5
        open(unit=ui, file=out_filename, status='REPLACE', action='WRITE', &
           & iostat=err, iomsg=err_msg)
            
            open_chk: if ( err /= 0 ) then
                 write(*, *) 'File OPEN failed: iostat = ', err
                 write(*, *) 'Error message = ', err_msg

            else open_chk
                write(unit=ui, fmt=10)

                write(unit=ui, fmt=25) sys%ntyp
                write(unit=ui, fmt=35) Natm
                write(unit=ui, fmt=10)

                write(unit=ui, fmt=11)
                species_lbl1: do nspecies = 1, sys%ntyp

                    pseudo_file = trim(adjustl(adjustr(sys%typ_lbl_unq(nspecies))))//'.psml'
                    write(unit=ui, fmt=42) nspecies, sys%atomic_no_unq(nspecies), &
                                         & sys%typ_lbl_unq(nspecies), pseudo_file
                end do species_lbl1
                write(unit=ui, fmt=12)
                write(unit=ui, fmt=10)

                write(unit=ui, fmt=13)
                species_lbl2: do nspecies = 1, sys%ntyp

                    write(unit=ui, fmt=43) nspecies, sys%mass_unq(nspecies)

                end do species_lbl2
                write(unit=ui, fmt=14)
                write(unit=ui, fmt=10)

                write(unit=ui, fmt=19)
                write(unit=ui, fmt=10)

                write(unit=ui, fmt=15)
                do ii = 1, 3
                    write(unit=ui, fmt=50) ( dble(sys%sup_cell(ii)) * &
                                           & sys%latvec(:, ii) )
                end do
                write(unit=ui, fmt=16)
                write(unit=ui, fmt=10)

                write(unit=ui, fmt=75)
                write(unit=ui, fmt=10)

                write(unit=ui, fmt=17)
                atm_loop: do na = 1, Natm

                    unq_indx = sys%typ_indx( cell_basis_record(4, na) )
                    atm_lable = sys%typ_lbl_unq( unq_indx )

                    write(unit=ui, fmt=100) snap_shot(:, na), unq_indx, na, atm_lable 

                end do atm_loop
                write(unit=ui, fmt=18)
                write(unit=ui, fmt=10)

            end if open_chk

        close(unit=ui, status='KEEP', iostat=err, iomsg=err_msg)

        close_chk: if ( err /= 0 ) then
             write(*, *) 'File close failed: iostat = ', err
             write(*, *) 'Error message = ', err_msg
        end if close_chk

        !10 FORMAT('(/)')
        10 FORMAT(" ")

        25 FORMAT("NumberOfSpecies          ", I3)
        35 FORMAT("NumberOfAtoms            ", I3)

        11 FORMAT("%block ChemicalSpeciesLabel")
        12 FORMAT("%endblock ChemicalSpeciesLabel")
        42 FORMAT(3X, I3, '   ', I4, '   ',  A4, '   ', A24)

        13 FORMAT("%block AtomicMass")
        14 FORMAT("%endblock AtomicMass")
        43 FORMAT(3X, I3, '   ', F15.9)

        19 FORMAT("LatticeConstant           1.000000 Ang")
        15 FORMAT("%block LatticeVectors")
        16 FORMAT("%endblock LatticeVectors")
        50 FORMAT(3X, 2(F17.11, ' '), F17.11)

        75 FORMAT("AtomicCoordinatesFormat Ang")

        17 FORMAT("%block AtomicCoordinatesAndAtomicSpecies")
        18 FORMAT("%endblock AtomicCoordinatesAndAtomicSpecies")
        100 FORMAT( 4X, 2(F17.11, ' '), F17.11, '  ', I3, '  ', I6, '  ', A4 )

    end subroutine WriteSnapshotSIESTA


    subroutine WriteSnapshotAbinit( sys, snap_shot, cell_basis_record, T, Lat_vec_inv, &
                                  & abinit_xred, Natm, num )

        implicit none

        real(dp), parameter                             :: AngToBohr = 1.889726125
        integer, parameter                              :: TYPAT_LEN = 25

        type(cell), intent(in)                          :: sys
        real(dp), dimension(:,:), intent(in)            :: snap_shot
        integer, dimension(:,:), intent(in)             :: cell_basis_record
        real(dp), intent(in)                            :: T
        real(dp), dimension(3, 3), intent(in)           :: Lat_vec_inv
        integer, intent(in)                             :: Natm, num
        logical, intent(in)                             :: abinit_xred

        ! ================================= Local variables ================================= !

        integer                                         :: na, typ_atm_line_len, unq_indx, &
                                                         & nchunks, rest, ii, ui, err
        integer, dimension(Natm)                        :: typat_arr !Automatic array

        character(len=4)                                :: atm_lable
        character(len=8)                                :: frmt1, frmt2
        character(len=24)                               :: Temp_char, num_char, ntyp_char, &
                                                         & typ_atm_line_ch
        character(len=64)                               :: frmt_typat_line1, frmt_typat_line, &
                                                         & frmt_natmtyp_len
        character(len=512)                              :: out_filename, err_msg

        ! ================================= Local variables ================================= !


        frmt1 = '(F6.1)'
        write(Temp_char, frmt1) T

        frmt2 = '(I2)'
        write(num_char, frmt2) num

        out_filename = './SnapShots/'//trim(sys%prefix)//'_T'//&
                      &trim(adjustl(adjustr(Temp_char)))//'K_'//&
                      &trim(adjustl(adjustr(num_char)))//'.abi'

        ui = 5
        open(unit=ui, file=out_filename, status='REPLACE', action='WRITE', &
           & iostat=err, iomsg=err_msg)
            
            open_chk: if ( err /= 0 ) then
                 write(*, *) 'File OPEN failed: iostat = ', err
                 write(*, *) 'Error message = ', err_msg

            else open_chk

                write(unit=ui, fmt=10)

                write(unit=ui, fmt=1001)
                write(unit=ui, fmt=1000)
                write(unit=ui, fmt=1002)
                write(unit=ui, fmt=1003)
                write(unit=ui, fmt=1004)
                write(unit=ui, fmt=1005)
                write(unit=ui, fmt=10)

                write(unit=ui, fmt=25) sys%ntyp
                write(unit=ui, fmt=35) Natm

                frmt2 = '(I6)'
                write(ntyp_char, frmt2) sys%ntyp
                frmt_natmtyp_len = "('znucl        ',"//trim(adjustl(adjustr(ntyp_char)))//"(I3, '  '))"
                write(unit=ui, fmt=frmt_natmtyp_len) sys%atomic_no_unq(:)
                write(unit=ui, fmt=10)

                write(unit=ui, fmt=19)
                write(unit=ui, fmt=10)

                write(unit=ui, fmt=49) ( dble(sys%sup_cell(1)) * &
                                       & sys%latvec(:, 1) )
                do ii = 2, 3
                    write(unit=ui, fmt=50) ( dble(sys%sup_cell(ii)) * &
                                           & sys%latvec(:, ii) )
                end do
                write(unit=ui, fmt=10)

                if ( abinit_xred ) then
                    write(unit=ui, fmt=76)
                else
                    write(unit=ui, fmt=75)
                end if
                atm_loop: do na = 1, Natm

                    unq_indx = sys%typ_indx( cell_basis_record(4, na) )
                    atm_lable = sys%typ_lbl_unq( unq_indx )
                    typat_arr(na) = unq_indx

                    check_xred: if ( abinit_xred ) then
                        write(unit=ui, fmt=100) matmul( Lat_vec_inv, snap_shot(:, na) ), atm_lable, na, unq_indx

                    else check_xred
                        write(unit=ui, fmt=100) snap_shot(:, na)*AngToBohr, atm_lable, na, unq_indx
                        !write(unit=ui, fmt=100) snap_shot(:, na), atm_lable, na, unq_indx
                    end if check_xred

                end do atm_loop
                write(unit=ui, fmt=10)

                ! --------------------------------------------------------------------------------------------- !
                frmt2 = '(I6)'
                check_numatoms: if (Natm < TYPAT_LEN ) then
                    typ_atm_line_len = Natm

                else check_numatoms
                    typ_atm_line_len = TYPAT_LEN

                end if check_numatoms

                write(typ_atm_line_ch, frmt2) typ_atm_line_len
                frmt_typat_line1 = "('typat        ',"//trim(adjustl(adjustr(typ_atm_line_ch)))//"(I3, ' '))"
                write(unit=ui, fmt=frmt_typat_line1) typat_arr(1:typ_atm_line_len)

                nchunks = Natm / TYPAT_LEN
                frmt_typat_line = "(13X, "//trim(adjustl(adjustr(typ_atm_line_ch)))//"(I3, ' '))"
                chunk_loop: do ii = 2, nchunks
                    write(unit=ui, fmt=frmt_typat_line) typat_arr( ((ii-1)*TYPAT_LEN+1) : (ii*TYPAT_LEN) )
                end do chunk_loop

                rest = mod(Natm, TYPAT_LEN)
                last_line: if ( (rest /= 0) .and. (rest /= Natm) ) then

                    write(typ_atm_line_ch, frmt2) rest
                    frmt_typat_line = "(13X, "//trim(adjustl(adjustr(typ_atm_line_ch)))//"(I3, ' '))"
                    write(unit=ui, fmt=frmt_typat_line) typat_arr( (Natm-rest+1):Natm )

                end if last_line
                ! --------------------------------------------------------------------------------------------- !
                write(unit=ui, fmt=10)

                write(unit=ui, fmt=11)
                write(unit=ui, fmt=12)
                write(unit=ui, fmt=10)

            end if open_chk

        close(unit=ui, status='KEEP', iostat=err, iomsg=err_msg)

        close_chk: if ( err /= 0 ) then
             write(*, *) 'File close failed: iostat = ', err
             write(*, *) 'Error message = ', err_msg
        end if close_chk

        !10 FORMAT('(/)')
        10 FORMAT(" ")

        1001 FORMAT("chkprim      0")
        1000 FORMAT("maxnsym      3860")
        1002 FORMAT("autoparal    1")
        1003 FORMAT("#paral_kgb    1")
        1004 FORMAT("#npband       $NPROCS # Some condition must satisfy")
        1005 FORMAT("#npfft        1")

        25 FORMAT("ntypat       ", I3)
        35 FORMAT("natom        ", I3)
        19 FORMAT("acell        3*1.889726125")
        49 FORMAT("rprim        ", 2(F17.11, ' '), F17.11)
        50 FORMAT(13X, 2(F17.11, ' '), F17.11)
        75 FORMAT("xcart        ")
        76 FORMAT("xred         ")
        100 FORMAT(13X, 2(F17.11, ' '), F17.11, '  # ', A4, '  ', I6, '  ', I3 )
        11 FORMAT("kptopt       0  # 0 => read directly nkpt, kpt, kptnrm and wtk")
        12 FORMAT("nkpt         1  # Taken by default to be 0.0 0.0 0.0")

    end subroutine WriteSnapshotAbinit

end module ThermalStochasticSnap

