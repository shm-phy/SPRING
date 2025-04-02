
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

module ConserveQuasiMomentum

    use kinds,          only : dp
    use constants,      only : cmp_prec, Sctr3THz2, &
                             & zero_prec, OneSixth
    use unit_cell,      only : cell
    use timer_class,    only : timer
    use hdf5,           only : hid_t, hsize_t, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER,&
                             & H5F_ACC_TRUNC_F, h5open_f, h5fcreate_f, &
                             & h5screate_simple_f, h5dcreate_f, h5dwrite_f, &
                             & h5dclose_f, h5sclose_f, h5fclose_f, h5close_f,&
                             & H5F_ACC_RDONLY_F, h5dget_type_f, h5tclose_f, &
                             & h5dget_space_f, h5sget_simple_extent_dims_f, h5dread_f, &
                             & h5dopen_f, h5fopen_f, h5gcreate_f, h5gopen_f, h5gclose_f

    use Irr_q_point,    only : q_points_data
    use FC3_mod,        only : FC3type
    use phonon_m,       only : Phon

    implicit none

    EXTERNAL                :: dlasrt

    private

    ! qConsrvEl is a datatype that stores the elemental information of the      ! 
    ! quasi-momentum conservation. The q2indx stores the q2 of q+q1-q2=G        !
    ! (q-q1-q2=G). The flag stores a logical value whether it can be obtained   !
    ! from permutaton symmetery between q1 and q2. If flag is false it can be   !
    ! obtained from permutaion symmetry, i.e, -q2(q2) and -q1(q1) already       !
    ! exists (flag is true for that) and we can find the scattering matrix      !
    ! element from them.                                                        !
    type, public            :: qConsrvEl

        logical                 :: flag
        integer                 :: q2indx
        integer                 :: SctrElPos

    end type qConsrvEl

    ! singleQ is a datatype that stores the quasi-momentum conservation information     !
    ! from a q of q+q1-q2=G (q-q1-q2=G) (here q is in irreducible BZ). The size of      !
    ! the allocatable array (with elements being qConsrvEl) is the number  of possible  !
    ! q1, i.e, the mesh size. The consrv_cnt stores the number of scattering matrix     !
    ! elements to be calculated and rest ( size of mesh - consrv_cnt ) can be obtained  !
    ! through permutation symmetry between q1 and q2.                                   !
    type, public            :: singleQ

        type(qConsrvEl), dimension(:), allocatable      :: ConsrvQ
        integer                                         :: consrv_cnt

    end type singleQ

    ! ********************        Scattering Matrix Elements       ******************** !
    type, public            :: Lmb

        real(dp), dimension(:,:,:), allocatable         :: ss1s2

    end type Lmb

    type, public            :: q1q2Lmb

        type(Lmb), dimension(:), allocatable            :: q1q2

    end type q1q2Lmb
    ! ********************        Scattering Matrix Elements       ******************** !

    ! AllQ is a datatype that stores the quasi-momentum conservationn information       !
    ! for all the q ( in irreducible BZ). ConsrvPlus stores quasi-momentum conservation !
    ! iformation for q+q1-q2=G and similarly ConsrvMinum stores quasi-momentum          !
    ! conservation information for q-q1-q2=G. The size of the allocatable array is the  !
    ! number of irreducible q, but it is devided between the images.                    !
    type, public            :: AllQ

        type(singleQ), dimension(:), allocatable      :: Consrvqq1q2

        type(singleQ), dimension(:), allocatable      :: ConsrvPlus
        type(singleQ), dimension(:), allocatable      :: ConsrvMinus

        type(q1q2Lmb), dimension(:), allocatable        :: WPlus        
        type(q1q2Lmb), dimension(:), allocatable        :: WMinus        

        contains
            procedure, public, pass                     :: QconsrvMomentum, WriteRestartFile, &
                                                         & TetrahedronWplus, TetrahedronWminus, &
                                                         & TetrahedronIterPlus, TetrahedronIterMinus

            procedure, private, pass                    :: ReadRestartFile

            procedure, private, nopass                  :: ScatterMatEl, ScatterProbability, &
                                                         & my_sort

    end type AllQ

contains


    subroutine QconsrvMomentum(this, time_passed, time_limit, Qpoints, mesh, my_Qsize, my_offset, & 
                             & sys, FC3, phonon_dat, RestartDir, ReadRestart, progress, STOP_FLAG)

        implicit none

        class(AllQ)                                                 :: this

        type(timer)                                                 :: time_passed
        real(dp), intent(in)                                        :: time_limit
        type(q_points_data), intent(in)                             :: Qpoints

        integer, dimension(3), intent(in)                           :: mesh
        integer, intent(in)                                         :: my_Qsize, my_offset

        type(cell), intent(in)                                      :: sys
        type(FC3type), intent(in)                                   :: FC3
        type(Phon), intent(in)                                      :: phonon_dat

        character(len=*), intent(in)                                :: RestartDir
        logical, intent(in)                                         :: ReadRestart
        integer, dimension(2), intent(out)                          :: progress !* Restart *!
        logical, intent(out)                                        :: STOP_FLAG

        !================================= Local variable ===================================!

        complex(dp), dimension(:, :), allocatable                   :: Evecq, Evecq1, Evecq2
        complex(dp), dimension(:, :, :), allocatable                :: Phiss1s2

        real(dp)                                                    :: elps_t, multiply_factor
        real(dp), dimension(3)                                      :: q1_float, q2_float
        real(dp), dimension(:), allocatable                         :: OnebyOmega_q, OnebyOmega_q1, OnebyOmega_q2
        real(dp), dimension(:), allocatable                         :: nBE_q1, nBE_q2
        real(dp), dimension(:,:,:), allocatable                     :: W_ss1s2_p, W_ss1s2_m

        integer, allocatable, dimension(:)                          :: sgn_mul
        integer, dimension(3)                                       :: mesh_mul, q1, q2, q, q2G, &
                                                                     & bound
        integer, allocatable, dimension(:, :)                       :: recp_latt

        integer                                                     :: num_grid_pnt, q2_unq, qindx
        integer                                                     :: mq1indx, mq2indx, PRINT_IMG
        integer                                                     :: i, i1, i2, iG, &
                                                                     & cnt, cnt_q1, num_cell
        integer                                                     :: Nbasis, qindx2, NumQ1

        integer, dimension(2)                                       :: old_progress !* Restart *!

        logical                                                     :: exist1, exist2, &
                                                                     & chk, old_chk
        logical, dimension(3)                                       :: q0chk

        !================================= Local variable ===================================!

        STOP_FLAG = .false.

        PRINT_IMG = ( 1 + num_images() ) / 2
        Nbasis = size( FC3%atmNum )

        allocate( this%WPlus( my_Qsize ) )
        allocate( this%WMinus( my_Qsize ) )

        allocate( Evecq(3*Nbasis, 3*Nbasis) )
        allocate( Evecq1(3*Nbasis, 3*Nbasis) )
        allocate( Evecq2(3*Nbasis, 3*Nbasis) )

        allocate( OnebyOmega_q(3*Nbasis) )
        allocate( OnebyOmega_q1(3*Nbasis) )
        allocate( OnebyOmega_q2(3*Nbasis) )

        allocate( nBE_q1(3*Nbasis) )
        allocate( nBE_q2(3*Nbasis) )

        allocate( Phiss1s2(3*Nbasis, 3*Nbasis, 3*Nbasis) )
        
        allocate( W_ss1s2_p(3*Nbasis, 3*Nbasis, 3*Nbasis) )
        allocate( W_ss1s2_m(3*Nbasis, 3*Nbasis, 3*Nbasis) )


        ! First allocate the ConsrvPlus with the number     !
        ! of irreducible qs the current image is given.     !
        allocate( this%Consrvqq1q2( my_Qsize ) )
        allocate( this%ConsrvPlus( my_Qsize ) )
        allocate( this%ConsrvMinus( my_Qsize ) )

        num_grid_pnt = product( mesh )
        NumQ1 = ( num_grid_pnt / 2 ) + 1

        multiply_factor = ( Sctr3THz2 / dble(num_grid_pnt) )

        !** Restart **!
        ReadRestartChk1: if ( ReadRestart ) then

            call this%ReadRestartFile( Nbasis, NumQ1, RestartDir, old_progress )

        else ReadRestartChk1

            old_progress = 0

        end if ReadRestartChk1
        !** Restart **!

        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)

        !=========== Find all the possible reciprocal lattice vectors (G) ===========!
        allocate( sgn_mul(3) )
        sgn_mul = (/0, 1, -1/)
        !sgn_mul = (/0, 1, -1, 2, -2/)
        num_cell = size(sgn_mul)

        allocate( recp_latt(3, num_cell**3) )

        cnt = 1
        i3_loop: do i2 = 1, num_cell
            i2_loop: do i1 = 1, num_cell
                i1_loop: do i = 1, num_cell

                    recp_latt(:, cnt) = (/sgn_mul(i)*mesh(1), sgn_mul(i1)*mesh(2), sgn_mul(i2)*mesh(3)/)
                    cnt = cnt + 1

                end do i1_loop
            end do i2_loop
        end do i3_loop
        !=========== Find all the possible reciprocal lattice vectors (G) ===========!

        bound = mesh / 2
        cnt_q1 = 0

        irr_q_loop: do i = 1, my_Qsize

            ! Form each q allocate the ConsrvQ of singleQ datatype with the size of     !
            ! mesh-grid points, i.e., the number of q1 possible.                        !
            allocate( this%Consrvqq1q2(i)%ConsrvQ( num_grid_pnt ) )
            allocate( this%ConsrvPlus(i)%ConsrvQ( num_grid_pnt ) )
            allocate( this%ConsrvMinus(i)%ConsrvQ( num_grid_pnt ) )

            ! Initially make the flag of all qConsrvEl to .false.    ! 
            this%Consrvqq1q2(i)%ConsrvQ(:)%flag = .false.
            this%ConsrvPlus(i)%ConsrvQ(:)%flag = .false.
            this%ConsrvMinus(i)%ConsrvQ(:)%flag = .false.

            qindx = my_offset + i
            q = Qpoints%irr_q_int(1:3, qindx)

            qindx2 = Qpoints%indx_map( dot_product( modulo( q, mesh ), mesh_mul ) + 1 )

            q0chk(1) = all( q == 0 )
            Evecq = phonon_dat%Evec(:, :, qindx2)
            OnebyOmega_q = phonon_dat%OnebyOmega(:, qindx2)

            old_progress_chk1: if ( i > old_progress(1) ) then

                allocate( this%WPlus(i)%q1q2( NumQ1 ) )
                allocate( this%WMinus(i)%q1q2( NumQ1 ) )

            end if old_progress_chk1

            progress(1) = i !* Restart *!

            if ( this_image() == PRINT_IMG ) write(*, 125) i, my_Qsize
            125 FORMAT(6X, 'Calculating Scattering Matrix-elements, q-index = ', I5, '/ ', I5, ' ...')

            q1_loop: do i1 = 1, num_grid_pnt
                q1 = Qpoints%q_pnt_int(1:3, i1)

                exist1 =  this%Consrvqq1q2(i)%ConsrvQ(i1)%flag

                q2G = q + q1 

                mq1indx = Qpoints%indx_map( dot_product( modulo( -1*q1, mesh ), mesh_mul ) + 1 )

                G_loop: do iG = 1, (num_cell**3)
                    q2 = recp_latt(:, iG) - q2G

                    out_bound_chk1: if ( .not. (any(q2 > bound) .or. any(q2 < -bound)) ) then

                        q2_unq = dot_product( modulo( q2, mesh ), mesh_mul ) + 1

                        out_bound_chk2: if ( (q2_unq <= num_grid_pnt) .and. (q2_unq >= 1) ) then

                            i2 = Qpoints%indx_map( q2_unq )

                            mq2indx = Qpoints%indx_map( dot_product( modulo( -1*q2, mesh ), mesh_mul ) + 1 )

                            exist2 =  this%Consrvqq1q2(i)%ConsrvQ(i2)%flag
                            !exist2 =  this%ConsrvPlus(i)%ConsrvQ(mq2indx)%flag

                            not_exists: if ( (.not. exist1) .and. (.not. exist2) ) then

                                cnt_q1 = cnt_q1 + 1

                                this%Consrvqq1q2(i)%ConsrvQ(i1)%flag = .true.
                                this%Consrvqq1q2(i)%ConsrvQ(i1)%q2indx = i2

                                this%ConsrvPlus(i)%ConsrvQ(i1)%flag = .true.
                                this%ConsrvPlus(i)%ConsrvQ(i1)%q2indx = mq2indx

                                this%ConsrvMinus(i)%ConsrvQ(mq1indx)%flag = .true.
                                this%ConsrvMinus(i)%ConsrvQ(mq1indx)%q2indx = mq2indx

                                old_chk = ( ( i >= old_progress(1) ) .and. ( cnt_q1 > old_progress(2) ) ) .or. &
                                        & ( i > old_progress(1) )

                                old_progress_chk2: if ( old_chk ) then

                                    q1_float = Qpoints%q_pnt(:, i1)
                                    Evecq1 = phonon_dat%Evec(:, :, i1)
                                    OnebyOmega_q1 = phonon_dat%OnebyOmega(:, i1)
                                    nBE_q1 = phonon_dat%nBE(:, i1)
                                    q0chk(2) = all( q1 == 0 )

                                    q2_float = Qpoints%q_pnt(:, i2)
                                    Evecq2 = phonon_dat%Evec(:, :, i2)
                                    OnebyOmega_q2 = phonon_dat%OnebyOmega(:, i2)
                                    nBE_q2 = phonon_dat%nBE(:, i2)
                                    q0chk(3) = all( q2 == 0 )

                                    call ScatterMatEl(q1_float, q2_float, q0chk, Nbasis, &
                                                    & FC3, sys, Evecq, Evecq1, Evecq2, Phiss1s2) 

                                    call ScatterProbability( q0chk, Nbasis, multiply_factor, Phiss1s2, &
                                                           & OnebyOmega_q, OnebyOmega_q1, OnebyOmega_q2, &
                                                           & nBE_q1, nBE_q2, W_ss1s2_p, W_ss1s2_m)

                                    allocate( this%WPlus(i)%q1q2(cnt_q1)%ss1s2(3*Nbasis, 3*Nbasis, 3*Nbasis) )
                                    this%WPlus(i)%q1q2(cnt_q1)%ss1s2 = W_ss1s2_p

                                    allocate( this%WMinus(i)%q1q2(cnt_q1)%ss1s2(3*Nbasis, 3*Nbasis, 3*Nbasis) )
                                    this%WMinus(i)%q1q2(cnt_q1)%ss1s2 = W_ss1s2_m


                                end if old_progress_chk2

                                this%ConsrvPlus(i)%ConsrvQ(i1)%SctrElPos = cnt_q1
                                this%ConsrvMinus(i)%ConsrvQ(mq1indx)%SctrElPos = cnt_q1

                                
                                !*! debug !*!
                                !*! write(*, 10)
                                !*! write(*, 12) cnt_q1
                                !*! write(*, 35) q, q1, q2, (q+q1+q2) !-recp_latt(:, iG))
                                !*! write(*, 45) q, q1, (-1*q2), (q+q1-(-1*q2)) !-recp_latt(:, iG))
                                !*! write(*, 55) q, (-1*q1), (-1*q2), (q-(-1*q1)-(-1*q2)) !-recp_latt(:, iG))
                                !*! debug !*!

                                chk = all( (q1 - q2) == 0 ) ! (q1 = 0 and q2 = 0) or (q2 = q1)

                                eq_chk: if ( (.not. chk) ) then

                                    this%Consrvqq1q2(i)%ConsrvQ(i2)%flag = .false.
                                    this%Consrvqq1q2(i)%ConsrvQ(i2)%q2indx = i1

                                    this%ConsrvPlus(i)%ConsrvQ(i2)%flag = .false.
                                    this%ConsrvPlus(i)%ConsrvQ(i2)%q2indx = mq1indx

                                    this%ConsrvMinus(i)%ConsrvQ(mq2indx)%flag = .false.
                                    this%ConsrvMinus(i)%ConsrvQ(mq2indx)%q2indx = mq1indx

                                    this%ConsrvPlus(i)%ConsrvQ(i2)%SctrElPos = cnt_q1
                                    this%ConsrvMinus(i)%ConsrvQ(mq2indx)%SctrElPos = cnt_q1

                                    !*! debug !*!
                                    !*! write(*, 35) q, q2, q1, (q+q2+q1) !-recp_latt(:, iG))
                                    !*! write(*, 45) q, (-1*(-1*q2)), (-1*q1), (q+(q2)-(-1*q1)) !-recp_latt(:, iG))
                                    !*! write(*, 55) q, (-1*q2), (-1*q1), (q-(-1*q2)-(-1*q1)) !-recp_latt(:, iG))
                                    !*! debug !*!

                                !*! debug !*!
                                !*! else eq_chk
                                !*!     write(*, 25) i1, i2, i
                                !*!     25 FORMAT("i1(", I5, ") - i2(", I5, ") = 0 found for i = ", I5)
                                !*! debug !*!

                                end if eq_chk

                                progress(2) = cnt_q1    !* Restart *!

                                elps_t = ( time_passed%elapsed_time() / 60.0_dp )

                                WriteAndStop: if ( elps_t > time_limit ) then

                                    call this%WriteRestartFile( Nbasis, NumQ1, progress, RestartDir )
                                    STOP_FLAG = .true.
                                    RETURN

                                end if WriteAndStop

                                !*! debug !*!
                                !*! write(*, 10)
                                !*! write(*, *)
                                !*! debug !*!

                            end if not_exists

                            !write(*, 44) q, q1, q2, (q+q1-q2), (q+q1-q2-recp_latt(:, iG))
                            
                        end if out_bound_chk2
                        
                    end if out_bound_chk1

                end do G_loop

            end do q1_loop

            this%ConsrvPlus(i)%consrv_cnt = cnt_q1
            !*! debug !*!
            !*! write(*, 15) cnt_q1, qindx
            !*! debug !*!
            cnt_q1 = 0

            !*! write(*, *) this%ConsrvPlus(i)%ConsrvQ(:)%flag
        end do irr_q_loop

        deallocate( Evecq, Evecq1, Evecq2 )
        deallocate( OnebyOmega_q, OnebyOmega_q1, OnebyOmega_q2 )
        deallocate( nBE_q1, nBE_q2 )
        deallocate( Phiss1s2 )
        deallocate( W_ss1s2_p, W_ss1s2_m )

        !*! debug !*!
        !*! write(*, 99) cnt_q1
        !*! 99 FORMAT('Number of quasi momentum conservation (q+q1+q2 = G ) ', I8)

        !*! 35 FORMAT('[ |', 3I3, '| + |', 3I3, '| + |', 3I3, '| = |', 3I3, '| ]')
        !*! 45 FORMAT('[ |', 3I3, '| + |', 3I3, '| - |', 3I3, '| = |', 3I3, '| ]')
        !*! 55 FORMAT('[ |', 3I3, '| - |', 3I3, '| - |', 3I3, '| = |', 3I3, '| ]')

        !*! 12 FORMAT('Count = ', I6)
        !*! 10 FORMAT(60('='))
        !*! 15 FORMAT("Total Count q = ", I6, ", i = ", I6)
        !*! debug !*!
        
        deallocate( recp_latt, sgn_mul )

    end subroutine QconsrvMomentum


    subroutine ScatterMatEl(q1, q2, q0chk, Nbasis, FC3, sys, &
                          & Evecq, Evecq1, Evecq2, Phi_ss1s2)

        implicit none

        real(dp), dimension(3), intent(in)                              :: q1, q2!, q
        logical, dimension(3), intent(in)                               :: q0chk
        integer, intent(in)                                             :: Nbasis
        type(FC3type), intent(in)                                       :: FC3
        type(cell), intent(in)                                          :: sys

        complex(dp), dimension(3*Nbasis, 3*Nbasis), intent(in)          :: Evecq, Evecq1, Evecq2

        complex(dp), dimension(3*Nbasis,3*Nbasis,3*Nbasis), intent(out) :: Phi_ss1s2

        !================================= Local variables ===================================!

        complex(dp)                                 :: expnq1R1, expnq2R2, expn_mass, phi3_mass_expn, &
                                                     & wqs_mu_alpha, wq1s1_nu_beta, wq2s2_eta_gama

        real(dp), dimension(3)                      :: N1R, N2R
        real(dp)                                    :: q1R1, q2R2, mass_fac

        integer                                     :: Natm3, Ndof, &
                                                     & rowno_mu_alpha, rowno_nu_beta, rowno_eta_gama
        integer                                     :: mu, N1, nu, N2, eta, &
                                                     & alpha, beta, gama, &
                                                     & s, s1, s2, ii
        integer, dimension(3)                       :: strt

        !================================= Local variables ===================================!

        Ndof = 3*Nbasis

        Phi_ss1s2 = dcmplx(0.0_dp, 0.0_dp)

        strt(:) = 1
        do ii = 1, 3
            if ( q0chk(ii) ) strt(ii) = 4
        end do

        mu_loop: do mu = 1, Nbasis

            Natm3 = FC3%atmNum(mu)

            N1_loop: do N1 = 1, Natm3

                nu = FC3%basisTyp(N1, mu)
                N1R = FC3%atmPosR(:, N1, mu)

                q1R1 = dot_product(q1, N1R)
                !expnq1R1 = cdexp( iu * q1R1 )

                expnq1R1 = dcmplx( dcos(q1R1), dsin(q1R1) )

                N2_loop: do N2 = 1, Natm3

                    eta = FC3%basisTyp(N2, mu)
                    N2R = FC3%atmPosR(:, N2, mu)

                    q2R2 = dot_product(q2, N2R)
                    !expnq2R2 = cdexp( iu * q2R2 )

                    expnq2R2 = dcmplx( dcos(q2R2), dsin(q2R2) )

                    mass_fac = ( sys%OnebySqrtMass(mu) * sys%OnebySqrtMass(nu) * sys%OnebySqrtMass(eta) )

                    expn_mass = expnq1R1 * expnq2R2 * mass_fac

                    alpha_loop: do alpha = 1, 3
                        rowno_mu_alpha = 3*(mu-1) + (alpha-1) + 1

                        s_loop: do s = strt(1), Ndof
                            wqs_mu_alpha = Evecq(rowno_mu_alpha, s)

                            beta_loop: do beta = 1, 3
                                rowno_nu_beta = 3*(nu-1) + (beta-1) + 1

                                s1_loop: do s1 = strt(2), Ndof
                                    wq1s1_nu_beta = Evecq1(rowno_nu_beta, s1)

                                    gama_loop: do gama = 1, 3
                                        rowno_eta_gama = 3*(eta-1) + (gama-1) + 1

                                        phi3_mass_expn = FC3%fC(gama, beta, alpha, N2, N1, mu) * expn_mass

                                        s2_loop: do s2 = strt(3), Ndof
                                            wq2s2_eta_gama = Evecq2(rowno_eta_gama, s2) 
                                            
                                            Phi_ss1s2(s2, s1, s) = Phi_ss1s2(s2, s1, s) + &
                                          & phi3_mass_expn * wqs_mu_alpha * wq1s1_nu_beta * wq2s2_eta_gama

                                        end do s2_loop

                                    end do gama_loop

                                end do s1_loop

                            end do beta_loop

                        end do s_loop

                    end do alpha_loop

                end do N2_loop

            end do N1_loop

        end do mu_loop

    end subroutine ScatterMatEl


    subroutine ScatterProbability( q0chk, Nbasis, multiply_factor, Phi_ss1s2, &
                                 & OnebyOmega_q, OnebyOmega_q1, OnebyOmega_q2, &
                                 & nBE_q1, nBE_q2, W_ss1s2_p, W_ss1s2_m )

        implicit none

        logical, dimension(3), intent(in)                                   :: q0chk
        integer, intent(in)                                                 :: Nbasis
        real(dp), intent(in)                                                :: multiply_factor

        complex(dp), dimension(3*Nbasis,3*Nbasis,3*Nbasis), intent(in)      :: Phi_ss1s2

        real(dp), dimension(3*Nbasis), intent(in)                           :: OnebyOmega_q, OnebyOmega_q1, &
                                                                             & OnebyOmega_q2
        real(dp), dimension(3*Nbasis), intent(in)                           :: nBE_q1, nBE_q2

        real(dp), dimension(3*Nbasis,3*Nbasis,3*Nbasis), intent(out)        :: W_ss1s2_p, W_ss1s2_m

        !! ===================================== Local variables =================================== !!

        real(dp), dimension(2)                      :: ReIm
        real(dp)                                    :: abs2, elmnt

        integer                                     :: Ndof, ii, s, s1, s2
        integer, dimension(3)                       :: strt

        !! ===================================== Local variables =================================== !!

        Ndof = 3 * Nbasis

        strt(:) = 1
        do ii = 1, 3
            if ( q0chk(ii) ) strt(ii) = 4
        end do

        W_ss1s2_p = 0.0_dp
        W_ss1s2_m = 0.0_dp

        s_loop2: do s = strt(1), Ndof
            s1_loop2: do s1 = strt(2), Ndof
                s2_loop2: do s2 = strt(3), Ndof

                    ReIm = transfer( Phi_ss1s2(s2, s1, s), ReIm )
                    abs2 = dot_product( ReIm, ReIm )

                    elmnt = abs2 * multiply_factor * ( OnebyOmega_q(s) * OnebyOmega_q1(s1) * OnebyOmega_q2(s2) )

                    W_ss1s2_p(s2, s1, s) = elmnt * ( nBE_q1(s1) - nBE_q2(s2) )
                    W_ss1s2_m(s2, s1, s) = elmnt * ( 1.0_dp + nBE_q1(s1) + nBE_q2(s2) )

                end do s2_loop2
            end do s1_loop2
        end do s_loop2

    end subroutine ScatterProbability


    subroutine my_sort(Omega_vert, tetra_ver_q1, tetra_ver_q2, sctrElPos, directFlag)

        implicit none

        real(dp), dimension(4), intent(inout)                       :: Omega_vert
        integer, dimension(4), intent(inout)                        :: tetra_ver_q1, tetra_ver_q2, sctrElPos
        logical, dimension(4), intent(inout)                        :: directFlag

        !! ===================================== Local variables =================================== !!

        real(dp), dimension(4)                      :: Omega_vertc
        integer, dimension(4)                       :: tetra_ver_q1c, tetra_ver_q2c, sctrElPosc
        logical, dimension(4)                       :: directFlagc

        integer, dimension(1)                       :: indx
        integer, dimension(1)                       :: indx_back
        integer                                     :: info, ii, strt_pnt, ActualIndx

        logical, dimension(4)                       :: covered
        logical                                     :: back_chk

        !! ===================================== Local variables =================================== !!

        covered(:) = .false.
        back_chk = .false.

        Omega_vertc = Omega_vert
        tetra_ver_q1c = tetra_ver_q1
        tetra_ver_q2c = tetra_ver_q2
        sctrElPosc = sctrElPos
        directFlagc = directFlag

        call dlasrt( 'I', 4, Omega_vert, info )

        vert_loop: do ii = 1, 4

            strt_pnt = 0
            loop: do

                strt_pnt = mod( (strt_pnt + 1), 4)
                if (strt_pnt == 0) strt_pnt = 4

                indx = findloc( Omega_vertc(strt_pnt:4), Omega_vert(ii) )
                ActualIndx = strt_pnt + (indx(1) - 1)

                if ( .not. covered(ActualIndx) ) EXIT loop

            end do loop

            covered( ActualIndx ) = .true.

            tetra_ver_q1(ii) = tetra_ver_q1c( ActualIndx )
            tetra_ver_q2(ii) = tetra_ver_q2c( ActualIndx )
            sctrElPos(ii) = sctrElPosc( ActualIndx )
            directFlag(ii) = directFlagc( ActualIndx )

            debug: if ( back_chk ) then

                indx_back = findloc(Omega_vertc, Omega_vert(ii), back=back_chk)

                if ( ActualIndx /= indx_back(1) ) then

                    write(*, 80)
                    write(*, 90) Omega_vert

                    80 FORMAT("WARNING: {-/+ Omega(q1) + Omega(q2)} is same in two vertices of tetrahedron")
                    90 FORMAT("         The values of sorted Mq1q2s at vertices: [", 3(F12.5, ', '), F12.5, "] THz")

                end if

            end if debug

        end do vert_loop

    end subroutine my_sort


    subroutine TetrahedronWplus(this, Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                              & Ndof, num_irr_q, Inv_tauqs)

        implicit none

        class(AllQ), intent(in)                                     :: this

        type(q_points_data), intent(in)                             :: Qpoints
        type(Phon), intent(in)                                      :: phonon_dat

        integer, dimension(3), intent(in)                           :: mesh
        integer, intent(in)                                         :: my_Qsize, my_offset, Ndof, num_irr_q


        real(dp), dimension(Ndof, num_irr_q), intent(inout)         :: Inv_tauqs

        !================================= Local variable ===================================!

        integer                                             :: NumTetra

        integer                                             :: NT, qi, qindx, qindx2, &
                                                             & s, ii, s1, s2, kk, vv

        integer, dimension(3)                               :: mesh_mul, q_int
        logical, dimension(3)                               :: q0chk

        real(dp), dimension(:), allocatable                 :: Mq

        real(dp)                                            :: Mqs, g, tetra_val
        logical                                             :: enterCase

        real(dp), dimension(4)                              :: Iver

        integer, dimension(4)                               :: ver_q1, ver_q2, sctrPos, &
                                                             & ver_q1_srt, ver_q2_srt, sctrPos_srt !,ver_q1_unq

        logical, dimension(4)                               :: drctFlag, drctFlag_srt

        real(dp), dimension(:,:), allocatable               :: Mvrtq1, Mvrtq2

        real(dp), dimension(4)                              :: Mq1q2s

        !================================= Local variable ===================================!

        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)

        NumTetra = size( Qpoints%tetrahedrons, 2 )

        allocate( Mq(Ndof) )
        allocate( Mvrtq1(Ndof, 4) )
        allocate( Mvrtq2(Ndof, 4) )

        irr_q_loop: do qi = 1, my_Qsize

            qindx = my_offset + qi

            q_int = Qpoints%irr_q_int(1:3, qindx)
            qindx2 = Qpoints%indx_map( dot_product( modulo( q_int, mesh ), mesh_mul ) + 1 )

            Mq = phonon_dat%omega(:, qindx2)

            q0chk(1) = all( q_int == 0 )

            !-!s_loop: do s = 1, Ndof

            !-!    Mqs = Mq(s)

            loop_Tetra: do NT = 1, NumTetra

                !ver_q1_unq = Qpoints%tetrahedrons(:, NT)
                ver_q1 = Qpoints%tetrahedrons(:, NT)

                vert_loop1: do ii = 1, 4

                    !ver_q1(ii) = Qpoints%indx_map( ver_q1_unq(ii) )
                    ver_q1(ii) = Qpoints%indx_map( ver_q1(ii) )

                    Mvrtq1(:, ii) = phonon_dat%omega(:, ver_q1(ii))
                    sctrPos(ii) = this%ConsrvPlus(qi)%ConsrvQ( ver_q1(ii) )%SctrElPos
                    drctFlag(ii) = this%ConsrvPlus(qi)%ConsrvQ( ver_q1(ii) )%flag

                    ver_q2(ii) = this%ConsrvPlus(qi)%ConsrvQ( ver_q1(ii) )%q2indx
                    Mvrtq2(:, ii) = phonon_dat%omega(:, ver_q2(ii))

                end do vert_loop1

                s_loop: do s = 1, Ndof
                    Mqs = Mq(s)

                    s1_loop: do s1 = 1, Ndof
                        s2_loop: do s2 = 1, Ndof

                            vert_loop2: do kk = 1, 4

                                !Mq1q2s(kk) = ( -1.0_dp * Mvrtq1(s1, kk) ) + Mvrtq2(s2, kk)

                                Mq1q2s(kk) = ( Mvrtq2(s2, kk) - Mvrtq1(s1, kk) )

                            end do vert_loop2

                            ver_q1_srt = ver_q1
                            ver_q2_srt = ver_q2
                            sctrPos_srt = sctrPos
                            drctFlag_srt = drctFlag

                            call my_sort(Mq1q2s, ver_q1_srt, ver_q2_srt, sctrPos_srt, drctFlag_srt)

                            enterCase = .false.
                            tetra_val = 0.0_dp

                            tetra_coef_case: if ( (Mqs > Mq1q2s(1)) .and. (Mqs < Mq1q2s(2)) ) then

                                ! ------------------- value of r and values of I at vertices ------------------- !
                                g = ( 3.0_dp * ( (Mqs - Mq1q2s(1)) ** 2 ) ) / &
                                  & ( (Mq1q2s(2) - Mq1q2s(1)) * & 
                                  &   (Mq1q2s(3) - Mq1q2s(1)) * &
                                  &   (Mq1q2s(4) - Mq1q2s(1)) )

                                Iver(1) = ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(1) - Mq1q2s(2)) ) + &
                                          & ( (Mqs - Mq1q2s(3)) / (Mq1q2s(1) - Mq1q2s(3)) ) + &
                                          & ( (Mqs - Mq1q2s(4)) / (Mq1q2s(1) - Mq1q2s(4)) ) ) / 3.0_dp

                                do vv = 2, 4
                                    Iver(vv) = ( (Mqs - Mq1q2s(1)) / (Mq1q2s(vv) - Mq1q2s(1)) ) / 3.0_dp
                                end do
                                ! ------------------- value of r and values of I at vertices ------------------- !

                                enterCase = .true.


                            else if ( (Mqs > Mq1q2s(2)) .and. (Mqs < Mq1q2s(3)) ) then tetra_coef_case

                                ! ------------------- value of r and values of I at vertices ------------------- !
                                g = ( 3.0_dp / ( (Mq1q2s(2) - Mq1q2s(3)) * &
                                  &              (Mq1q2s(4) - Mq1q2s(1)) ) ) * &
                                  & ( ( ( (Mqs - Mq1q2s(1)) * (Mqs - Mq1q2s(3)) ) / &
                                  &       (Mq1q2s(3) - Mq1q2s(1)) ) - &
                                  &   ( ( (Mqs - Mq1q2s(2)) * (Mqs - Mq1q2s(4)) ) / &
                                  &       (Mq1q2s(2) - Mq1q2s(4)) ) )

                                Iver(1) = ( ( (Mqs - Mq1q2s(4)) / (Mq1q2s(1) - Mq1q2s(4)) ) / 3.0_dp ) + &
                                        & ( ( ( (Mqs - Mq1q2s(3)) / (Mq1q2s(1) - Mq1q2s(3)) ) * &
                                        &   ( (Mqs - Mq1q2s(1)) / (Mq1q2s(3) - Mq1q2s(1)) ) * &
                                        &   ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) ) / &
                                        &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )

                                Iver(2) = ( ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) / 3.0_dp ) + &
                                        & ( ( ( ((Mqs - Mq1q2s(4)) / (Mq1q2s(2) - Mq1q2s(4))) ** 2 ) * &
                                        &   ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) ) / &
                                        &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )

                                Iver(3) = ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) / 3.0_dp ) + &
                                        & ( ( ( ((Mqs - Mq1q2s(1)) / (Mq1q2s(3) - Mq1q2s(1))) ** 2 ) * &
                                        &   ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) ) / &
                                        &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )
                            
                                Iver(4) = ( ( (Mqs - Mq1q2s(1)) / (Mq1q2s(4) - Mq1q2s(1)) ) / 3.0_dp ) + &
                                        & ( ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(4) - Mq1q2s(2)) ) * &
                                        &   ( (Mqs - Mq1q2s(4)) / (Mq1q2s(2) - Mq1q2s(4)) ) * &
                                        &   ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) ) / &
                                        &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )
                                ! ------------------- value of r and values of I at vertices ------------------- !

                                enterCase = .true.


                            else if ( (Mqs > Mq1q2s(3)) .and. (Mqs < Mq1q2s(4)) ) then tetra_coef_case

                                ! ------------------- value of r and values of I at vertices ------------------- !
                                g = ( 3.0_dp * ((Mqs - Mq1q2s(4)) ** 2) ) / &
                                  & ( ( Mq1q2s(4) - Mq1q2s(1) ) * ( Mq1q2s(4) - Mq1q2s(2) ) * &
                                  &   ( Mq1q2s(4) - Mq1q2s(3) ) )

                                do vv = 1, 3
                                    Iver(vv) = ( (Mqs - Mq1q2s(4)) / (Mq1q2s(vv) - Mq1q2s(4)) ) / 3.0_dp
                                end do

                                Iver(4) = ( ( (Mqs - Mq1q2s(1)) / (Mq1q2s(4) - Mq1q2s(1)) ) + &
                                        &   ( (Mqs - Mq1q2s(2)) / (Mq1q2s(4) - Mq1q2s(2)) ) + &
                                        &   ( (Mqs - Mq1q2s(3)) / (Mq1q2s(4) - Mq1q2s(3)) ) ) / 3.0_dp
                                ! ------------------- value of r and values of I at vertices ------------------- !

                                enterCase = .true.

                            else tetra_coef_case

                                does_not_enters: if ( (Mqs >= Mq1q2s(1)) .and. &
                                                    & (Mqs <= Mq1q2s(4)) .and. (.not. enterCase) ) then
                                    
                                    write(*, 55) Mqs
                                    write(*, 65) Mq1q2s

                                    55 FORMAT("WARNING: Mqs = ", F12.5, "THz. It does not enters any cases.")
                                    65 FORMAT("         The values of Mq1q2s at vertices: [", 3(F12.5, ', '), F12.5, "] THz")

                                end if does_not_enters

                                enterCase = .false.

                            end if tetra_coef_case


                            chk_enterCase: if ( enterCase ) then

                                chk_sumeq1: if ( dabs(sum(Iver) - 1.0_dp) < zero_prec ) then

                                    ! ----------------------------- sum over the 4 vertices ------------------------------ !
                                    vert_loop3: do vv = 1, 4

                                        per_symm: if ( drctFlag_srt(vv) ) then

                                            tetra_val = tetra_val + ( Iver(vv) * &
                                                      & this%WPlus(qi)%q1q2( sctrPos_srt(vv) )%ss1s2(s2, s1, s) )

                                        else per_symm

                                            !*! tetra_val = tetra_val + ( Iver(vv) * &
                                            !*!           & (-1.0_dp) * this%WPlus(qi)%q1q2( sctrPos_srt(vv) )%ss1s2(s1, s2, s) )

                                            tetra_val = tetra_val - ( Iver(vv) * &
                                                      & this%WPlus(qi)%q1q2( sctrPos_srt(vv) )%ss1s2(s1, s2, s) )

                                        end if per_symm

                                    end do vert_loop3
                                    ! ----------------------------- sum over the 4 vertices ------------------------------ !
                                    tetra_val = OneSixth * g * tetra_val

                                    Inv_tauqs(s, qindx) = Inv_tauqs(s, qindx) + tetra_val ! ** !

                                else chk_sumeq1
                                    
                                    tetra_val = 0.0_dp
                                    write(*, 45) Iver

                                end if chk_sumeq1

                            else chk_enterCase

                                tetra_val = 0.0_dp

                            end if chk_enterCase

                        end do s2_loop
                    end do s1_loop
                end do s_loop

            end do loop_Tetra

        end do irr_q_loop

        45 FORMAT("WARNING: sum( Iver(:) ) .ne. 1. Values of Iver(1:4): [", 3(F12.5, ', '), F12.5, "]")

        deallocate( Mq )
        deallocate( Mvrtq1, Mvrtq2 )

    end subroutine TetrahedronWplus


    subroutine TetrahedronWminus(this, Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                               & Ndof, num_irr_q, Inv_tauqs)

        implicit none

        class(AllQ), intent(in)                                     :: this

        type(q_points_data), intent(in)                             :: Qpoints
        type(Phon), intent(in)                                      :: phonon_dat

        integer, dimension(3), intent(in)                           :: mesh
        integer, intent(in)                                         :: my_Qsize, my_offset, Ndof, num_irr_q

        real(dp), dimension(Ndof, num_irr_q), intent(inout)         :: Inv_tauqs

        !================================= Local variable ===================================!

        integer                                             :: NumTetra

        integer                                             :: NT, qi, qindx, qindx2, &
                                                             & s, ii, s1, s2, kk, vv

        integer, dimension(3)                               :: mesh_mul, q_int
        logical, dimension(3)                               :: q0chk

        real(dp), dimension(:), allocatable                 :: Mq

        real(dp)                                            :: Mqs, g, tetra_val
        logical                                             :: enterCase

        real(dp), dimension(4)                              :: Iver

        integer, dimension(4)                               :: ver_q1, ver_q2, sctrPos, &
                                                             & ver_q1_srt, ver_q2_srt, sctrPos_srt

        logical, dimension(4)                               :: drctFlag, drctFlag_srt

        real(dp), dimension(:,:), allocatable               :: Mvrtq1, Mvrtq2

        real(dp), dimension(4)                              :: Mq1q2s

        !================================= Local variable ===================================!

        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)

        NumTetra = size( Qpoints%tetrahedrons, 2 )

        allocate( Mq(Ndof) )
        allocate( Mvrtq1(Ndof, 4) )
        allocate( Mvrtq2(Ndof, 4) )

        irr_q_loop: do qi = 1, my_Qsize

            qindx = my_offset + qi

            q_int = Qpoints%irr_q_int(1:3, qindx)
            qindx2 = Qpoints%indx_map( dot_product( modulo( q_int, mesh ), mesh_mul ) + 1 )

            Mq = phonon_dat%omega(:, qindx2)

            q0chk(1) = all( q_int == 0 )

            !-!s_loop: do s = 1, Ndof

            !-!    Mqs = Mq(s)

            loop_Tetra: do NT = 1, NumTetra

                ver_q1 = Qpoints%tetrahedrons(:, NT)

                vert_loop1: do ii = 1, 4

                    ver_q1(ii) = Qpoints%indx_map( ver_q1(ii) ) ! ** !

                    Mvrtq1(:, ii) = phonon_dat%omega(:, ver_q1(ii))
                    sctrPos(ii) = this%ConsrvMinus(qi)%ConsrvQ( ver_q1(ii) )%SctrElPos
                    drctFlag(ii) = this%ConsrvMinus(qi)%ConsrvQ( ver_q1(ii) )%flag

                    ver_q2(ii) = this%ConsrvMinus(qi)%ConsrvQ( ver_q1(ii) )%q2indx
                    Mvrtq2(:, ii) = phonon_dat%omega(:, ver_q2(ii))

                end do vert_loop1

                s_loop: do s = 1, Ndof
                    Mqs = Mq(s)

                    s1_loop: do s1 = 1, Ndof
                        s2_loop: do s2 = 1, Ndof

                            vert_loop2: do kk = 1, 4

                                Mq1q2s(kk) = ( Mvrtq1(s1, kk) + Mvrtq2(s2, kk) )

                            end do vert_loop2

                            ver_q1_srt = ver_q1
                            ver_q2_srt = ver_q2
                            sctrPos_srt = sctrPos
                            drctFlag_srt = drctFlag

                            call my_sort(Mq1q2s, ver_q1_srt, ver_q2_srt, sctrPos_srt, drctFlag_srt)

                            enterCase = .false.
                            tetra_val = 0.0_dp

                            tetra_coef_case: if ( (Mqs > Mq1q2s(1)) .and. (Mqs < Mq1q2s(2)) ) then

                                ! ------------------- value of r and values of I at vertices ------------------- !
                                g = ( 3.0_dp * ( (Mqs - Mq1q2s(1)) ** 2 ) ) / &
                                  & ( (Mq1q2s(2) - Mq1q2s(1)) * & 
                                  &   (Mq1q2s(3) - Mq1q2s(1)) * &
                                  &   (Mq1q2s(4) - Mq1q2s(1)) )

                                Iver(1) = ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(1) - Mq1q2s(2)) ) + &
                                          & ( (Mqs - Mq1q2s(3)) / (Mq1q2s(1) - Mq1q2s(3)) ) + &
                                          & ( (Mqs - Mq1q2s(4)) / (Mq1q2s(1) - Mq1q2s(4)) ) ) / 3.0_dp

                                do vv = 2, 4
                                    Iver(vv) = ( (Mqs - Mq1q2s(1)) / (Mq1q2s(vv) - Mq1q2s(1)) ) / 3.0_dp
                                end do
                                ! ------------------- value of r and values of I at vertices ------------------- !

                                enterCase = .true.


                            else if ( (Mqs > Mq1q2s(2)) .and. (Mqs < Mq1q2s(3)) ) then tetra_coef_case

                                ! ------------------- value of r and values of I at vertices ------------------- !
                                g = ( 3.0_dp / ( (Mq1q2s(2) - Mq1q2s(3)) * &
                                  &              (Mq1q2s(4) - Mq1q2s(1)) ) ) * &
                                  & ( ( ( (Mqs - Mq1q2s(1)) * (Mqs - Mq1q2s(3)) ) / &
                                  &       (Mq1q2s(3) - Mq1q2s(1)) ) - &
                                  &   ( ( (Mqs - Mq1q2s(2)) * (Mqs - Mq1q2s(4)) ) / &
                                  &       (Mq1q2s(2) - Mq1q2s(4)) ) )

                                Iver(1) = ( ( (Mqs - Mq1q2s(4)) / (Mq1q2s(1) - Mq1q2s(4)) ) / 3.0_dp ) + &
                                        & ( ( ( (Mqs - Mq1q2s(3)) / (Mq1q2s(1) - Mq1q2s(3)) ) * &
                                        &   ( (Mqs - Mq1q2s(1)) / (Mq1q2s(3) - Mq1q2s(1)) ) * &
                                        &   ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) ) / &
                                        &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )

                                Iver(2) = ( ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) / 3.0_dp ) + &
                                        & ( ( ( ((Mqs - Mq1q2s(4)) / (Mq1q2s(2) - Mq1q2s(4))) ** 2 ) * &
                                        &   ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) ) / &
                                        &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )

                                Iver(3) = ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) / 3.0_dp ) + &
                                        & ( ( ( ((Mqs - Mq1q2s(1)) / (Mq1q2s(3) - Mq1q2s(1))) ** 2 ) * &
                                        &   ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) ) / &
                                        &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )
                            
                                Iver(4) = ( ( (Mqs - Mq1q2s(1)) / (Mq1q2s(4) - Mq1q2s(1)) ) / 3.0_dp ) + &
                                        & ( ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(4) - Mq1q2s(2)) ) * &
                                        &   ( (Mqs - Mq1q2s(4)) / (Mq1q2s(2) - Mq1q2s(4)) ) * &
                                        &   ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) ) / &
                                        &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )
                                ! ------------------- value of r and values of I at vertices ------------------- !

                                enterCase = .true.


                            else if ( (Mqs > Mq1q2s(3)) .and. (Mqs < Mq1q2s(4)) ) then tetra_coef_case

                                ! ------------------- value of r and values of I at vertices ------------------- !
                                g = ( 3.0_dp * ((Mqs - Mq1q2s(4)) ** 2) ) / &
                                  & ( ( Mq1q2s(4) - Mq1q2s(1) ) * ( Mq1q2s(4) - Mq1q2s(2) ) * &
                                  &   ( Mq1q2s(4) - Mq1q2s(3) ) )

                                do vv = 1, 3
                                    Iver(vv) = ( (Mqs - Mq1q2s(4)) / (Mq1q2s(vv) - Mq1q2s(4)) ) / 3.0_dp
                                end do

                                Iver(4) = ( ( (Mqs - Mq1q2s(1)) / (Mq1q2s(4) - Mq1q2s(1)) ) + &
                                        &   ( (Mqs - Mq1q2s(2)) / (Mq1q2s(4) - Mq1q2s(2)) ) + &
                                        &   ( (Mqs - Mq1q2s(3)) / (Mq1q2s(4) - Mq1q2s(3)) ) ) / 3.0_dp
                                ! ------------------- value of r and values of I at vertices ------------------- !

                                enterCase = .true.

                            else tetra_coef_case

                                does_not_enters: if ( (Mqs >= Mq1q2s(1)) .and. &
                                                    & (Mqs <= Mq1q2s(4)) .and. (.not. enterCase) ) then
                                    
                                    write(*, 55) Mqs
                                    write(*, 65) Mq1q2s

                                    55 FORMAT("WARNING: Mqs = ", F12.5, "THz. It does not enters any cases.")
                                    65 FORMAT("         The values of Mq1q2s at vertices: [", 3(F12.5, ', '), F12.5, "] THz")

                                end if does_not_enters

                                enterCase = .false.

                            end if tetra_coef_case


                            chk_enterCase: if ( enterCase ) then

                                chk_sumeq1: if ( dabs(sum(Iver) - 1.0_dp) < zero_prec ) then

                                    ! ----------------------------- sum over the 4 vertices ------------------------------ !
                                    vert_loop3: do vv = 1, 4

                                        per_symm: if ( drctFlag_srt(vv) ) then

                                            tetra_val = tetra_val + ( Iver(vv) * &
                                                      & this%WMinus(qi)%q1q2( sctrPos_srt(vv) )%ss1s2(s2, s1, s) )

                                        else per_symm

                                            tetra_val = tetra_val + ( Iver(vv) * &
                                                      & this%WMinus(qi)%q1q2( sctrPos_srt(vv) )%ss1s2(s1, s2, s) )

                                        end if per_symm

                                    end do vert_loop3
                                    ! ----------------------------- sum over the 4 vertices ------------------------------ !
                                    tetra_val = 0.5_dp * OneSixth * g * tetra_val

                                    Inv_tauqs(s, qindx) = Inv_tauqs(s, qindx) + tetra_val ! ** !

                                else chk_sumeq1
                                    
                                    tetra_val = 0.0_dp
                                    write(*, 45) Iver

                                end if chk_sumeq1

                            else chk_enterCase

                                tetra_val = 0.0_dp

                            end if chk_enterCase

                        end do s2_loop
                    end do s1_loop
                end do s_loop

            end do loop_Tetra

        end do irr_q_loop

        45 FORMAT("WARNING: sum( Iver(:) ) .ne. 1. Values of Iver(1:4): [", 3(F12.5, ', '), F12.5, "]")

        deallocate( Mq )
        deallocate( Mvrtq1, Mvrtq2 )

    end subroutine TetrahedronWminus

    include "RestartReadWrite.f90"
    include "TetrahedronDeltaSum.f90"

end module ConserveQuasiMomentum

