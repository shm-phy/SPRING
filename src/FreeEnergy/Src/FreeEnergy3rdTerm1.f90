
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

module FreeEnergy3rdTerm1

    use kinds,          only : dp
    use constants,      only : eVConv3_t1, OneThird
    use unit_cell,      only : cell

    use Irr_q_point,    only : q_points_data
    use FC3_mod,        only : FC3type
    use phonon_m,       only : Phon

    implicit none
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

    ! AllQ is a datatype that stores the quasi-momentum conservationn information       !
    ! for all the q ( in irreducible BZ). ConsrvPlus stores quasi-momentum conservation !
    ! iformation for q+q1+q2=G and similarly ConsrvMinum stores quasi-momentum          !
    ! conservation information for q-q1-q2=G. The size of the allocatable array is the  !
    ! number of irreducible q, but it is devided between the images.                    !
    type, public            :: AllQ

        type(singleQ), dimension(:), allocatable      :: Consrvqq1q2

        contains
            procedure, public, pass                     :: QconsrvMomentum 

            procedure, private, nopass                  :: ScatterMatEl, ScatterProbability, &
                                                         & FreeEnergy3rd_term1

    end type AllQ

contains


    subroutine QconsrvMomentum(this, Qpoints, mesh, my_Qsize, my_offset, & 
                             & sys, FC3, phonon_dat, FreeEng3rd_t1)

        implicit none

        class(AllQ)                                                 :: this

        type(q_points_data), intent(in)                             :: Qpoints
        integer, dimension(3), intent(in)                           :: mesh
        integer, intent(in)                                         :: my_Qsize, my_offset

        type(cell), intent(in)                                      :: sys
        type(FC3type), intent(in)                                   :: FC3
        type(Phon), intent(in)                                      :: phonon_dat

        real(dp), intent(out)                                       :: FreeEng3rd_t1

        !================================= Local variable ===================================!

        complex(dp), dimension(:, :), allocatable                   :: Evecq, Evecq1, Evecq2
        complex(dp), dimension(:, :, :), allocatable                :: Phiss1s2

        real(dp)                                                    :: omega_cutoff, FreeEng
        real(dp), dimension(3)                                      :: q1_float, q2_float
        real(dp), dimension(:), allocatable                         :: omega_q, omega_q1, omega_q2
        real(dp), dimension(:), allocatable                         :: nBE_q, nBE_q1, nBE_q2
        real(dp), dimension(:,:,:), allocatable                     :: W_ss1s2

        integer, allocatable, dimension(:)                          :: sgn_mul
        integer, dimension(3)                                       :: mesh_mul, q1, q2, q, q2G, &
                                                                     & bound
        integer, allocatable, dimension(:, :)                       :: recp_latt

        integer                                                     :: num_grid_pnt, q2_unq, qindx
        integer                                                     :: i, i0, i1, i2, iG, & !,nn
                                                                     & cnt, cnt_q1, num_cell

        integer                                                     :: Nbasis, NumQ1, degeneracy, PRINT_IMG

        logical                                                     :: exist1, exist2, &
                                                                     & chk
        logical, dimension(3)                                       :: q0chk

        !================================= Local variable ===================================!

        FreeEng3rd_t1 = 0.0_dp

        PRINT_IMG = ( 1 + num_images() ) / 2

        Nbasis = size( FC3%atmNum )

        allocate( Evecq(3*Nbasis, 3*Nbasis) )
        allocate( Evecq1(3*Nbasis, 3*Nbasis) )
        allocate( Evecq2(3*Nbasis, 3*Nbasis) )

        allocate( omega_q(3*Nbasis) )
        allocate( omega_q1(3*Nbasis) )
        allocate( omega_q2(3*Nbasis) )

        allocate( nBE_q(3*Nbasis) )
        allocate( nBE_q1(3*Nbasis) )
        allocate( nBE_q2(3*Nbasis) )

        omega_cutoff = phonon_dat%omega_min

        ! First allocate the ConsrvPlus with the number     !
        ! of irreducible qs the current image is given.     !
        allocate( this%Consrvqq1q2( my_Qsize ) )

        num_grid_pnt = product( mesh )
        NumQ1 = ( num_grid_pnt / 2 ) + 1

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

            ! Initially make the flag of all qConsrvEl to .false.    ! 
            this%Consrvqq1q2(i)%ConsrvQ(:)%flag = .false.
            !-! do nn = 1, num_grid_pnt
            !-!     this%Consrvqq1q2(i)%ConsrvQ(nn)%flag = .false.
            !-! end do

            qindx = my_offset + i
            q = Qpoints%irr_q_int(1:3, qindx)

            degeneracy = Qpoints%irr_q_int(4, qindx)

            i0 = Qpoints%indx_map( dot_product( modulo( q, mesh ), mesh_mul ) + 1 )

            q0chk(1) = all( q == 0 )
            Evecq = phonon_dat%Evec(:, :, i0)
            omega_q = phonon_dat%omega(:, i0)
            nBE_q = phonon_dat%nBE(:, i0)

            if ( this_image() == PRINT_IMG ) write(*, 125) i, my_Qsize
            125 FORMAT( 10X, 'Calculating scattering matrix-elements for free-energy 3rd (term 1), q-index = ', &
                      & I5, '/ ', I5, ' ...')

            q1_loop: do i1 = 1, num_grid_pnt

                q1 = Qpoints%q_pnt_int(1:3, i1)

                exist1 =  this%Consrvqq1q2(i)%ConsrvQ(i1)%flag

                q2G = q + q1 

                G_loop: do iG = 1, (num_cell**3)
                    q2 = recp_latt(:, iG) - q2G

                    out_bound_chk1: if ( .not. (any(q2 > bound) .or. any(q2 < -bound)) ) then

                        q2_unq = dot_product( modulo( q2, mesh ), mesh_mul ) + 1

                        out_bound_chk2: if ( (q2_unq <= num_grid_pnt) .and. (q2_unq >= 1) ) then

                            i2 = Qpoints%indx_map( q2_unq )

                            exist2 =  this%Consrvqq1q2(i)%ConsrvQ(i2)%flag

                            not_exists: if ( (.not. exist1) .and. (.not. exist2) ) then

                                cnt_q1 = cnt_q1 + 1

                                this%Consrvqq1q2(i)%ConsrvQ(i1)%flag = .true.
                                this%Consrvqq1q2(i)%ConsrvQ(i1)%q2indx = i2

                                q1_float = Qpoints%q_pnt(:, i1)
                                Evecq1 = phonon_dat%Evec(:, :, i1)
                                omega_q1 = phonon_dat%omega(:, i1)
                                nBE_q1 = phonon_dat%nBE(:, i1)
                                q0chk(2) = all( q1 == 0 )

                                q2_float = Qpoints%q_pnt(:, i2)
                                Evecq2 = phonon_dat%Evec(:, :, i2)
                                omega_q2 = phonon_dat%omega(:, i2)
                                nBE_q2 = phonon_dat%nBE(:, i2)
                                q0chk(3) = all( q2 == 0 )

                                call ScatterMatEl(q1_float, q2_float, q0chk, Nbasis, &
                                                & FC3, sys, Evecq, Evecq1, Evecq2, Phiss1s2) 

                                call ScatterProbability( Phiss1s2, q0chk, Nbasis, num_grid_pnt, &
                                                       & omega_q, omega_q1, omega_q2, W_ss1s2 )

                                chk = .not. ( all( (q1 - q2) == 0 ) ) ! (q1 = 0 and q2 = 0) or (q2 = q1)

                                call FreeEnergy3rd_term1( Nbasis, W_ss1s2, omega_cutoff, &
                                                       & omega_q, omega_q1, omega_q2, &
                                                       & nBE_q, nBE_q1, nBE_q2, chk, FreeEng )

                                FreeEng3rd_t1 = FreeEng3rd_t1 + dble( degeneracy ) * (-1.0_dp * eVConv3_t1) * FreeEng

                                this%Consrvqq1q2(i)%ConsrvQ(i1)%SctrElPos = cnt_q1

                                deallocate( Phiss1s2 )
                                deallocate( W_ss1s2 )
                                
                                !*! debug !*!
                                !*! write(*, 10)
                                !*! write(*, 12) cnt_q1
                                !*! write(*, 35) q, q1, q2, (q+q1+q2) !-recp_latt(:, iG))
                                !*! debug !*!

                                eq_chk: if ( chk ) then

                                    this%Consrvqq1q2(i)%ConsrvQ(i2)%flag = .false.
                                    this%Consrvqq1q2(i)%ConsrvQ(i2)%q2indx = i1

                                    this%Consrvqq1q2(i)%ConsrvQ(i2)%SctrElPos = cnt_q1

                                    !*! debug !*!
                                    !*! write(*, 35) q, q2, q1, (q+q2+q1) !-recp_latt(:, iG))
                                    !*! debug !*!

                                !*! debug !*!
                                !*! else eq_chk
                                !*!     write(*, 25) i1, i2, i
                                !*!     25 FORMAT("i1(", I5, ") - i2(", I5, ") = 0 found for i = ", I5)
                                !*! debug !*!

                                end if eq_chk

                                !*! debug !*!
                                !*! write(*, 10)
                                !*! write(*, *)
                                !*! debug !*!

                            end if not_exists

                        end if out_bound_chk2
                        
                    end if out_bound_chk1

                end do G_loop

            end do q1_loop

            this%Consrvqq1q2(i)%consrv_cnt = cnt_q1

            !*! debug !*!
            !*! write(*, 15) cnt_q1, qindx
            !*! debug !*!
            cnt_q1 = 0

            !*! write(*, *) this%ConsrvPlus(i)%ConsrvQ(:)%flag
        end do irr_q_loop

        deallocate( Evecq, Evecq1, Evecq2 )
        deallocate( omega_q, omega_q1, omega_q2 )
        deallocate( nBE_q1, nBE_q2 )

        !*! debug !*!
        !*! write(*, 99) cnt_q1
        !*! 99 FORMAT('Number of quasi momentum conservation (q+q1+q2 = G ) ', I8)

        !*! 35 FORMAT('[ |', 3I3, '| + |', 3I3, '| + |', 3I3, '| = |', 3I3, '| ]')

        !*! 12 FORMAT('Count = ', I6)
        !*! 10 FORMAT(60('='))
        !*! 15 FORMAT("Total Count q = ", I6, ", i = ", I6)
        !*! debug !*!
        
        deallocate( recp_latt, sgn_mul )

    end subroutine QconsrvMomentum


    subroutine ScatterMatEl(q1, q2, q0chk, Nbasis, FC3, sys, &
                          & Evecq, Evecq1, Evecq2, Phi_ss1s2)

        implicit none

        real(dp), dimension(3), intent(in)                          :: q1, q2!, q
        logical, dimension(3), intent(in)                           :: q0chk
        integer, intent(in)                                         :: Nbasis
        type(FC3type), intent(in)                                   :: FC3
        type(cell), intent(in)                                      :: sys
        complex(dp), dimension(:, :), intent(in)                    :: Evecq, Evecq1, Evecq2

        complex(dp), dimension(:,:,:), allocatable, intent(out)     :: Phi_ss1s2

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

        allocate( Phi_ss1s2(Ndof, Ndof, Ndof) )
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

                    mass_fac = 1.0_dp / dsqrt(  sys%mass(mu) *  sys%mass(nu) * sys%mass(eta) )

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


    subroutine ScatterProbability( Phi_ss1s2, q0chk, Nbasis, Nq, &
                                 & omega_q, omega_q1, omega_q2, W_ss1s2 )

        implicit none

        complex(dp), dimension(:,:,:), intent(in)                   :: Phi_ss1s2
        logical, dimension(3), intent(in)                           :: q0chk
        integer, intent(in)                                         :: Nbasis, Nq
        real(dp), dimension(:), intent(in)                          :: omega_q, omega_q1, omega_q2

        real(dp), dimension(:,:,:), allocatable, intent(out)        :: W_ss1s2

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

        allocate( W_ss1s2(Ndof, Ndof, Ndof) )

        W_ss1s2 = 0.0_dp

        s_loop2: do s = strt(1), Ndof
            s1_loop2: do s1 = strt(2), Ndof
                s2_loop2: do s2 = strt(3), Ndof

                    ReIm = transfer( Phi_ss1s2(s2, s1, s), ReIm )
                    abs2 = dot_product( ReIm, ReIm )

                    elmnt = abs2 / ( omega_q(s) * omega_q1(s1) * omega_q2(s2) * dble(Nq) )

                    W_ss1s2(s2, s1, s) = elmnt 

                end do s2_loop2
            end do s1_loop2
        end do s_loop2

    end subroutine ScatterProbability


    subroutine FreeEnergy3rd_term1( Nbasis, W_ss1s2, omega_cutoff, &
                                  & omega_q, omega_q1, omega_q2, &
                                  & nBE_q, nBE_q1, nBE_q2, per_chk, FreeEng )

        implicit none

        integer, intent(in)                                         :: Nbasis
        real(dp), dimension(:,:,:), intent(in)                      :: W_ss1s2
        real(dp), intent(in)                                        :: omega_cutoff
        real(dp), dimension(:), intent(in)                          :: omega_q, omega_q1, omega_q2
        real(dp), dimension(:), intent(in)                          :: nBE_q, nBE_q1, nBE_q2
        !logical, dimension(3), intent(in)                           :: q0chk
        logical, intent(in)                                         :: per_chk

        real(dp), intent(out)                                       :: FreeEng

        !! ===================================== Local variables =================================== !!

        real(dp)                                    :: omega_s, omega_s1, omega_s2, &
                                                     & n_s, n_s1, n_s2
        real(dp)                                    :: denom1, denom2, denom2_per, & 
                                                     & mul_term1, mul_term2, mul_term

        integer                                     :: Ndof, s, s1, s2
        
        logical                                     :: not_zero1, not_zero2, not_zero2_per

        !! ===================================== Local variables =================================== !!

        FreeEng = 0.0_dp

        Ndof = 3 * Nbasis

        s_loop: do s = 1, Ndof
            omega_s = omega_q(s)
            n_s = nBE_q(s)

            s1_loop: do s1 = 1, Ndof
                omega_s1 = omega_q1(s1)
                n_s1 = nBE_q1(s1)

                s2_loop: do s2 = 1, Ndof
                    omega_s2 = omega_q2(s2)
                    n_s2 = nBE_q2(s2)

                    mul_term1 = 0.0_dp
                    mul_term2 = 0.0_dp
                    mul_term = 0.0_dp
                    
                    denom1 = omega_s + omega_s1 + omega_s2
                    denom2 = omega_s + omega_s1 - omega_s2
                    denom2_per = omega_s + omega_s2 - omega_s1

                    not_zero1 = ( dabs( denom1 ) > omega_cutoff )
                    not_zero2 = ( dabs( denom2 ) > omega_cutoff )
                    not_zero2_per = ( dabs( denom2_per ) > omega_cutoff )

                    AvoidDivg1: if ( not_zero1 ) then
                        mul_term1 = (n_s * n_s1 + n_s + OneThird) / denom1

                        PermuteChk1: if ( per_chk ) then
                        mul_term1 = mul_term1 + (n_s * n_s2 + n_s + OneThird) / denom1

                        end if PermuteChk1

                    end if AvoidDivg1

                    AvoidDivg2: if ( not_zero2 ) then
                        mul_term2 = ( 2.0_dp * n_s * n_s2 - n_s * n_s1 + n_s2 ) / denom2

                    end if AvoidDivg2

                    PermuteChk2: if ( per_chk .and. not_zero2_per ) then
                        mul_term2 = mul_term2 + ( 2.0_dp * n_s * n_s1 - n_s * n_s2 + n_s1 ) / denom2_per

                    end if PermuteChk2

                    mul_term = mul_term1 + mul_term2

                    FreeEng = FreeEng + ( W_ss1s2(s2, s1, s) * mul_term )

                end do s2_loop
            end do s1_loop
        end do s_loop

    end subroutine FreeEnergy3rd_term1

end module FreeEnergy3rdTerm1

