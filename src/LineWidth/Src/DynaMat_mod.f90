
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

module DynaMat

    use kinds,          only : dp
    use constants,      only : PI, eV_A2, wTHz, vg_ms, BEConv, &
                             & EPS, cmp_prec
    use unit_cell,      only : cell
    use FC2_mod,        only : FC2type
    use EwaldMod,       only : EwaldParam
    use EwaldSum,       only : RealSpaceEW, KSpaceEW, CombineZStar, &
                             & RealSpaceEW_gv, KSpaceEW_gv, CombineZStar_gv
    use mklWrap,        only : SolveEigen

    implicit none
    private

    public  :: get_phonon, get_phonon_singleq

contains


subroutine DynMat_SR(sys, FC2, q, Dyn)

    implicit none
    complex(dp), parameter      :: iu = dcmplx(0.0_dp, 1.0_dp)

    type(cell), intent(in)                                  :: sys
    type(FC2type), intent(in)                               :: FC2
    real(dp), dimension(3), intent(in)                      :: q
    complex(dp), dimension(:,:), allocatable, intent(out)   :: Dyn

    !! ================================= Local variables ================================= !!
    complex(dp)                             :: expnt, dm_el
    real(dp)                                :: qdotR, N2xyz(3)
    integer                                 :: N2crt(3), my_nu
    integer                                 :: basis, mu, nu, alpha, beta, &
                                               Natm, N2, col_no, row_no
    !! ================================= Local variables ================================= !!

    basis = sys%natm

    allocate(Dyn(3*basis, 3*basis))
    Dyn = dcmplx(0.0_dp, 0.0_dp)

    do nu=1, basis
        do beta=1, 3
            col_no = 3*(nu-1) + (beta-1) + 1

            do mu=1, basis
                Natm = FC2%atmNum(mu)
                do alpha=1, 3
                    dm_el = dcmplx(0.0_dp, 0.0_dp)

                    do N2=1, Natm
                       N2crt = FC2%atmIndx(1:3, N2, mu)
                       my_nu = FC2%atmIndx(4, N2, mu)
                       if (my_nu == nu) then

                           N2xyz = matmul(sys%latvec, N2crt)
                           qdotR = dot_product(q, N2xyz)
                           expnt = cdexp(iu * qdotR)
                           !expnt = dcmplx( dcos(qdotR), dsin(qdotR) )
                           dm_el = dm_el + ( FC2%fC(1, beta, alpha, N2, mu) * expnt ) 

                       end if
                    end do

                    row_no = 3*(mu-1) + (alpha-1)  + 1
                    Dyn(row_no, col_no) = dm_el / ( dsqrt(sys%mass(mu) * sys%mass(nu)) )
                end do
            end do

        end do
    end do

end subroutine DynMat_SR


subroutine DynMatgv_SR(sys, FC2, q, Dyn, dDyn_dq)

    implicit none
    complex(dp), parameter      :: iu = dcmplx(0.0_dp, 1.0_dp)

    type(cell), intent(in)                                  :: sys
    type(FC2type), intent(in)                               :: FC2
    real(dp), dimension(3), intent(in)                      :: q
    complex(dp), dimension(:,:), allocatable, intent(out)   :: Dyn
    !! ****************** Extra addition for group velocity ****************** !!
    complex(dp), dimension(:,:,:), allocatable, intent(out) :: dDyn_dq
    !! ****************** Extra addition for group velocity ****************** !!

    ! ================================ Local variables ================================ !

    integer                                 :: N2crt(3), my_nu
    real(dp)                                :: qdotR, N2xyz(3), mass_fact
    complex(dp)                             :: expnt, dm_el, fc2_expnt
    !! ****************** Extra addition for group velocity ****************** !!
    complex(dp), dimension(3)               :: ddmel_dq
    integer                                 :: kk, ii
    !! ****************** Extra addition for group velocity ****************** !!
    integer                                 :: basis, mu, nu, alpha, beta, &
                                               Natm, N2, col_no, row_no

    ! ================================ Local variables ================================ !

    basis = sys%natm

    allocate(Dyn(3*basis, 3*basis))
    Dyn = dcmplx(0.0_dp, 0.0_dp)

    !! ****************** Extra addition for group velocity ****************** !!
    allocate(dDyn_dq(3*basis, 3*basis, 3))
    dDyn_dq = dcmplx(0.0_dp, 0.0_dp)
    !! ****************** Extra addition for group velocity ****************** !!

    nu_loop: do nu=1, basis
        beta_loop: do beta=1, 3
            col_no = 3*(nu-1) + (beta-1) + 1

            mu_loop: do mu=1, basis
                Natm = FC2%atmNum(mu)
                mass_fact = 1.0_dp / dsqrt(sys%mass(mu) * sys%mass(nu))

                alpha_loop: do alpha=1, 3

                    dm_el = dcmplx(0.0_dp, 0.0_dp)
                    ddmel_dq = dcmplx(0.0_dp, 0.0_dp)

                    N2_loop: do N2=1, Natm
                       N2crt = FC2%atmIndx(1:3, N2, mu)
                       my_nu = FC2%atmIndx(4, N2, mu)

                       nu_chk: if (my_nu == nu) then

                           N2xyz = matmul(sys%latvec, N2crt)
                           qdotR = dot_product(q, N2xyz)
                           expnt = cdexp(iu * qdotR)
                           !expnt = dcmplx( dcos(qdotR), dsin(qdotR) )
                           fc2_expnt = FC2%fC(1, beta, alpha, N2, mu) * expnt

                           dm_el = dm_el + fc2_expnt

                           !! ****************** Extra addition for group velocity ****************** !!
                           drv_loop1: do kk = 1, 3

                                ddmel_dq(kk) = ddmel_dq(kk) + ( fc2_expnt * iu * N2xyz(kk) )

                           end do drv_loop1
                           !! ****************** Extra addition for group velocity ****************** !!

                       end if nu_chk
                    end do N2_loop

                    row_no = 3*(mu-1) + (alpha-1)  + 1
                    Dyn(row_no, col_no) = dm_el * mass_fact

                    !! ****************** Extra addition for group velocity ****************** !!
                    drv_loop2: do ii = 1, 3

                        dDyn_dq(row_no, col_no, ii) = ddmel_dq(ii) * mass_fact

                    end do drv_loop2
                    !! ****************** Extra addition for group velocity ****************** !!

                end do alpha_loop
            end do mu_loop

        end do beta_loop
    end do nu_loop

end subroutine DynMatgv_SR


subroutine DynMat_Ew(sys, EwaldConst, q, DynEw)

    implicit none

    type(cell), intent(in)                                      :: sys
    type(EwaldParam), intent(in)                                :: EwaldConst
    real(dp), dimension(3), intent(in)                          :: q

    complex(dp), allocatable, intent(out)                       :: DynEw(:, :)


    !! ================================= Local variables ================================= !!
    complex(dp), allocatable, dimension(:,:,:,:)    :: CRSpace, CKSpace &
                                                     &, Cq, CqZ
    integer                                         :: Nbasis, mu, nu, &
                                                     & alpha, beta
    integer                                         :: col_no, row_no
    !! ================================= Local variables ================================= !!

    Nbasis = sys%natm
    allocate(Cq(3, 3, Nbasis, Nbasis))
    Cq = dcmplx(0.0_dp, 0.0_dp)

    call RealSpaceEW(q, EwaldConst%Rpoints, sys%tau, EwaldConst%eps_inv, &
                    &EwaldConst%det, EwaldConst%Lmb, CRSpace)

    call KSpaceEW(q, EwaldConst%Gpoints, sys%tau, sys%eps, sys%vol, EwaldConst%Lmb, CKSpace)
    
    Cq = CKSpace - CRSpace - EwaldConst%CConst
    call CombineZStar(Cq, sys%Zstar, CqZ)
        
    mu_loop: do mu = 1, Nbasis
        a_loop: do alpha = 1, 3
            b_loop: do beta = 1, 3

                CqZ(beta, alpha, mu, mu) = &
                & CqZ(beta, alpha, mu, mu) - &
                & EwaldConst%CSumq0(beta, alpha, mu)
                
            end do b_loop
        end do a_loop
    end do mu_loop

    CqZ = CqZ * eV_A2

    allocate(DynEw(3*Nbasis, 3*Nbasis))
    DynEw = dcmplx(0.0_dp, 0.0_dp)

    nu_loop: do nu = 1, Nbasis
        beta_loop: do beta = 1, 3
                col_no = 3*(nu-1) + (beta-1) + 1

                mu_loop2: do mu = 1, Nbasis
                    alpha_loop: do alpha = 1, 3
                        row_no = 3*(mu-1) + (alpha-1)  + 1

                        DynEw(row_no, col_no) = &
                        & CqZ(beta, alpha, nu, mu) / ( dsqrt(sys%mass(mu) * sys%mass(nu)) )

                    end do alpha_loop
                end do mu_loop2
        end do beta_loop
    end do nu_loop

    deallocate(CRSpace)
    deallocate(CKSpace)
    deallocate(Cq)
    deallocate(CqZ)

end subroutine DynMat_Ew


subroutine DynMatgv_Ew(sys, EwaldConst, q, DynEw, dDynEw_dq)

    implicit none

    type(cell), intent(in)                                      :: sys
    type(EwaldParam), intent(in)                                :: EwaldConst
    real(dp), dimension(3), intent(in)                          :: q

    complex(dp), allocatable, intent(out)                       :: DynEw(:, :)
    !! ****************** Extra addition for group velocity ****************** !!
    complex(dp), allocatable, intent(out)                       :: dDynEw_dq(:, :, :)
    !! ****************** Extra addition for group velocity ****************** !!

    ! ================================ Local variables ================================ !

    complex(dp), allocatable, dimension(:,:,:,:)    :: CRSpace, CKSpace &
                                                     &, Cq, CqZ

    !! ****************** Extra addition for group velocity ****************** !!
    complex(dp), allocatable, dimension(:,:,:,:,:)  :: dCR_dq, dCK_dq, &
                                                     & dCq_dq, dCqZ_dq
    !! ****************** Extra addition for group velocity ****************** !!

    real(dp)                                        :: mass_fac
    integer                                         :: Nbasis, mu, nu, &
                                                     & alpha, beta, ii
    integer                                         :: col_no, row_no

    ! ================================ Local variables ================================ !

    Nbasis = sys%natm
    allocate(Cq(3, 3, Nbasis, Nbasis))
    Cq = dcmplx(0.0_dp, 0.0_dp)

    !! ****************** Extra addition for group velocity ****************** !!
    allocate(dCq_dq(3, 3, Nbasis, Nbasis, 3))
    dCq_dq = dcmplx(0.0_dp, 0.0_dp)
    !! ****************** Extra addition for group velocity ****************** !!

    call RealSpaceEW_gv(q, EwaldConst%Rpoints, sys%tau, EwaldConst%eps_inv, &
                     & EwaldConst%det, EwaldConst%Lmb, CRSpace, dCR_dq)

    call KSpaceEW_gv(q, EwaldConst%Gpoints, sys%tau, sys%eps, sys%vol, EwaldConst%Lmb, &
                   & CKSpace, dCK_dq)

    Cq = CKSpace - CRSpace - EwaldConst%CConst

    !! ****************** Extra addition for group velocity ****************** !!
    dCq_dq = dCK_dq - dCR_dq
    !! ****************** Extra addition for group velocity ****************** !!

    call CombineZStar_gv(Cq, dCq_dq, sys%Zstar, CqZ, dCqZ_dq)
        
    mu_loop: do mu = 1, Nbasis
        a_loop: do alpha = 1, 3
            b_loop: do beta = 1, 3

                CqZ(beta, alpha, mu, mu) = &
                & CqZ(beta, alpha, mu, mu) - &
                & EwaldConst%CSumq0(beta, alpha, mu)
                
            end do b_loop
        end do a_loop
    end do mu_loop

    CqZ = CqZ * eV_A2

    !! ****************** Extra addition for group velocity ****************** !!
    dCqZ_dq = dCqZ_dq * eV_A2
    !! ****************** Extra addition for group velocity ****************** !!

    allocate(DynEw(3*Nbasis, 3*Nbasis))
    DynEw = dcmplx(0.0_dp, 0.0_dp)

    !! ****************** Extra addition for group velocity ****************** !!
    allocate(dDynEw_dq(3*Nbasis, 3*Nbasis, 3))
    dDynEw_dq = dcmplx(0.0_dp, 0.0_dp)
    !! ****************** Extra addition for group velocity ****************** !!

    nu_loop: do nu = 1, Nbasis
        beta_loop: do beta = 1, 3
                col_no = 3*(nu-1) + (beta-1) + 1

                mu_loop2: do mu = 1, Nbasis
                    mass_fac = 1.0_dp / ( dsqrt(sys%mass(mu) * sys%mass(nu)) )

                    alpha_loop: do alpha = 1, 3
                        row_no = 3*(mu-1) + (alpha-1)  + 1

                        DynEw(row_no, col_no) = &
                        & CqZ(beta, alpha, nu, mu) * mass_fac
                        
                        !! ****************** Extra addition for group velocity ****************** !!
                        drv_loop: do ii = 1, 3

                            dDynEw_dq(row_no, col_no, ii) = &
                          & dCqZ_dq(beta, alpha, nu, mu, ii) * mass_fac 

                        end do drv_loop
                        !! ****************** Extra addition for group velocity ****************** !!

                    end do alpha_loop
                end do mu_loop2
        end do beta_loop
    end do nu_loop

    deallocate(CRSpace)
    deallocate(CKSpace)
    deallocate(Cq)
    deallocate(CqZ)

    !! ****************** Extra addition for group velocity ****************** !!
    deallocate(dCR_dq)
    deallocate(dCK_dq)
    deallocate(dCq_dq)
    deallocate(dCqZ_dq)
    !! ****************** Extra addition for group velocity ****************** !!

end subroutine DynMatgv_Ew


subroutine DynMat_NA(sys, EwaldConst, q, DynNA)

    implicit none

    type(cell), intent(in)                                      :: sys
    type(EwaldParam), intent(in)                                :: EwaldConst
    real(dp), dimension(3), intent(in)                          :: q

    complex(dp), allocatable, intent(out)                       :: DynNA(:, :)


    !! ====================================== Local variables ====================================== !!
    complex(dp), allocatable, dimension(:,:,:,:)    :: CRSpace, CKSpace &
                                                     &, Cq0, Cq0Z
    complex(dp), allocatable                        :: DynNA_q0(:, :)

    real(dp), dimension(3)                          :: eps_q
    real(dp), dimension(3)                          :: qzero
    real(dp)                                        :: denom, term1, term2, &
                                                     & mass_fac
    integer                                         :: Nbasis, mu, nu, &
                                                     & alpha, beta
    integer                                         :: col_no, row_no
    !! ====================================== Local variables ====================================== !!

    img1_chk1: if ( this_image() == 1) then
        write(*, 75) q
        75 FORMAT('Calculating non-analytic term at q = 0, with qhat = (', 2(E11.4, ', '), E11.4, ')')
    end if img1_chk1

    Nbasis = sys%natm
    allocate(Cq0(3, 3, Nbasis, Nbasis))
    Cq0 = dcmplx(0.0_dp, 0.0_dp)

    qzero = (/0.0_dp, 0.0_dp, 0.0_dp/)
    call RealSpaceEW(qzero, EwaldConst%Rpoints, sys%tau, EwaldConst%eps_inv, &
                    &EwaldConst%det, EwaldConst%Lmb, CRSpace)

    call KSpaceEW(qzero, EwaldConst%Gpoints, sys%tau, sys%eps, sys%vol, EwaldConst%Lmb, CKSpace)
    
    Cq0 = CKSpace - CRSpace - EwaldConst%CConst
    call CombineZStar(Cq0, sys%Zstar, Cq0Z)
        
    mu_loop: do mu = 1, Nbasis
        a_loop: do alpha = 1, 3
            b_loop: do beta = 1, 3
        
                Cq0Z(beta, alpha, mu, mu) = &
                & Cq0Z(beta, alpha, mu, mu) - &
                & EwaldConst%CSumq0(beta, alpha, mu)
        
            end do b_loop
        end do a_loop
    end do mu_loop

    Cq0Z = Cq0Z * eV_A2

    allocate(DynNA_q0(3*Nbasis, 3*Nbasis))
    DynNA_q0 = dcmplx(0.0_dp, 0.0_dp)

    allocate(DynNA(3*Nbasis, 3*Nbasis))
    DynNA = dcmplx(0.0_dp, 0.0_dp)

    nu_loop: do nu = 1, Nbasis
        beta_loop: do beta = 1, 3
                col_no = 3*(nu-1) + (beta-1) + 1
                term2 = dot_product(q, sys%Zstar(beta, :, nu))

                mu_loop2: do mu = 1, Nbasis
                    alpha_loop: do alpha = 1, 3
                        row_no = 3*(mu-1) + (alpha-1)  + 1
                        term1 = dot_product(q, sys%Zstar(alpha, :, mu))
                        mass_fac = dsqrt(sys%mass(mu) * sys%mass(nu))

                        DynNA_q0(row_no, col_no) = &
                      & Cq0Z(beta, alpha, nu, mu) / mass_fac

                        DynNA(row_no, col_no) = &
                      & (term1 * term2) / mass_fac

                    end do alpha_loop
                end do mu_loop2
        end do beta_loop
    end do nu_loop

    eps_q = matmul(q, sys%eps)
    denom = dot_product(q, eps_q)
    DynNA = (DynNA * 4.0_dp * PI * eV_A2) / (sys%vol * denom)

    DynNA = DynNA_q0 + DynNA

    deallocate(CRSpace)
    deallocate(CKSpace)
    deallocate(Cq0)
    deallocate(Cq0Z)
    deallocate(DynNA_q0)

end subroutine DynMat_NA


subroutine get_phonon(sys, FC2, EwaldConst, Nq, Ndof, myQsize, myOffset, &
                    & qpoints, T, LongEw, Evec, freq, OnebyFreq, nBE, grp_vel, &
                    & freq_min, freq_max)

    implicit none

    type(cell), intent(in)                                      :: sys
    type(FC2type), intent(in)                                   :: FC2
    type(EwaldParam), intent(in)                                :: EwaldConst

    integer, intent(in)                                         :: Nq, Ndof, myQsize, myOffset

    real(dp), dimension(3, Nq), intent(in)                      :: qpoints
    real(dp), intent(in)                                        :: T
    logical, intent(in)                                         :: LongEW

    complex(dp), dimension(Ndof, Ndof, Nq), intent(inout)       :: Evec
    real(dp), dimension(Ndof, Nq), intent(inout)                :: freq, OnebyFreq
    real(dp), dimension(Ndof, Nq), intent(inout)                :: nBE
    real(dp), dimension(3, Ndof, Nq), intent(inout)             :: grp_vel
    real(dp), intent(out)                                       :: freq_min, freq_max

    !! ====================================== Local variables ====================================== !!

    complex(dp), dimension(:,:), allocatable                :: Dyn
    complex(dp), dimension(:,:), allocatable                :: DynEw

    !! ****************** Extra addition for group velocity ****************** !!
    complex(dp), dimension(:,:), allocatable                :: DynEwhp
    complex(dp), dimension(:,:), allocatable                :: DynEwhm

    complex(dp), dimension(:), allocatable                  :: eg_vec
    complex(dp), dimension(:), allocatable                  :: dDyndq_ev
    complex(dp), dimension(3)                               :: cmplx_vel

    complex(dp), dimension(:,:,:), allocatable              :: dDyn_dq
    complex(dp), dimension(:,:,:), allocatable              :: dDynEw_dq

    real(dp), dimension(3)                                  :: q_h, dq
    real(dp)                                                :: h, omega
    integer                                                 :: kk, df, strt
    !! ****************** Extra addition for group velocity ****************** !!
    real(dp), dimension(3)                                  :: q, qhat
    real(dp), dimension(:), allocatable                     :: W

    integer                                                 :: ii, jj, my_i

    logical, dimension(:), allocatable                      :: neg_check

    !! ====================================== Local variables ====================================== !!


    allocate(neg_check(Ndof))

    allocate( eg_vec(Ndof) )
    allocate( dDyndq_ev(Ndof) )
    eg_vec = dcmplx(0.0_dp, 0.0_dp)
    dDyndq_ev = dcmplx(0.0_dp, 0.0_dp)

    !**! ToDo : Add OpenMP directives

    q_loop: do my_i = 1, myQsize

        ii = myOffset + my_i

        q = qpoints(:, ii)

        call DynMatgv_SR(sys, FC2, q, Dyn, dDyn_dq)

        Ewald_chk: if ( LongEw ) then

            q0_chk1: if ( all(dabs(q) < cmp_prec) ) then

                endPnt_chk: if ( ii /= Nq ) then
                    qhat = qpoints(:, ii+1) - qpoints(:, ii)
    
                else endPnt_chk
                    qhat = qpoints(:, ii) - qpoints(:, ii-1)

                end if endPnt_chk
    
                h = norm2( qhat ) / 8.0_dp

                call DynMat_NA(sys, EwaldConst, qhat, DynEw)

                !! ********** Switch to finite difference method for dDynEw_dq at q = 0 ********** !!

                allocate( dDynEw_dq(Ndof, Ndof, 3) )
                dDynEw_dq = dcmplx(0.0_dp, 0.0_dp)
                
                drv_loop1: do kk = 1, 3

                    dq = 0.0_dp
                    dq(kk) = 0.5_dp*h

                    q_h = q + dq
                    call DynMat_Ew(sys, EwaldConst, q_h, DynEwhp)

                    q_h = q - dq
                    call DynMat_Ew(sys, EwaldConst, q_h, DynEwhm)

                    dDynEw_dq(:, :, kk) = (DynEwhp - DynEwhm) / h

                    deallocate( DynEwhp, DynEwhm )

                end do drv_loop1

                !! ********** Switch to finite difference method for dDynEw_dq at q = 0 ********** !!
    
            else q0_chk1

                qhat = q
                call DynMatgv_Ew(sys, EwaldConst, qhat, DynEw, dDynEw_dq)
    
            end if q0_chk1
    
            Dyn = Dyn + DynEW
            !! ****************** Extra addition for group velocity ****************** !!
            dDyn_dq = dDyn_dq + dDynEw_dq
            !! ****************** Extra addition for group velocity ****************** !!
            
        end if Ewald_chk

        call SolveEigen(Dyn, W)

        neg_check = ( (dabs(W) > EPS) .and. (W < 0.0_dp) )
        W = dsqrt(dabs(W)) !* wTHz

        do jj = 1, size(neg_check)

            warning: if ( neg_check(jj) ) then 
                W(jj) = W(jj) * (-1.0_dp)

                write(*, 54) q, (W(jj) * wTHz)
                54 FORMAT("WARNING: Imaginary frequency at q = (", &
                         & 2(F9.4, ', '), F9.4, "). Freqency (omega) = ", F9.4, "THz")
            end if warning

        end do

        !! ****************** Extra addition for group velocity ****************** !!
        cmplx_vel = dcmplx(0.0_dp, 0.0_dp)

        q0_chk2: if ( all(dabs(q) < cmp_prec) ) then

            !! *************** Record maximum and minimum frequency for cut-off purpose *************** !!
            freq_max = maxval( W ) * wTHz

            freq_min = W(3) * wTHz !! Manual says it is ascending order. So it may be safe
            neg_chk: if ( freq_min < 0.0_dp ) then
                freq_min = dabs( freq_min )
            end if neg_chk

            freq_min = freq_min + 0.05_dp * freq_min
            !! *************** Record maximum and minimum frequency for cut-off purpose *************** !!

            grp_vel(:, 1, ii) = 0.0_dp
            grp_vel(:, 2, ii) = 0.0_dp
            grp_vel(:, 3, ii) = 0.0_dp

            nBE(1, ii) = 0.0_dp
            nBE(2, ii) = 0.0_dp
            nBE(3, ii) = 0.0_dp

            OnebyFreq(1, ii) = 0.0_dp !**_**!   
            OnebyFreq(2, ii) = 0.0_dp !**_**!
            OnebyFreq(3, ii) = 0.0_dp !**_**!

            strt = 4

        else q0_chk2

            strt = 1

        end if q0_chk2

        df_loop: do df = strt, Ndof

            omega = W(df)
            eg_vec = Dyn(:, df)

            drv_loop2: do kk = 1, 3

                dDyndq_ev = matmul( dDyn_dq(:,:,kk), eg_vec )
                cmplx_vel(kk)  = dot_product( eg_vec, dDyndq_ev ) / ( 2.0_dp * omega )

            end do drv_loop2

            grp_vel(:, df, ii) = dble( cmplx_vel ) * vg_ms

            nBE(df, ii) = 1.0_dp / ( dexp(BEConv * dabs(omega) / T) - 1.0_dp )

            OnebyFreq(df, ii) = 1.0_dp / ( omega * wTHz ) !**_**!

            !*! debug !*!
            !*! img1_chk2: if ( this_image() == 1) then
            !*!     write(*, *) "Img_vel: ", (dimag( cmplx_vel ) * vg_ms)
            !*!     write(*, *) "Real_vel: ", (dble( cmplx_vel ) * vg_ms)
            !*! end if img1_chk2
            !*! debug !*!

        end do df_loop

        !! ****************** Extra addition for group velocity ****************** !!

        freq(:, ii) = W * wTHz

        Evec(:, :, ii) = Dyn(:, :)

        deallocate( Dyn, dDyn_dq )
        deallocate( W )

        !*! debug !*!
        !*! write(*, *) "q = ", q
        !*! write(*, *)
        !*! debug !*!

        if ( LongEW ) deallocate(DynEw, dDynEw_dq)

    end do q_loop

    deallocate( neg_check, eg_vec, dDyndq_ev )

end subroutine get_phonon


subroutine get_phonon_singleq(sys, FC2, q, EwaldConst, T, LongEw, &
                            & Evec, freq, nBE, grp_vel)

    implicit none

    type(cell), intent(in)                                      :: sys
    type(FC2type), intent(in)                                   :: FC2
    real(dp), dimension(3), intent(in)                          :: q
    type(EwaldParam), intent(in)                                :: EwaldConst
    real(dp), intent(in)                                        :: T
    logical, intent(in)                                         :: LongEW

    complex(dp), dimension(:, :), allocatable, intent(out)      :: Evec
    real(dp), dimension(:), allocatable, intent(out)            :: freq
    real(dp), dimension(:), allocatable, intent(out)            :: nBE
    real(dp), dimension(:, :), allocatable, intent(out)         :: grp_vel

    !! ====================================== Local variables ====================================== !!

    complex(dp), dimension(:,:), allocatable                :: Dyn
    complex(dp), dimension(:,:), allocatable                :: DynEw

    !! ****************** Extra addition for group velocity ****************** !!
    complex(dp), dimension(:,:), allocatable                :: DynEwhp
    complex(dp), dimension(:,:), allocatable                :: DynEwhm

    complex(dp), dimension(:), allocatable                  :: eg_vec
    complex(dp), dimension(:), allocatable                  :: dDyndq_ev
    complex(dp), dimension(3)                               :: cmplx_vel

    complex(dp), dimension(:,:,:), allocatable              :: dDyn_dq
    complex(dp), dimension(:,:,:), allocatable              :: dDynEw_dq
    real(dp), dimension(3)                                  :: q_h, dq
    real(dp)                                                :: h, omega
    integer                                                 :: kk, df, strt
    !! ****************** Extra addition for group velocity ****************** !!

    real(dp), dimension(3)                                  :: qhat, qdir
    real(dp), dimension(:), allocatable                     :: W
    integer                                                 :: basis, jj
    logical, dimension(:), allocatable                      :: neg_check

    !! ====================================== Local variables ====================================== !!

    basis = sys%natm

    allocate(neg_check(3*basis))

    allocate( freq(3*basis) )
    allocate( nBE(3*basis) )
    nBE = 0.0_dp
    allocate( Evec(3*basis, 3*basis) )
    allocate( grp_vel(3, 3*basis) )

    allocate( eg_vec(3*basis) )
    allocate( dDyndq_ev(3*basis) )
    eg_vec = dcmplx(0.0_dp, 0.0_dp)
    dDyndq_ev = dcmplx(0.0_dp, 0.0_dp)

    qdir = (/0.0034_dp, 0.00_dp, 0.001234_dp/)

    call DynMatgv_SR(sys, FC2, q, Dyn, dDyn_dq)

    Ewald_chk: if ( LongEw ) then

        q0_chk1: if ( all(dabs(q) < cmp_prec) ) then

            qhat = q + qdir

            h = norm2( qhat ) / 8.0_dp

            call DynMat_NA(sys, EwaldConst, qhat, DynEw)

            !! ********** Switch to finite difference method for dDynEw_dq at q = 0 ********** !!

            allocate( dDynEw_dq(3*basis, 3*basis, 3) )
            dDynEw_dq = dcmplx(0.0_dp, 0.0_dp)
            
            drv_loop1: do kk = 1, 3

                dq = 0.0_dp
                dq(kk) = 0.5_dp*h

                q_h = q + dq
                call DynMat_Ew(sys, EwaldConst, q_h, DynEwhp)

                q_h = q - dq
                call DynMat_Ew(sys, EwaldConst, q_h, DynEwhm)

                dDynEw_dq(:, :, kk) = (DynEwhp - DynEwhm) / h

                deallocate( DynEwhp, DynEwhm )

            end do drv_loop1

            !! ********** Switch to finite difference method for dDynEw_dq at q = 0 ********** !!
    
        else q0_chk1

            qhat = q
            call DynMatgv_Ew(sys, EwaldConst, qhat, DynEw, dDynEw_dq)
    
        end if q0_chk1
    
        Dyn = Dyn + DynEW
        !! ****************** Extra addition for group velocity ****************** !!
        dDyn_dq = dDyn_dq + dDynEw_dq
        !! ****************** Extra addition for group velocity ****************** !!
        
    end if Ewald_chk

    call SolveEigen(Dyn, W)

    neg_check = ( (dabs(W) > EPS) .and. (W < 0.0_dp) )
    W = dsqrt(dabs(W)) !* wTHz

    do jj = 1, size(neg_check)

        warning: if ( neg_check(jj) ) then 
            W(jj) = W(jj) * (-1.0_dp)

            write(*, 54) q, jj, (W(jj) * wTHz)
            54 FORMAT("WARNING: Imaginary frequency at q = (", &
                     & 2(F9.4, ', '), F9.4, "), branch = ", I4, ". Freqency (omega) = ", F9.4, "rad.THz")
        end if warning

    end do

    !! ****************** Extra addition for group velocity ****************** !!
    cmplx_vel = dcmplx(0.0_dp, 0.0_dp)

    q0_chk2: if ( all(dabs(q) < cmp_prec) ) then

        grp_vel(:, 1) = 0.0_dp
        grp_vel(:, 2) = 0.0_dp
        grp_vel(:, 3) = 0.0_dp

        nBE(1) = 0.0_dp
        nBE(2) = 0.0_dp
        nBE(3) = 0.0_dp

        strt = 4

    else q0_chk2

        strt = 1

    end if q0_chk2

    df_loop: do df = strt, 3*basis 

        omega = W(df)
        eg_vec = Dyn(:, df)

        drv_loop2: do kk = 1, 3

            dDyndq_ev = matmul( dDyn_dq(:,:,kk), eg_vec )
            cmplx_vel(kk)  = dot_product( eg_vec, dDyndq_ev ) / ( 2.0_dp * omega )

        end do drv_loop2

        grp_vel(:, df) = dble( cmplx_vel ) * vg_ms

        nBE(df) = 1.0_dp / ( dexp(BEConv * dabs(omega) / T) - 1.0_dp )

    end do df_loop

    !! ****************** Extra addition for group velocity ****************** !!

    freq(:) = W * wTHz

    Evec(:, :) = Dyn(:, :)

    deallocate( Dyn, dDyn_dq )
    deallocate( W )

    if ( LongEW ) deallocate(DynEw, dDynEw_dq)

    deallocate( neg_check, eg_vec, dDyndq_ev )

end subroutine get_phonon_singleq


end module DynaMat

