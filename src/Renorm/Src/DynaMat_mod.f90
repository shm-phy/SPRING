
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
    use constants,      only : PI, eV_A2, THzConv, EPS, cmp_prec
    use unit_cell,      only : cell
    use FC2_mod,        only : FC2type
    use EwaldMod,       only : EwaldParam
    use EwaldSum,       only : RealSpaceEW, KSpaceEW, CombineZStar
    use mklWrap,        only : SolveEigen
    implicit none
    private

    public  :: get_phonon_disp, SetDynEw_allq, get_freq_ev !DynMat_Ew, DynMat_NA

contains

subroutine DynMat_SR(sys, FC2, q, Dyn)

    implicit none
    complex(dp), parameter      :: iu = dcmplx(0.0_dp, 1.0_dp)

    type(cell), intent(in)                                  :: sys
    type(FC2type), intent(in)                               :: FC2
    real(dp), dimension(3), intent(in)                      :: q
    complex(dp), dimension(:,:), allocatable, intent(out)   :: Dyn

    integer                                 :: N2crt(3), my_nu
    real(dp)                                :: qdotR, N2xyz(3)
    complex(dp)                             :: expnt, dm_el
    integer                                 :: basis, mu, nu, alpha, beta, &
                                               Natm, N2, col_no, row_no

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

subroutine DynMat_Ew(sys, EwaldConst, q, DynEw)

    implicit none

    type(cell), intent(in)                                      :: sys
    type(EwaldParam), intent(in)                                :: EwaldConst
    real(dp), dimension(3), intent(in)                          :: q

    complex(dp), allocatable, intent(out)                       :: DynEw(:, :)


    complex(dp), allocatable, dimension(:,:,:,:)    :: CRSpace, CKSpace &
                                                     &, Cq, CqZ
    integer                                         :: Nbasis, mu, nu, &
                                                     & alpha, beta
    integer                                         :: col_no, row_no

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


subroutine DynMat_NA(sys, EwaldConst, q, DynNA)

    implicit none

    type(cell), intent(in)                                      :: sys
    type(EwaldParam), intent(in)                                :: EwaldConst
    real(dp), dimension(3), intent(in)                          :: q

    complex(dp), allocatable, intent(out)                       :: DynNA(:, :)


    real(dp), dimension(3)                          :: qzero
    complex(dp), allocatable, dimension(:,:,:,:)    :: CRSpace, CKSpace &
                                                     &, Cq0, Cq0Z
    complex(dp), allocatable                        :: DynNA_q0(:, :)

    integer                                         :: Nbasis, mu, nu, &
                                                     & alpha, beta
    integer                                         :: col_no, row_no

    real(dp), dimension(3)                          :: eps_q
    real(dp)                                        :: denom, term1, term2, &
                                                     & mass_fac

    write(*, 75) q
    75 FORMAT('Calculating non-analytic term at q = 0, with qhat = (', 3E11.4 ,')')

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


subroutine SetDynEw_allq(sys, q_points, EwaldConst)

    implicit none

    type(cell), intent(in)                                      :: sys
    real(dp), dimension(:, :), intent(in)                       :: q_points
    type(EwaldParam), intent(inout)                             :: EwaldConst

    real(dp), dimension(3)                                      :: q, qhat
    complex(dp), dimension(:,:), allocatable                    :: DynEw
    complex(dp), dimension(:,:,:), allocatable                  :: DynEwq

    integer                                                     :: nbasis, ii, Nq

    nbasis = sys%natm
    Nq = size(q_points, 2)
    allocate( DynEwq(3*nbasis, 3*nbasis, Nq) )

    loop_q: do ii = 1, Nq

        q = q_points(:, ii)

        q0_chk: if ( all(dabs(q) < cmp_prec) ) then

            endPnt_chk: if ( ii /= Nq ) then
                qhat = q_points(:, ii+1) - q_points(:, ii)
    
            else endPnt_chk
                qhat = q_points(:, ii) - q_points(:, ii-1)

            end if endPnt_chk
    
            call DynMat_NA(sys, EwaldConst, qhat, DynEw)
    
        else q0_chk
            qhat = q
            call DynMat_Ew(sys, EwaldConst, qhat, DynEw)
    
        end if q0_chk
    
        DynEwq(:, :, ii) = DynEw(:, :)
        deallocate(DynEw)

    end do loop_q

    call EwaldConst%set_DynEw(DynEwq)
    deallocate(DynEwq)

end subroutine SetDynEw_allq


subroutine get_phonon_disp(sys, FC2, qpoints, EwaldConst, LongEw, freq_dat)

    implicit none

    type(cell), intent(in)                                      :: sys
    type(FC2type), intent(in)                                   :: FC2
    real(dp), dimension(:, :), intent(in)                       :: qpoints
    type(EwaldParam), intent(in)                                :: EwaldConst
    logical, intent(in)                                         :: LongEW

    real(dp), dimension(:, :), allocatable, intent(out)         :: freq_dat

    real(dp), dimension(3)                                      :: q
    complex(dp), dimension(:,:), allocatable                    :: Dyn
    complex(dp), dimension(:,:), allocatable                    :: DynEw
    real(dp), dimension(:), allocatable                         :: W
    logical, dimension(:), allocatable                          :: neg_check

    integer                                                     :: basis, ii, jj, Nq

    basis = sys%natm
    Nq = size(qpoints, 2)

    allocate(freq_dat((3*basis+1), Nq))
    allocate(neg_check(3*basis))

    q_loop: do ii = 1, Nq

        q = qpoints(2:4, ii)
        call DynMat_SR(sys, FC2, q, Dyn)

        Ewald_chk: if ( LongEw ) then

            DynEw = EwaldConst%DynEwq(:, :, ii)
            Dyn = Dyn + DynEW

        end if Ewald_chk

        call SolveEigen(Dyn, W)

        neg_check = ( (dabs(W) > EPS) .and. (W < 0.0) )
        W = dsqrt(dabs(W)) * THzConv
        do jj = 1, size(neg_check)
            if ( neg_check(jj) ) W(jj) = W(jj) * (-1.0_dp)
        end do

        freq_dat(1, ii) = qpoints(1, ii)
        freq_dat(2:, ii) = W

        deallocate(Dyn)
        deallocate(W)
        if ( LongEW ) deallocate(DynEw)

    end do q_loop

end subroutine get_phonon_disp


subroutine get_freq_ev(sys, FC2, qpoints, EwaldConst, LongEw, freq_dat, Evec)

    implicit none

    type(cell), intent(in)                                      :: sys
    type(FC2type), intent(in)                                   :: FC2
    real(dp), dimension(:, :), intent(in)                       :: qpoints
    type(EwaldParam), intent(in)                                :: EwaldConst
    logical, intent(in)                                         :: LongEW

    real(dp), dimension(:, :), allocatable, intent(out)         :: freq_dat
    complex(dp), dimension(:, :, :), allocatable, intent(out)   :: Evec

    real(dp), dimension(3)                                      :: q
    complex(dp), dimension(:,:), allocatable                    :: Dyn
    complex(dp), dimension(:,:), allocatable                    :: DynEw
    real(dp), dimension(:), allocatable                         :: W
    logical, dimension(:), allocatable                          :: neg_check

    integer                                                     :: basis, ii, jj, Nq


    basis = sys%natm
    Nq = size(qpoints, 2)

    allocate( freq_dat(3*basis, Nq) )
    allocate( neg_check(3*basis) )

    allocate( Evec(3*basis, 3*basis, Nq) ) 

    q_loop: do ii = 1, Nq

        q = qpoints(:, ii)
        call DynMat_SR(sys, FC2, q, Dyn)

        Ewald_chk: if ( LongEw ) then

            DynEw = EwaldConst%DynEwq(:, :, ii)
            Dyn = Dyn + DynEW

        end if Ewald_chk

        call SolveEigen(Dyn, W)

        Evec(:, :, ii) = Dyn(:, :)

        neg_check = ( (dabs(W) > EPS) .and. (W < 0.0_dp) )
        W = dsqrt(dabs(W)) !* THzConv

        do jj = 1, size(neg_check)
            !write(*, *) W(jj)
            if ( neg_check(jj) ) then
                !write(*, *) W(jj)
                W(jj) = W(jj) * (-1.0_dp)
            end if
        end do

        freq_dat(:, ii) = W

        deallocate(Dyn)
        deallocate(W)
        if ( LongEW ) deallocate(DynEw)

    end do q_loop

end subroutine get_freq_ev

end module DynaMat

