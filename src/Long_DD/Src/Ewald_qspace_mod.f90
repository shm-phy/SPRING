
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

module EwaldSum

    use kinds,          only : dp
    use constants,      only : PI, cmp_prec
    use unit_cell,      only : cell
    use CreateMesh,     only : mesh_points
    use IndxChng,       only : fold_indx
    use mklWrap,        only : dinverse, ddet
    implicit none
    private

    public  :: DoEwaldSum, Inverse_FT

contains

    subroutine DoEwaldSum(sys, qmesh, Rmesh, Gmesh, Lmb, decide_EwParam, prec_Ew, qpoints, CEw)
    
        implicit none
    
        type(cell), intent(in)                          :: sys
        integer, dimension(3), intent(in)               :: qmesh, Rmesh, Gmesh
        real(dp), intent(in)                            :: Lmb, prec_Ew
        logical, intent(in)                             :: decide_EwParam
        real(dp), allocatable, intent(out)              :: qpoints(:, :)
        complex(dp), allocatable, intent(out)           :: CEw(:, :, :, :, :)
    
        real(dp), dimension(3, 3)                       :: eps_inv
        real(dp)                                        :: det
        real(dp), dimension(3)                          :: q
        real(dp), allocatable, dimension(:, :)          :: Rpoints, Gpoints
        complex(dp), allocatable, dimension(:,:,:,:)    :: CConst, CRSpace, CKSpace &
                                                         &, Cq, CqZ
        complex(dp), allocatable, dimension(:,:,:)      :: CSumq0
    
        integer                                         :: Nbasis, mu, alpha, &
                                                         & beta, Nq, qp
        real(dp)                                        :: Ew_alpha
        integer                                         :: NKcut
        integer, dimension(3)                           :: Ew_Gmesh
    
        eps_inv = sys%eps
    
        eps_inv = transpose(eps_inv)
        call dinverse(eps_inv)
        eps_inv = transpose(eps_inv)
    
        chk_EwPrm: if ( decide_EwParam ) then
    
            call DecideEwaldParameters(sys%eps, eps_inv, sys%latvec, sys%G, sys%sup_cell, &
                                     & prec_Ew, Ew_alpha, NKcut)
            Ew_Gmesh(:) = 2*NKcut + 1
    
        else chk_EwPrm
            Ew_alpha = Lmb
            Ew_Gmesh = Gmesh
    
        end if chk_EwPrm
    
        write(*, 150) Ew_alpha
        150 FORMAT('The parameter alpha for Ewald Sum: ', F7.4)
        write(*, 250) Ew_Gmesh
        250 FORMAT('The NKcut (K-space cutoff) for Ewald Sum: ', 3I3)
    
        call ddet(sys%eps, det)
    
        Nbasis = sys%natm
        allocate(Cq(3, 3, Nbasis, Nbasis))
        Cq = dcmplx(0.0_dp, 0.0_dp)
    
        call mesh_points(cell_vec = sys%G, mesh = qmesh, div_cell = .true., mesh_cord = qpoints)
        call mesh_points(sys%latvec, Rmesh, .false., Rpoints)
        call mesh_points(sys%G, Ew_Gmesh, .false., Gpoints)
    
        call ConstTerm(eps_inv, det, Ew_alpha, Nbasis, CConst)
    
        call SumRule(sys, Rpoints, Gpoints, eps_inv, det, Ew_alpha, CConst, CSumq0)
    
        Nq = size(qpoints, 2)
        allocate(CEw(3, 3, Nbasis, Nbasis, Nq))
    
        !******* Loop over all q *********!
        q_loop: do qp = 1, Nq
    
            q = qpoints(:, qp)
            
            call RealSpaceEW(q, Rpoints, sys%tau, eps_inv, det, Ew_alpha, CRSpace)
    
            call KSpaceEW(q, Gpoints, sys%tau, sys%eps, sys%vol, Ew_alpha, CKSpace)
            
            Cq = CKSpace - CRSpace - CConst
            call CombineZStar(Cq, sys%Zstar, CqZ)
            
            mu_loop: do mu = 1, Nbasis
                a_loop: do alpha = 1, 3
                    b_loop: do beta = 1, 3
                    
                          CqZ(beta, alpha, mu, mu) = &
                        & CqZ(beta, alpha, mu, mu) - &
                        & CSumq0(beta, alpha, mu)
                    
                    end do b_loop
                end do a_loop
            end do mu_loop
    
            CEw(:, :, :, :, qp) = CqZ
    
            deallocate(CRSpace)
            deallocate(CKSpace)
            Cq = dcmplx(0.0_dp, 0.0_dp)
            deallocate(CqZ)
    
        end do q_loop
        !******* Loop over all q *********!
    
        deallocate(CConst)
        deallocate(CSumq0)
        deallocate(Gpoints)
        deallocate(Rpoints)
    
    end subroutine DoEwaldSum
    
    
    subroutine DecideEwaldParameters(eps, eps_inv, lat_vec, G, SupDim, prec_Ew, Ew_alpha, NKcut)
    
        implicit none
    
        real(dp), dimension(3, 3), intent(in)                   :: eps, eps_inv, &
                                                                 & lat_vec, G
        integer, dimension(3), intent(in)                       :: SupDim
        real(dp), intent(in)                                    :: prec_Ew
    
        real(dp), intent(out)                                   :: Ew_alpha
        integer, intent(out)                                    :: NKcut
    
        integer, dimension(3)           :: Ncut
        real(dp), dimension(3)          :: Rcut, dV, delV
        real(dp)                        :: p, D, Kcut, normG
    
        p = (-1.0_dp * dlog(prec_Ew))
    
        Ncut = (SupDim / 2) + 1
        Ncut(2:) = 0 !****!
        Rcut = matmul(lat_vec, dble(Ncut)) 
    
        dV = Rcut 
        delV = matmul(dV, eps_inv)
        D = dot_product(delV, dV)
        D = dsqrt(D)
        Ew_alpha = dsqrt(p) / D
    
        normG = norm2(G(:,1))
        Kcut = 2.0_dp * p / (D * dsqrt( (eps(1,1) + eps(2,2) + eps(3,3))/3.0_dp ))
        NKcut = ceiling(KCut/normG)
    
    end subroutine DecideEwaldParameters
    
    
    subroutine ConstTerm(eps_inv, det_eps, Lmb, Nbasis, CConst)
    
        implicit none
    
        real(dp), dimension(3, 3), intent(in)               :: eps_inv
        real(dp), intent(in)                                :: det_eps, Lmb
        integer, intent(in)                                 :: Nbasis
        complex(dp), allocatable, intent(out)               :: CConst(:, :, :, :)
    
        real(dp)                                            :: cnst
        integer                                             :: mu, nu, alpha, beta
    
        allocate(CConst(3, 3, Nbasis, Nbasis))
        CConst = dcmplx(0.0_dp, 0.0_dp)
    
        do mu = 1, Nbasis
            do nu = 1, Nbasis
    
                if (mu == nu) then
                    do alpha = 1, 3
                        do beta = 1, 3
    
                            CConst(beta, alpha, nu, mu) = eps_inv(beta, alpha)
    
                        end do
                    end do
                end if
    
            end do
        end do
    
        cnst = 4.0_dp * (Lmb ** 3) / (3.0_dp * dsqrt(PI) * dsqrt(det_eps))
    
        CConst = CConst * cnst
    
    end subroutine ConstTerm
    
    
    subroutine  RealSpaceEW(q, Ra, tau, eps_inv, det_eps, Lmb, CRSpace)
    
        implicit none
        complex(dp), parameter                              :: iu = dcmplx(0.0_dp, 1.0_dp)
    
        real(dp), dimension(3), intent(in)                  :: q
        real(dp), contiguous, dimension(:, :), intent(in)   :: Ra
        real(dp), contiguous, dimension(:, :), intent(in)   :: tau
        real(dp), contiguous, dimension(:, :), intent(in)   :: eps_inv
        real(dp), intent(in)                                :: det_eps, Lmb
        complex(dp), allocatable, intent(out)               :: CRSpace(:, :, :, :)
    
        real(dp), dimension(3)                              :: R, dV, delV
        real(dp)                                            :: q_dot_R, DC, hmn1, &
                                                             & hmn2, hab
        complex(dp)                                         :: expn
    
        integer                                             :: Nbasis, mu, nu, &
                                                             & alpha, beta, NR, rp
    
        Nbasis = size(tau, 2)
        NR = size(Ra, 2)
    
        allocate(CRSpace(3, 3, Nbasis, Nbasis))
        CRSpace = dcmplx(0.0_dp, 0.0_dp)
    
        do rp = 1, NR
    
            R = Ra(:, rp)
            q_dot_R = dot_product(q, R)
            expn = cdexp(iu * q_dot_R)
    
            do mu = 1, Nbasis
                do nu = 1, Nbasis
    
                    if ( all(dabs(R) < cmp_prec) .and. (mu == nu) ) then 
    
                        delV = (/0.0_dp, 0.0_dp, 0.0_dp/)
                        hmn1 = 0.0_dp
                        hmn2 = 0.0_dp
    
                    else
                        dV = R + tau(:, nu) - tau(:, mu)
                        delV = matmul(dV, eps_inv)
                        DC = dot_product(delV, dV)
                        DC = dsqrt(DC)
    
                        delV = Lmb * delV
    
                        DC = Lmb * DC
    
                        call H_func(DC, hmn1, hmn2)
    
                    end if
    
                    do alpha = 1, 3
                        do beta = 1, 3
    
                            hab = (delV(alpha) * delV(beta) * hmn1) - &
                                & (eps_inv(beta, alpha) * hmn2)
    
                            CRSpace(beta, alpha, nu, mu) = &
                          & CRSpace(beta, alpha, nu, mu) + expn * hab 
    
                        end do
                    end do
    
                end do
            end do
        end do
    
        CRSpace = CRSpace * (Lmb ** 3) / dsqrt(det_eps)
    
    end subroutine RealSpaceEW
    
    
    subroutine H_func(ld, h1, h2)
    
        implicit none
    
        real(dp), intent(in)            :: ld
        real(dp), intent(out)           :: h1, h2
    
        real(dp)                        :: cerfn, y2, invy2, invy3, pi_expnt
        
        cerfn = derfc(ld)
    
        y2 = ld ** 2
        invy2 = 1.0_dp / y2
        invy3 = invy2 / ld
    
        pi_expnt = (2.0_dp / dsqrt(PI)) * dexp(-1.0_dp * y2)
    
        h1 = (3.0_dp * invy3  * cerfn + pi_expnt * (3.0_dp * invy2 + 2.0_dp)) * invy2
    
        h2 = cerfn * invy3  +  pi_expnt * invy2
    
    end subroutine H_func
    
    
    subroutine KSpaceEW(q, Ga, tau, eps, vol, Lmb, Cq)
    
        implicit none
        complex(dp), parameter      :: iu = dcmplx(0.0_dp, 1.0_dp)
    
        real(dp), dimension(3), intent(in)                  :: q
        real(dp), contiguous, dimension(:, :), intent(in)   :: Ga
        real(dp), contiguous, dimension(:, :), intent(in)   :: tau
        real(dp), dimension(3, 3), intent(in)               :: eps
        real(dp), intent(in)                                :: vol, Lmb
        complex(dp), allocatable, intent(out)               :: Cq(:, :, :, :)
    
        real(dp), dimension(3)                              :: G, K, eK, tmn
    
        real(dp)                                            :: KeK, KeKL, &
                                                             & KG, K_dot_tau
    
        complex(dp)                                         :: expnt, Kmn, Kab
    
        integer                                             :: Nbasis, mu, nu, &
                                                             & alpha, beta, NG, gp
    
        Nbasis = size(tau, 2)
        NG = size(Ga, 2)
    
        allocate(Cq(3, 3, Nbasis, Nbasis))
        Cq = dcmplx(0.0_dp, 0.0_dp)
    
        do gp = 1, NG
    
            G = Ga(:, gp)
    
            ne0chk: if ( (any(dabs(q) > cmp_prec)) .or. &
                       & (any(dabs(G) > cmp_prec)) ) then
    
                K = G + q
    
                eK = matmul(K, eps)
                KeK = dot_product(K, eK)
    
                KeKL = -1.0_dp * KeK / (4.0_dp * (Lmb ** 2))
                KG = dexp(KeKL) / KeK
    
                do mu = 1, Nbasis
                    do nu = 1, Nbasis
    
                        tmn = tau(:, mu) - tau(:, nu)
                        !tmn = tau(:, nu) - tau(:, mu)
                        K_dot_tau = dot_product(K, tmn)
                        expnt = cdexp(iu * K_dot_tau)
    
                        Kmn = KG * expnt
    
                        do alpha = 1, 3
                            do beta = 1, 3
    
                                Kab = K(alpha) * K(beta) * Kmn
    
                                Cq(beta, alpha, nu, mu) =   &
                                & Cq(beta, alpha, nu, mu) + &
                                & Kab
    
                            end do
                        end do
    
                    end do
                end do
    
            end if ne0chk
        end do
    
        Cq = Cq * (4.0_dp * PI / vol)
    
    end subroutine KSpaceEW
    
    
    subroutine CombineZStar(Cq, Zstar, CqZ)
    
        implicit none
    
        complex(dp), contiguous, dimension(:,:,:,:), intent(in)     :: Cq
        real(dp), contiguous, dimension(:,:,:), intent(in)          :: Zstar
        complex(dp), allocatable, dimension(:,:,:,:), intent(out)   :: CqZ
    
        integer                                                     :: Nbasis, mu, &
                                                                     & nu, alpha,  &
                                                                     & beta, ap, bp
        Nbasis = size(Cq, 4)
        allocate(CqZ(3, 3, Nbasis, Nbasis))
        CqZ = dcmplx(0.0_dp, 0.0_dp)
    
        mu_loop: do mu = 1, Nbasis
            nu_loop: do nu = 1, Nbasis
                a_loop: do alpha = 1, 3
                    b_loop: do beta = 1, 3
    
                        ap_loop: do ap = 1, 3
                            bp_loop: do bp = 1, 3
    
                                  CqZ(beta, alpha, nu, mu) = &
                                & CqZ(beta, alpha, nu, mu) + &
                                & ( Zstar(ap, alpha, mu) *   &
                                &   Zstar(bp, beta, nu) *    &
                                &   Cq(bp, ap, nu, mu) )
    
                            end do bp_loop
                        end do ap_loop
    
                    end do b_loop
                end do a_loop
            end do nu_loop
        end do mu_loop
    
    end subroutine CombineZstar
    
    
    subroutine SumRule(sys, Rpoints, Gpoints, eps_inv, det, Lmb, CConst, Csumq0)
    
        implicit none
    
        type(cell), intent(in)                                  :: sys
        real(dp), contiguous, dimension(:, :), intent(in)       :: Rpoints, Gpoints
        real(dp), dimension(3, 3), intent(in)                   :: eps_inv
        real(dp), intent(in)                                    :: det, Lmb
        complex(dp), contiguous, dimension(:,:,:,:), intent(in) :: CConst
        complex(dp), allocatable, intent(out)                   :: CSumq0(:, :, :)
    
        complex(dp), allocatable, dimension(:,:,:,:)            :: CRq0, CKq0, Cq0 &
                                                                &, Cq0Z
    
        real(dp), dimension(3)                                  :: q
        integer                                                 :: Nbasis, mu, nu, &
                                                                 & alpha, beta
    
        Nbasis = sys%natm
        allocate(Cq0(3, 3, Nbasis, Nbasis))
        Cq0 = dcmplx(0.0_dp, 0.0_dp)
    
        q = (/0.0_dp, 0.0_dp, 0.0_dp/)
        call RealSpaceEW(q, Rpoints, sys%tau, eps_inv, det, Lmb, CRq0)
        call KSpaceEW(q, Gpoints, sys%tau, sys%eps, sys%vol, Lmb, CKq0)
        Cq0 = CKq0 - CRq0 - CConst
        call CombineZstar(Cq0, sys%Zstar, Cq0Z)
    
        allocate(CSumq0(3, 3, Nbasis))
        CSumq0 = dcmplx(0.0_dp, 0.0_dp)
    
        mu_loop: do mu = 1, Nbasis
            a_loop: do alpha = 1, 3
                b_loop: do beta = 1, 3
    
                    nu_loop: do nu = 1, Nbasis
                          CSumq0(beta, alpha, mu) = &
                        & CSumq0(beta, alpha, mu) + &
                        & Cq0Z(beta, alpha, nu, mu)
                    end do nu_loop
    
                end do b_loop
            end do a_loop
        end do mu_loop
    
        deallocate(CRq0)
        deallocate(CKq0)
        deallocate(Cq0)
        deallocate(Cq0Z)
    
    end subroutine SumRule
    
    
    subroutine GetRSpacePos(sys, RSup, RCb)
    
        implicit none
    
        type(cell), intent(in)                                  :: sys
        integer, dimension(3), intent(in)                       :: RSup
        integer, dimension(:,:,:), allocatable, intent(out)     :: RCb
    
        integer                                                 :: num_col, Nbasis,&  
                                                                 & mu, ix, iy, iz, &
                                                                 & nu, num
        integer, dimension(3)                                   :: indx_C, fold_C
    
        Nbasis = sys%natm
        num_col = product(RSup)
    
        allocate(RCb(4, num_col*Nbasis, Nbasis))
    
        mu_loop: do mu = 1, Nbasis
            num = 1
            x_loop: do ix = 0, (Rsup(1)-1)
                y_loop: do iy = 0, (Rsup(2)-1)
                    z_loop: do iz = 0, (Rsup(3)-1)
                                
                        indx_C = (/ix, iy, iz/)
                        call fold_indx(Rsup, indx_C, fold_C)
    
                        nu_loop: do nu = 1, Nbasis
                                    
                            RCb(1:3, num, mu) = fold_C
                            RCb(4, num, mu) = nu
                            num = num + 1
                                    
                        end do nu_loop
    
                    end do z_loop
                end do y_loop
            end do x_loop
        end do mu_loop
    
    end subroutine GetRSpacePos
    
    
    subroutine get_index_num(Rcb, mu, indx, N)
    
        implicit none
    
        integer, dimension(:,:,:), intent(in)       :: RCb
        integer, intent(in)                         :: mu
        integer, dimension(4), intent(in)           :: indx
        integer, intent(out)                        :: N
    
        integer, dimension(4)                       :: tmp_arr
        integer                                     :: nn, num_col
    
        num_col = size(RCb, 2)
    
        N = 0
        loop_atm: do nn = 1, num_col
    
            tmp_arr = Rcb(:, nn, mu)
    
            chk: if ( all((tmp_arr-indx) .eq. 0) ) then
                N = nn
                exit loop_atm
            end if chk
    
        end do loop_atm
    
        if ( N == 0 ) then
            write(*, *) "Error: Index not found"
        end if
    
    end subroutine get_index_num
    

    subroutine Inverse_FT(sys, CEw, qpoints, IFC, RCb)
    
        implicit none
        real(dp), parameter         :: eV_A2 = ((1.602176634_dp * 2.99792458_dp)**2 &
                                             & * 6.241509074460763_dp * 0.1_dp)
        complex(dp), parameter      :: iu = dcmplx(0.0_dp, 1.0_dp)
    
        type(cell), intent(in)                                      :: sys
        complex(dp), contiguous, dimension(:,:,:,:,:), intent(in)   :: CEw
        real(dp), contiguous, dimension(:,:), intent(in)            :: qpoints 
        complex(dp), allocatable, dimension(:,:,:,:), intent(out)   :: IFC
        integer, allocatable, dimension(:,:,:), intent(out)         :: RCb
    
        integer, dimension(:,:), allocatable                        :: Rc
        integer, dimension(3)                                       :: Nb
        real(dp), dimension(3)                                      :: q, Rb
        real(dp)                                                    :: qdotR
        complex(dp)                                                 :: expn
        integer, dimension(4)                                       :: indx
        integer                                                     :: qp, Rp, mu, nu&
                                                                    &,N2, alpha, beta&
                                                                    &,Nbasis, NC, Nq
    
        Nbasis = sys%natm
        NC = product(sys%sup_cell)
    
        call GetRSpacePos(sys, sys%sup_cell, RCb)
    
        allocate(Rc(3, NC))
        Rc = RCb(1:3, ::Nbasis, 1)
    
        Nq = size(qpoints, 2)
    
        allocate(IFC(3, 3, NC*Nbasis, Nbasis))
        IFC = dcmplx(0.0_dp, 0.0_dp)
    
        q_sum: do qp = 1, Nq
            q = qpoints(:, qp)
    
            !********************** Sum over all q ***************************!
            R_loop: do Rp = 1, NC
    
                Nb = Rc(:, Rp)
                indx(1:3) = Nb
    
                Rb = matmul(sys%latvec, Nb)
                qdotR = dot_product(q, Rb)
                expn = cdexp(-1.0_dp * iu * qdotR)
                !expn = dcmplx( dcos(qdotR), (-1.0_dp * dsin(qdotR)) )
    
                mu_loop: do mu = 1, Nbasis
                    nu_loop: do nu = 1, Nbasis
    
                        indx(4) = nu
                        call get_index_num(RCb, mu, indx, N2)
    
                        a_loop: do alpha = 1, 3
                            b_loop: do beta = 1, 3
    
                                  IFC(beta, alpha, N2, mu) = &
                                & IFC(beta, alpha, N2, mu) + &
                                & CEw(beta, alpha, nu, mu, qp) * expn 
    
                            end do b_loop
                        end do a_loop
    
                    end do nu_loop
                end do mu_loop
    
            end do R_loop
            !********************** Sum over all q ***************************!
    
        end do q_sum
    
        IFC = (IFC * eV_A2) / Nq
    
        deallocate(Rc)
    
    end subroutine Inverse_FT

end module EwaldSum

