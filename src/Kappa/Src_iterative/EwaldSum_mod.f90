
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

    implicit none
    private

    public  :: RealSpaceEW, RealSpaceEW_gv, &
             & KSpaceEW, KSpaceEW_gv, &
             & CombineZStar, CombineZStar_gv


contains

subroutine  RealSpaceEW(q, Ra, tau, eps_inv, det_eps, Lmb, CRSpace)

    implicit none
    complex(dp), parameter                              :: iu = dcmplx(0.0_dp, 1.0_dp)

    real(dp), dimension(3), intent(in)                  :: q
    real(dp), contiguous, dimension(:, :), intent(in)   :: Ra
    real(dp), contiguous, dimension(:, :), intent(in)   :: tau
    real(dp), contiguous, dimension(:, :), intent(in)   :: eps_inv
    real(dp), intent(in)                                :: det_eps, Lmb
    complex(dp), allocatable, intent(out)               :: CRSpace(:, :, :, :)

    ! ================================== Local Variables ================================== !
    complex(dp)                                         :: expn

    real(dp), dimension(3)                              :: R, dV, delV
    real(dp)                                            :: q_dot_R, DC, hmn1, &
                                                         & hmn2, hab

    integer                                             :: Nbasis, mu, nu, &
                                                         & alpha, beta, NR, rp
    ! ================================== Local Variables ================================== !

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


subroutine  RealSpaceEW_gv(q, Ra, tau, eps_inv, det_eps, Lmb, CRSpace, &
                         & dCR_dq)

    implicit none
    complex(dp), parameter                              :: iu = dcmplx(0.0_dp, 1.0_dp)

    real(dp), dimension(3), intent(in)                  :: q
    real(dp), contiguous, dimension(:, :), intent(in)   :: Ra
    real(dp), contiguous, dimension(:, :), intent(in)   :: tau
    real(dp), contiguous, dimension(:, :), intent(in)   :: eps_inv
    real(dp), intent(in)                                :: det_eps, Lmb
    complex(dp), allocatable, intent(out)               :: CRSpace(:, :, :, :)
    !! ****************** Extra addition for group velocity ****************** !!
    complex(dp), allocatable, intent(out)               :: dCR_dq(:, :, :, :, :)
    !! ****************** Extra addition for group velocity ****************** !!

    ! ================================ Local variables ================================ !
    complex(dp)                                         :: expn

    real(dp), dimension(3)                              :: R, dV, delV
    real(dp)                                            :: q_dot_R, DC, hmn1, &
                                                         & hmn2, hab

    integer                                             :: Nbasis, mu, nu, &
                                                         & alpha, beta, NR, rp, ii
    ! ================================ Local variables ================================ !

    Nbasis = size(tau, 2)
    NR = size(Ra, 2)

    allocate(CRSpace(3, 3, Nbasis, Nbasis))
    CRSpace = dcmplx(0.0_dp, 0.0_dp)

    !! ****************** Extra addition for group velocity ****************** !!
    allocate(dCR_dq(3, 3, Nbasis, Nbasis, 3))
    dCR_dq = dcmplx(0.0_dp, 0.0_dp)
    !! ****************** Extra addition for group velocity ****************** !!

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
                        
                        !! ****************** Extra addition for group velocity ****************** !!
                        drv_loop: do ii = 1, 3

                            dCR_dq(beta, alpha, nu, mu, ii) = &
                          & dCR_dq(beta, alpha, nu, mu, ii) + &
                          & iu * R(ii) * expn * hab

                        end do drv_loop
                        !! ****************** Extra addition for group velocity ****************** !!

                    end do
                end do

            end do
        end do
    end do

    CRSpace = CRSpace * (Lmb ** 3) / dsqrt(det_eps)

    !! ****************** Extra addition for group velocity ****************** !!
    dCR_dq = dCR_dq * (Lmb ** 3) / dsqrt(det_eps)
    !! ****************** Extra addition for group velocity ****************** !!

end subroutine RealSpaceEW_gv


subroutine H_func(ld, h1, h2)

    implicit none

    real(dp), intent(in)            :: ld
    real(dp), intent(out)           :: h1, h2

    ! ================================== Local Variables ================================== !
    real(dp)                        :: cerfn, y2, invy2, invy3, pi_expnt
    ! ================================== Local Variables ================================== !
    
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

    ! ================================== Local Variables ================================== !
    complex(dp)                                         :: expnt, Kmn, Kab

    real(dp), dimension(3)                              :: G, K, eK, tmn

    real(dp)                                            :: KeK, KeKL, &
                                                         & KG, K_dot_tau

    integer                                             :: Nbasis, mu, nu, &
                                                         & alpha, beta, NG, gp
    ! ================================== Local Variables ================================== !

    Nbasis = size(tau, 2)
    NG = size(Ga, 2)

    allocate(Cq(3, 3, Nbasis, Nbasis))
    Cq = dcmplx(0.0_dp, 0.0_dp)

    do gp = 1, NG

        G = Ga(:, gp)

        ne0chk: if ( (any(dabs(q) > cmp_prec)) .or. (any(dabs(G) > cmp_prec)) ) then

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


subroutine KSpaceEW_gv(q, Ga, tau, eps, vol, Lmb, Cq, &
                       dCq_dq)

    implicit none
    complex(dp), parameter                              :: iu = dcmplx(0.0_dp, 1.0_dp)

    real(dp), dimension(3), intent(in)                  :: q
    real(dp), contiguous, dimension(:, :), intent(in)   :: Ga
    real(dp), contiguous, dimension(:, :), intent(in)   :: tau
    real(dp), dimension(3, 3), intent(in)               :: eps
    real(dp), intent(in)                                :: vol, Lmb
    complex(dp), allocatable, intent(out)               :: Cq(:, :, :, :)

    !! ****************** Extra addition for group velocity ****************** !!
    complex(dp), allocatable, intent(out)               :: dCq_dq(:, :, :, :, :)
    !! ****************** Extra addition for group velocity ****************** !!

    ! ================================ Local variables ================================ !
    complex(dp)                                         :: expnt, Kmn, Kab

    real(dp), dimension(3)                              :: G, K, eK, tmn

    !! ****************** Extra addition for group velocity ****************** !!
    real(dp), dimension(3, 3, 3)                        :: K_Mat
    real(dp), dimension(3)                              :: epsK
    real(dp)                                            :: KpLmb
    integer                                             :: ii
    !! ****************** Extra addition for group velocity ****************** !!

    real(dp)                                            :: KeK, KeKL, &
                                                         & KG, K_dot_tau

    integer                                             :: Nbasis, mu, nu, &
                                                         & alpha, beta, NG, gp
    ! ================================ Local variables ================================ !

    Nbasis = size(tau, 2)
    NG = size(Ga, 2)

    allocate(Cq(3, 3, Nbasis, Nbasis))
    Cq = dcmplx(0.0_dp, 0.0_dp)

    !! ****************** Extra addition for group velocity ****************** !!
    allocate( dCq_dq(3, 3, Nbasis, Nbasis, 3) )
    dCq_dq = dcmplx(0.0_dp, 0.0_dp)

    K_Mat = 0.0_dp
    epsK = 0.0_dp
    !! ****************** Extra addition for group velocity ****************** !!

    do gp = 1, NG

        G = Ga(:, gp)

        ne0chk: if ( (any(dabs(q) > cmp_prec)) .or. (any(dabs(G) > cmp_prec)) ) then

            K = G + q

            !! ****************** Extra addition for group velocity ****************** !!

            K_Mat(:, 1, 1) = (/2.0_dp*K(1), K(2),   K(3)/)
            K_Mat(:, 2, 1) = (/   K(2),     0.0_dp, 0.0_dp/)
            K_Mat(:, 3, 1) = (/   K(3),     0.0_dp, 0.0_dp/)

            K_Mat(:, 1, 2) = (/0.0_dp,   K(1),    0.0_dp/)
            K_Mat(:, 2, 2) = (/K(1), 2.0_dp*K(2), K(3)/)
            K_Mat(:, 3, 2) = (/0.0_dp,   K(3),    0.0_dp/)

            K_Mat(:, 1, 3) = (/0.0_dp, 0.0_dp,   K(1)/)
            K_Mat(:, 2, 3) = (/0.0_dp, 0.0_dp,   K(2)/)
            K_Mat(:, 3, 3) = (/K(1),   K(2),   2.0_dp*K(3)/)

            epsK(1) = dot_product( eps(:,1), K ) + dot_product( K, eps(1, :) )
            epsK(2) = dot_product( eps(:,2), K ) + dot_product( K, eps(2, :) )
            epsK(3) = dot_product( eps(:,3), K ) + dot_product( K, eps(3, :) )

            !! ****************** Extra addition for group velocity ****************** !!

            eK = matmul(K, eps)
            KeK = dot_product(K, eK)

            !! ****************** Extra addition for group velocity ****************** !!
            KpLmb = ( (1.0_dp / KeK ) + (0.25_dp / (Lmb ** 2)) )
            !! ****************** Extra addition for group velocity ****************** !!

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

                            !! ****************** Extra addition for group velocity ****************** !!
                            drv_loop: do ii = 1, 3

                                dCq_dq(beta, alpha, nu, mu, ii) = &
                             &  dCq_dq(beta, alpha, nu, mu, ii) + &
                             &  Kab * ( iu * tmn(ii) - epsK(ii) * KpLmb ) + &
                             &  Kmn * K_Mat(beta, alpha, ii)

                            end do drv_loop
                            !! ****************** Extra addition for group velocity ****************** !!

                        end do
                    end do

                end do
            end do

        end if ne0chk
    end do

    Cq = Cq * (4.0_dp * PI / vol)

    !! ****************** Extra addition for group velocity ****************** !!
    dCq_dq = dCq_dq * (4.0_dp * PI / vol)
    !! ****************** Extra addition for group velocity ****************** !!

end subroutine KSpaceEW_gv


subroutine CombineZStar(Cq, Zstar, CqZ)

    implicit none

    complex(dp), contiguous, dimension(:,:,:,:), intent(in)     :: Cq
    real(dp), contiguous, dimension(:,:,:), intent(in)          :: Zstar
    complex(dp), allocatable, dimension(:,:,:,:), intent(out)   :: CqZ

    ! ================================== Local Variables ================================== !
    integer                                                     :: Nbasis, mu, &
                                                                 & nu, alpha,  &
                                                                 & beta, ap, bp
    ! ================================== Local Variables ================================== !

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


subroutine CombineZStar_gv(Cq, dCq_dq, Zstar, CqZ, dCqZ_dq)

    implicit none

    complex(dp), contiguous, dimension(:,:,:,:), intent(in)     :: Cq
    !! ****************** Extra addition for group velocity ****************** !!
    complex(dp), contiguous, dimension(:,:,:,:,:), intent(in)   :: dCq_dq
    !! ****************** Extra addition for group velocity ****************** !!
    real(dp), contiguous, dimension(:,:,:), intent(in)          :: Zstar
    complex(dp), allocatable, dimension(:,:,:,:), intent(out)   :: CqZ
    !! ****************** Extra addition for group velocity ****************** !!
    complex(dp), allocatable, dimension(:,:,:,:,:), intent(out)   :: dCqZ_dq
    !! ****************** Extra addition for group velocity ****************** !!

    ! ================================ Local variables ================================ !

    integer                                                     :: Nbasis, mu, &
                                                                 & nu, alpha,  &
                                                                 & beta, ap, bp, ii

    ! ================================ Local variables ================================ !

    Nbasis = size(Cq, 4)
    allocate(CqZ(3, 3, Nbasis, Nbasis))
    CqZ = dcmplx(0.0_dp, 0.0_dp)

    !! ****************** Extra addition for group velocity ****************** !!
    allocate(dCqZ_dq(3, 3, Nbasis, Nbasis, 3))
    dCqZ_dq = dcmplx(0.0_dp, 0.0_dp)
    !! ****************** Extra addition for group velocity ****************** !!

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
                            
                            !! ****************** Extra addition for group velocity ****************** !!
                            drv_loop: do ii = 1, 3

                                dCqZ_dq(beta, alpha, nu, mu, ii) = &
                              & dCqZ_dq(beta, alpha, nu, mu, ii) + &
                              & ( Zstar(ap, alpha, mu) * Zstar(bp, beta, nu) * &
                              &  dCq_dq(bp, ap, nu, mu, ii) )

                            end do drv_loop
                            !! ****************** Extra addition for group velocity ****************** !!

                        end do bp_loop
                    end do ap_loop
                
                end do b_loop
            end do a_loop
        end do nu_loop
    end do mu_loop

end subroutine CombineZstar_gv

end module EwaldSum

