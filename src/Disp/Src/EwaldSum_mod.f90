
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

    public  :: RealSpaceEW, KSpaceEW, CombineZStar


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

end module EwaldSum

