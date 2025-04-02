
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


module Renorm

    use kinds,          only : dp
    use constants,      only : zero_prec, cmp_prec, MAX_ITER, &
                             & BEConv, A2Conv, THzConv
    use unit_cell,      only : cell
    use FC2_mod,        only : FC2type
    use FC4_mod,        only : FC4type
    use EwaldMod,       only : EwaldParam
    use DynaMat,        only : get_freq_ev

    implicit none
    private

    public  :: Renormalize

contains

subroutine Renormalize(sys, T, Renorm_accu, EwaldConst, &
                     & LongEw, FC2_old, FC4, q_points, FC2_new)

    implicit none
    !*! complex(dp), parameter                          :: iu = dcmplx(0.0_dp, 1.0_dp)

    type(cell), intent(in)                          :: sys
    real(dp), intent(in)                            :: T
    real(dp), intent(in)                            :: Renorm_accu
    type(EwaldParam), intent(in)                    :: EwaldConst
    logical, intent(in)                             :: LongEw
    type(FC2type), intent(in)                       :: FC2_old
    type(FC4type), intent(in)                       :: FC4
    real(dp), dimension(:, :), intent(in)           :: q_points

    type(FC2type), intent(inout)                    :: FC2_new

    ! ======================================== Local Variables ======================================== !
    complex(dp), dimension(:, :, :), allocatable    :: Evec
    real(dp), dimension(:, :), allocatable          :: freq_old, freq_new
    complex(dp), dimension(:, :), allocatable       :: Evec_q
    real(dp), dimension(:), allocatable             :: W_q

    integer, dimension(5)                           :: fC2_shp
    complex(dp), dimension(:,:,:,:,:), allocatable  :: dFC2

    integer, dimension(3)                           :: N3indx, N4indx

    real(dp), dimension(3)                          :: q, R3_R4

    real(dp)                                        :: q_dotR, mass_fac, Omega, &
                                                     & nqs, BE_fac, freq_acc

    complex(dp)                                     :: expnt, eV_cmp1, eV_cmp2, mul_fac

    integer                                         :: eta, ro, &
                                                     & row_no1, row_no2
    integer                                         :: Nq, Nbasis, dof, fC2_N2, Natm4

    integer                                         :: nn, ss, mu, N2nu, &
                                                     & N3eta, N4ro, aa,  &
                                                     & bb, cc, dd, iter_no, strt, omp_chunk

    logical                                         :: IsZero
    ! ======================================== Local Variables ======================================== !


    FC2_new = FC2_old
    call get_freq_ev(sys, FC2_old, q_points, EwaldConst, LongEw, freq_new, Evec)

    Nq = size(q_points, 2)
    Nbasis = sys%natm
    dof = 3*Nbasis

    allocate( Evec_q(dof, dof) )
    allocate( W_q(dof) )

    fC2_shp = shape(FC2_old%fC)
    allocate( dFC2(fC2_shp(1), fC2_shp(2), &
                 & fC2_shp(3), fC2_shp(4), fC2_shp(5)) )

    write(*, *)
    write(*, 129)
    129 FORMAT('Starting self consistence iteration loop')
    write(*, *)

    omp_chunk = 12
    iter_no = 1
    scf_iter_loop: do 

        write(*, 175) iter_no
        write(*, 595)

        dFC2 = dcmplx(0.0_dp, 0.0_dp)
        Evec_q = dcmplx(0.0_dp, 0.0_dp)
        W_q = 0.0_dp

        !$omp parallel do default(shared) &
                    !$omp schedule(dynamic, omp_chunk) &
                    !$omp reduction(+: dFC2) &
                    !$omp shared(FC4) &
                    !$omp shared(T, q_points) & 
                    !$omp shared(omp_chunk, Nq, dof, Nbasis) & 
                    !$omp firstprivate(Evec_q) &
                    !$omp private(eV_cmp1, eV_cmp2, mul_fac) &
                    !$omp firstprivate(W_q) &
                    !$omp private(q, Omega, nqs, BE_fac, R3_R4, q_dotR, expnt, mass_fac) &
                    !$omp private(nn, strt, ss, mu, Natm4, N2nu, fC2_N2, N3eta, N3indx, eta) &
                    !$omp private(N4ro, N4indx, ro, aa, bb, cc, dd, row_no1, row_no2) & 
                    !$omp private(IsZero) &
                    !$omp collapse(1) 

        q_sum: do nn = 1, Nq

            q = q_points(:, nn)

            strt = 1
            IsZero = all( dabs(q) < cmp_prec )
            if (IsZero) strt = 4

            Evec_q = Evec(:, :, nn)
            W_q = freq_new(:, nn)

            dof_sum: do ss = strt, dof

                Omega = W_q(ss)

                ne0_chk: if ( dabs(Omega) < zero_prec ) then
                    write(*, 255) q
                    write(*, 145) Omega

                else ne0_chk

                    neg_chk: if ( Omega < 0.0_dp ) then
                        write(*, 355) Omega, q

                    end if neg_chk

                    nqs = BEdist(Omega, T)
                    BE_fac = (2.0_dp * nqs + 1.0_dp)

                    mu_loop: do mu = 1, Nbasis
                        Natm4 = FC4%atmNum(mu)

                        N3_sum: do N3eta = 1, Natm4
                            N3indx = FC4%atmIndx(1:3, N3eta, mu)
                            eta = FC4%atmIndx(4, N3eta, mu)
                            
                            N4_sum: do N4ro = 1, Natm4
                                N4indx = FC4%atmIndx(1:3, N4ro, mu)
                                ro = FC4%atmIndx(4, N4ro, mu)

                                R3_R4 = matmul( sys%latvec, (N3indx-N4indx) )
                                q_dotR = dot_product( q, R3_R4 )
                                !expnt = cdexp( iu * q_dotR )
                                expnt = dcmplx( dcos(q_dotR), dsin(q_dotR) )

                                mass_fac = dsqrt( sys%mass(eta) * sys%mass(ro) )

                                gamma_sum: do cc = 1, 3
                                    row_no1 = 3*(eta-1) + (cc-1) + 1
                                    eV_cmp1 = Evec_q(row_no1, ss)


                                    delta_sum: do dd = 1, 3
                                        row_no2 = 3*(ro-1) + (dd-1) + 1
                                        eV_cmp2 = dconjg( Evec_q(row_no2, ss) )

                                        mul_fac = eV_cmp2 * eV_cmp1 * BE_fac * expnt / ( Omega * mass_fac )

                                        N2_loop: do N2nu = 1, Natm4
                                            fC2_N2 = FC4%atmIndx(5, N2nu, mu)

                                            alpha_loop: do aa = 1, 3
                                                beta_loop: do bb = 1, 3

                                                    !*!#!$omp atomic
                                                    dFC2(1, bb, aa, fC2_N2, mu) = &
                                                    dFC2(1, bb, aa, fC2_N2, mu) + &
                                                    FC4%fC(dd, cc, bb, aa, N4ro, N3eta, N2nu, mu) * mul_fac

                                                end do beta_loop
                                            end do alpha_loop
                                        end do N2_loop

                                    end do delta_sum
                                end do gamma_sum

                            end do N4_sum
                        end do N3_sum

                    end do mu_loop
                end if ne0_chk
            end do dof_sum
        end do q_sum

        !$omp end parallel do

        dFC2 = dFC2 * A2Conv / (4.0_dp * Nq)
        write(*, 855) maxval(dabs(dimag(dFC2)))

        FC2_new%fC = dble(dFC2)
        call AcousticSumRule(Nbasis, FC4, FC2_new)

        FC2_new = FC2_old + FC2_new

        freq_old = freq_new

        deallocate( freq_new, Evec )
        call get_freq_ev(sys, FC2_new, q_points, EwaldConst, LongEw, freq_new, Evec)

        freq_acc = maxval(dabs(freq_new-freq_old)) * THzConv
        write(*, 955) freq_acc

        exit_chk: if ( freq_acc < Renorm_accu ) then
            write(*, 333) iter_no

            write(*, 595)
            write(*, 975) iter_no
            write(*, *)

            exit scf_iter_loop

        else if ( iter_no > MAX_ITER ) then exit_chk
            write(*, 777) iter_no

            write(*, 595)
            write(*, 975) iter_no
            write(*, *)

            exit scf_iter_loop

        else exit_chk
            write(*, 595)
            write(*, 975) iter_no

            iter_no = iter_no+1
            write(*, *)

        end if exit_chk

    end do scf_iter_loop

    175 FORMAT('|', 37('='), ' Current iteration: ', &
             & I4, ' ', 37('='), '|')
    595 FORMAT('|', 99X, '|')

    255 FORMAT('|**', 16X, 'Zero frequency found at q = (', 3E11.4 ,')', 16X, '**|')
    145 FORMAT('|**', 19X, 'Frequency = ',E12.4, ' . Ignoring to avoid divergence.', 20X, '**|')

    355 FORMAT('|**', 4X, 'Warning: Imaginary frequency = ', E12.4, ', at q = (', 3E11.4 ,')', 4X, '**|')

    855 FORMAT('|', 12X, 'Maximum value of the imaginary part (ideally 0) of the IFC: ', E14.6, 13X, '|')

    955 FORMAT('|', 27X, 'Accuracy of the freuency: ', E14.6, ' THz', 28X, '|')

    975 FORMAT('|', 38('='), ' Iteration: ', I4, &
             & ' Ended ', 38('='), '|')

    333 FORMAT('|**', 16X, 'Convergence in Self Consistence loop achieved at iter no: ', I4, 17X, '**|')

    777 FORMAT('|', 23X, 'Convergence is not achieved, MAX_ITER exceeded: ', I4, 24X, '|')

end subroutine Renormalize


Function BEdist(w, T) RESULT(nqs)

    implicit none
    real(dp), intent(in)        :: w, T
    real(dp)                    :: nqs  !Result

    nqs = 1.0_dp / ( dexp(BEConv*w/T) - 1.0_dp )

end Function BEdist


subroutine AcousticSumRule(Nbasis, FC4, FC2)

    implicit none

    integer, intent(in)                     :: Nbasis
    type(FC4type), intent(in)               :: FC4
    type(FC2type), intent(inout)            :: FC2

    ! ================================== Local Variables ================================== !

    real(dp)                                :: sum_val, sub_term

    integer                                 :: mu, Natm4, alpha, beta, N4, N2

    ! ================================== Local Variables ================================== !

    mu_loop: do mu = 1, Nbasis
        Natm4 = FC4%atmNum(mu)

        alpha_loop: do alpha = 1, 3
            beta_loop: do beta = 1, 3

                sum_val = 0.0_dp

                ! --------------------------------------------------------------------------- !
                N2nu_loop1: do N4 = 1, Natm4

                    N2 = FC4%atmIndx(5, N4, mu)

                    sum_val = sum_val + FC2%fC(1, beta, alpha, N2, mu)

                end do N2nu_loop1

                sub_term = sum_val / dble( Natm4 )

                N2nu_loop2: do N4 = 1, Natm4

                    N2 = FC4%atmIndx(5, N4, mu)

                    FC2%fC(1, beta, alpha, N2, mu) = FC2%fC(1, beta, alpha, N2, mu) - sub_term

                end do N2nu_loop2
                ! --------------------------------------------------------------------------- !

            end do beta_loop
        end do alpha_loop

    end do mu_loop

end subroutine AcousticSumRule


end module Renorm

