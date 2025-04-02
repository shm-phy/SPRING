
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

module PointDefect

    use kinds,              only : dp
    use constants,          only : PI, wTHz
    use mklWrap,            only : zinverse
    use unit_cell,          only : cell
    use FC2_mod,            only : FC2type
    use phonon_m,           only : Phon
    use TetrahedronQ,       only : q_tetrahedron
    use GFIntegrationCoeff, only : CoefAtVertOfTetrhd, &
                                 & my_sort_GF, DiracDeltaCoeffVertices_1984

    implicit none

    PRIVATE
    PUBLIC      :: ScatterRatePointDefect, PerturbationMatrix

contains

    subroutine PerturbationMatrix( Ndof, Nbasis, M_defAtom, sys, FC2, FC2_per, V_pertbFC, V_pertbM )

        implicit none

        integer, intent(in)                             :: Ndof, Nbasis
        real(dp), dimension(Nbasis), intent(in)         :: M_defAtom
        type(cell), intent(in)                          :: sys
        type(FC2type), intent(in)                       :: FC2, FC2_per
        real(dp), dimension(Ndof, Ndof), intent(out)    :: V_pertbFC, V_pertbM

        !! ===================================== Local variables =================================== !!

        real(dp)                                :: ifc2_pristine, ifc2_perturbed

        integer, dimension(3)                   :: cell
        integer                                 :: nu, mu, beta, alpha, col_no, row_no

        !! ===================================== Local variables =================================== !!

        cell = (/0, 0, 0/)
        V_pertbFC(:, :) = 0.0_dp
        V_pertbM(:, :) = 0.0_dp

        nu_loop: do nu = 1, Nbasis
            beta_loop: do beta = 1, 3

                col_no = 3*(nu-1) + (beta-1) + 1

                mu_loop: do mu = 1, Nbasis
                    alpha_loop: do alpha = 1, 3

                        row_no = 3*(mu-1) + (alpha-1)  + 1

                        ifc2_pristine = IFC2( cell, mu, nu, alpha, beta, FC2 )
                        ifc2_perturbed = IFC2( cell, mu, nu, alpha, beta, FC2_per )

                        V_pertbFC(row_no, col_no) = (ifc2_perturbed - ifc2_pristine) / dsqrt(sys%mass(mu) * sys%mass(nu))

                        if ((mu == nu) .and. (alpha == beta)) V_pertbM(row_no, col_no) = (M_defAtom(mu)-sys%mass(mu)) / sys%mass(mu)

                    end do alpha_loop
                end do mu_loop
                
            end do beta_loop
        end do nu_loop

        V_pertbFC = V_pertbFC * (wTHz**2)

    end subroutine PerturbationMatrix


    Function IFC2( cell, mu, nu, alpha, beta, FC2 ) Result( fc )

        implicit none
        integer, dimension(3), intent(in)       :: cell
        integer, intent(in)                     :: mu, nu, alpha, beta
        type(FC2type), intent(in)               :: FC2

        real(dp)                                :: fc !Result

        !! ===================================== Local variables =================================== !!

        integer, dimension(3)       :: N2crt
        integer                     :: N2, Natm, my_nu

        !! ===================================== Local variables =================================== !!

        Natm = FC2%atmNum(mu)
        N2_loop: do N2=1, Natm

            N2crt = FC2%atmIndx(1:3, N2, mu)
            my_nu = FC2%atmIndx(4, N2, mu)

            if ( (my_nu == nu) .and. all( (N2crt-cell) == 0 ) ) then
                fc =  FC2%fC(1, beta, alpha, N2, mu) 
                EXIT N2_loop
            end if

        end do N2_loop

    end Function IFC2


    subroutine ScatterRatePointDefect( i0, Ndof, NTtrhdrn_pd, my_Qsize_pd, my_offsetq_pd, vol, rho_imp, &
                                     & V_pertbFC, V_pertbM, qPnts, qPnts_pd, ph_q0, ph_pd, &
                                     & MatInvFull, OptcTherm, tau_q0pd_inv )

        implicit none

        integer, intent(in)                             :: i0, Ndof, NTtrhdrn_pd, &
                                                         & my_Qsize_pd, my_offsetq_pd
        real(dp), intent(in)                            :: vol, rho_imp
        real(dp), dimension(Ndof, Ndof), intent(in)     :: V_pertbFC, V_pertbM
        type(q_tetrahedron), intent(in)                 :: qPnts, qPnts_pd
        type(Phon), intent(in)                          :: ph_q0, ph_pd
        logical, intent(in)                             :: MatInvFull, OptcTherm

        real(dp), dimension(Ndof), intent(out)          :: tau_q0pd_inv

        !! ===================================== Local variables =================================== !!

        character(len=128)                                          :: msg

        complex(dp), allocatable, dimension(:, :), codimension[:]   :: G
        complex(dp), dimension(:, :), allocatable                   :: I_VGinv, T
        complex(dp), dimension(Ndof)                                :: ev_q0s0, Txq0s0
        complex(dp)                                                 :: q0s0xTxq0s0

        !-!real(dp), allocatable, dimension(:), codimension[:]         :: tau_q0s0inv !gfrotran breaks when -finit-real=snan,
        !-!                                                                        !so define a rank=1, dim=1 instead of scalar
        real(dp), dimension(:, :), allocatable                      :: V
        real(dp), dimension(Ndof)                                   :: omega_q0, omega_q02
        real(dp)                                                    :: Wq0s02, tau_q0s0inv

        integer, dimension(3)                                       :: q0_int
        integer                                                     :: s0, strt, ii, istat

        logical                                                     :: not_invertible

        !! ===================================== Local variables =================================== !!

        tau_q0pd_inv(:) = 0.0_dp

        allocate( G(Ndof, Ndof)[*] )
        SYNC ALL

        allocate( I_VGinv(Ndof, Ndof), T(Ndof, Ndof) )
        allocate( V(Ndof, Ndof) )

        omega_q0 = ph_q0%omega(:, i0)
        omega_q02 = omega_q0 ** 2
        q0_int = qPnts%irr_q_pnt_int(:, i0)

        strt = 1
        if ( all(q0_int == 0) ) strt = 4

        Polarization0: do s0 = strt, Ndof

            Wq0s02 = omega_q02( s0 )
            ev_q0s0 = ph_q0%Evec(:, s0, i0)

            call GreensFunction( Ndof, NTtrhdrn_pd, my_Qsize_pd, my_offsetq_pd, Wq0s02, qPnts_pd, ph_pd, G )

            SYNC ALL
            call co_sum( G, stat=istat, ERRMSG=msg )
            if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

            V = V_pertbFC - ( Wq0s02 * V_pertbM )

            I_VGinv = matmul( V, G )

            FullMatInv: if ( MatInvFull ) then

                I_VGinv = -I_VGinv

                diagonal: do ii = 1, Ndof
                    I_VGinv(ii, ii) = 1.0_dp + I_VGinv(ii, ii)
                end do diagonal

                call zinverse( I_VGinv, not_invertible )

                checkInv: if ( not_invertible ) then

                    I_VGinv = matmul( V, G ) !Reuse variable
                    T = V + matmul( I_VGinv, V )

                else checkInv

                    T = matmul( I_VGinv, V )

                end if checkInv

            else FullMatInv

                T = V + matmul( I_VGinv, V ) !** (I - VG)^-1 x V ~ V + VGV **!

            end if FullMatInv

            Txq0s0 = matmul( T, ev_q0s0 )

            OpticalTheorem: if ( OptcTherm ) then

                !** Following with optical theorem **!
                q0s0xTxq0s0 = dot_product( ev_q0s0, Txq0s0 )
                tau_q0s0inv = -1.0_dp * DIMAG( q0s0xTxq0s0 ) / omega_q0( s0 )
                tau_q0pd_inv( s0 ) = rho_imp * vol * tau_q0s0inv
                !** Following with optical theorem **!

            else OpticalTheorem

                !** Following is for full Tetrahedron integration **!
                tau_q0s0inv = TetrahedronIntrg( Ndof, NTtrhdrn_pd, my_Qsize_pd, my_offsetq_pd, &
                                              & Wq0s02, Txq0s0, qPnts_pd, ph_pd ) 
                SYNC ALL
                call co_sum( tau_q0s0inv, stat=istat, ERRMSG=msg )
                if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
                tau_q0pd_inv( s0 ) = rho_imp * vol * tau_q0s0inv / omega_q0( s0 )
                !** Following is for full Tetrahedron integration **!

            end if OpticalTheorem

        end do Polarization0

        deallocate( G )
        deallocate( I_VGinv, T )
        deallocate( V )

    end subroutine ScatterRatePointDefect


    subroutine GreensFunction( Ndof, NTtrhdrn_pd, my_Qsize_pd, my_offsetq_pd, Wq0s02, qPnts_pd, ph_pd, G )

        implicit none
        integer, intent(in)                             :: Ndof, NTtrhdrn_pd, &
                                                         & my_Qsize_pd, my_offsetq_pd
        real(dp), intent(in)                            :: Wq0s02
        type(q_tetrahedron), intent(in)                 :: qPnts_pd
        type(Phon), intent(in)                          :: ph_pd

        complex(dp), dimension(Ndof, Ndof), intent(out) :: G

        !! ===================================== Local variables =================================== !!

        complex(dp), dimension(:, :), allocatable       :: F_v
        complex(dp), dimension(4)                       :: Coeff_v

        real(dp), dimension(4)                          :: W2AtVer_s

        integer, dimension(4)                           :: TetraVertQ_us, TetraVertQ_s
        integer                                         :: ii, nT, s, vv

        !! ===================================== Local variables =================================== !!

        G(:, :) = dcmplx( 0.0_dp, 0.0_dp )

        allocate( F_v(Ndof, Ndof) )

        TetrahedronLoop: do ii = 1, my_Qsize_pd

            nT = my_offsetq_pd + ii
            TetraVertQ_us = qPnts_pd%tetrahedrons(:, nT)

            Polarization: do s = 1, Ndof

                vertc_loop: do vv = 1, 4

                    W2AtVer_s( vv ) = ph_pd%omega( s, TetraVertQ_us(vv) ) !It will be eventually sorted, so _s

                end do vertc_loop

                TetraVertQ_s(:) = TetraVertQ_us(:) !Copy it, kepp the unsorted vertices untouched
                call CoefAtVertOfTetrhd( Wq0s02, W2AtVer_s, TetraVertQ_s, Coeff_v )

                SortedVertices: do vv = 1, 4

                    F_v(:, :) = ph_pd%SpcProjOp( :, :, s, TetraVertQ_s(vv) )
                    G = G  + ( Coeff_v(vv) * F_v )

                end do SortedVertices

            end do Polarization

        end do TetrahedronLoop

        G = G / NTtrhdrn_pd ! 1/NTtrhdrn_pd = 1/(6N) = v/Omega

        deallocate( F_v )

    end subroutine GreensFunction

    
    Function TetrahedronIntrg( Ndof, NTtrhdrn_pd, my_Qsize_pd, my_offsetq_pd, &
                             & Wq0s02, Txq0s0, qPnts_pd, ph_pd ) Result( tau_q0s0inv )

        implicit none
        integer, intent(in)                             :: Ndof, NTtrhdrn_pd, &
                                                         & my_Qsize_pd, my_offsetq_pd
        real(dp), intent(in)                            :: Wq0s02
        complex(dp), dimension(Ndof), intent(in)        :: Txq0s0
        type(q_tetrahedron), intent(in)                 :: qPnts_pd
        type(Phon), intent(in)                          :: ph_pd

        real(dp)                                        :: tau_q0s0inv ! Result

        !! ===================================== Local variables =================================== !!

        complex(dp), dimension(Ndof)    :: ev_qs
        complex(dp)                     :: qsxTxq0s0

        real(dp), dimension(2)          :: ReIm
        real(dp), dimension(4)          :: Fv_s, Coeff_v, W2AtVer_s

        integer, dimension(4)           :: TetraVertQ_us, TetraVertQ_s
        integer                         :: ii, nT, s, vv

        !! ===================================== Local variables =================================== !!

        tau_q0s0inv = 0.0_dp

        TetrahedronLoop: do ii = 1, my_Qsize_pd

            nT = my_offsetq_pd + ii
            TetraVertQ_us = qPnts_pd%tetrahedrons(:, nT)

            Polarization: do s = 1, Ndof

                vertc_loop: do vv = 1, 4

                    W2AtVer_s( vv ) = ph_pd%omega( s, TetraVertQ_us(vv) ) !It will be eventually sorted, so _s

                    ev_qs = ph_pd%Evec(:, s, TetraVertQ_us(vv) )
                    qsxTxq0s0 = dot_product( ev_qs, Txq0s0 )
                    ReIm = transfer( qsxTxq0s0, ReIm )
                    Fv_s(vv) = dot_product( ReIm, ReIm ) !It will be eventually sorted, so _s

                end do vertc_loop

                TetraVertQ_s(:) = TetraVertQ_us(:) !Copy it, kepp the unsorted vertices untouched

                call my_sort_GF( W2AtVer_s, TetraVertQ_s, F_vert=Fv_s )
                Coeff_v = DiracDeltaCoeffVertices_1984( W2AtVer_s, Wq0s02 )

                !tau_qOs0inv = tau_qOs0inv + dot_product( Coeff_v, Fv_s )
                SortedVertices: do vv = 1, 4
                    tau_q0s0inv = tau_q0s0inv  + ( Coeff_v(vv) * Fv_s(vv) )
                end do SortedVertices

            end do Polarization

        end do TetrahedronLoop

        tau_q0s0inv = PI * tau_q0s0inv / dble(NTtrhdrn_pd) ! 1/NTtrhdrn_pd = 1/(6N) = v/Omega

    end Function TetrahedronIntrg

end module PointDefect

