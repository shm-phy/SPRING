
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

module GFIntegrationCoeff

    use kinds,      only : dp
    use constants,  only : zero_prec, cmp_prec, PI

    implicit none

    EXTERNAL        :: dlasrt

    private
    public :: CoefAtVertOfTetrhd, my_sort_GF, DiracDeltaCoeffVertices_1984

contains

    subroutine CoefAtVertOfTetrhd( Omega_q0, Omega_vert, tetra_vert_q, Coeff_v, F_vert )

        implicit none

        real(dp), intent(in)                                :: Omega_q0
        real(dp), dimension(4), intent(inout)               :: Omega_vert   !* It will be sorted at exit *!
        integer, dimension(4), intent(inout)                :: tetra_vert_q !* It will be arranged according to sorted 
                                                                            !* Omega_vert at exit *!
        complex(dp), dimension(4), intent(out)              :: Coeff_v
        real(dp), dimension(4), intent(inout), OPTIONAL     :: F_vert  !* It will be arranged according to sorted
                                                                       !* Omega_vert at exit *!
        !! ===================================== Local variables =================================== !!

        real(dp), dimension(4)                      :: d_is, c_is
        integer                                     :: vv

        !! ===================================== Local variables =================================== !!

        if ( PRESENT(F_vert) ) then
            call my_sort_GF( Omega_vert, tetra_vert_q, F_vert )
        else
            call my_sort_GF( Omega_vert, tetra_vert_q )
        end if

        !* Real Part *!
        call GFcoeffv_RealPart( Omega_vert, Omega_q0, d_is )
        !* Real Part *!

        !* Imaginary Part *!
        c_is = DiracDeltaCoeffVertices_1984( Omega_vert, Omega_q0 )    !**!
        !c_is = DiracDeltaCoeffVertices_1979( Omega_vert, Omega_q0 )
        !* Imaginary Part *!

        !* Just to avoid NaN *!
        WHERE( ISNAN(d_is) ) d_is = 0.0_dp
        WHERE( ISNAN(c_is) ) c_is = 0.0_dp
        !* Just to avoid NaN *!

        do vv = 1, 4
            Coeff_v(vv) = dcmplx( d_is(vv), (-PI * c_is(vv)) )
        end do

    end subroutine CoefAtVertOfTetrhd


    subroutine my_sort_GF( Omega_vert, tetra_vert_q, F_vert )

        implicit none

        real(dp), dimension(4), intent(inout)                       :: Omega_vert
        integer, dimension(4), intent(inout)                        :: tetra_vert_q
        real(dp), dimension(4), intent(inout), OPTIONAL             :: F_vert

        !! ===================================== Local variables =================================== !!

        real(dp), dimension(4)                      :: Omega_vertc
        real(dp), dimension(4)                      :: F_vertc

        integer, dimension(4)                       :: tetra_vert_qc
        integer, dimension(1)                       :: indx
        integer, dimension(1)                       :: indx_back
        integer                                     :: info, ii, strt_pnt, ActualIndx

        logical, dimension(4)                       :: covered
        logical                                     :: back_chk

        !! ===================================== Local variables =================================== !!

        covered(:) = .false.
        back_chk = .false.

        Omega_vertc = Omega_vert
        tetra_vert_qc = tetra_vert_q
        if ( PRESENT(F_vert) ) F_vertc = F_vert

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

            tetra_vert_q(ii) = tetra_vert_qc( ActualIndx )
            if ( PRESENT(F_vert) ) F_vert(ii) = F_vertc( ActualIndx )

            debug: if ( back_chk ) then

                indx_back = findloc(Omega_vertc, Omega_vert(ii), back=back_chk)

                if ( ActualIndx /= indx_back(1) ) then

                    write(*, 80)
                    write(*, 90) Omega_vert

                    80 FORMAT("WARNING: Omega(q) is same in two vertices of tetrahedron")
                    90 FORMAT("         The values of sorted Omega(q) at vertices: [", 3(F12.5, ', '), F12.5, "] THz")

                end if

            end if debug

        end do vert_loop

    end subroutine my_sort_GF


    subroutine GFcoeffv_RealPart( Omega_vert, Omega_q0, d_is )

        implicit none

        real(dp), dimension(4), intent(in)          :: Omega_vert
        real(dp), intent(in)                        :: Omega_q0
        real(dp), dimension(4), intent(out)         :: d_is

        !! ===================================== Local variables =================================== !!

        real(dp)                    :: vert_el

        integer                     :: i, j, same_count, typeOfequal
        integer, dimension(2, 3)    :: same_pairs
        integer, dimension(2)       :: SameIndx

        logical, dimension(4)       :: check_vert
        logical                     :: only_2_same, first_time_en

        !! ===================================== Local variables =================================== !!

        check_vert(:) = .true.
        only_2_same = .false.

        same_count = 0
        typeOfequal = 0
        
        same_pairs(:, :) = 0

        vert_loop: do i = 1, 4

            vert_el = Omega_vert( i )
            first_time_en = .true.

            CheckToEnter: if ( check_vert(i) ) then

                loop2: do j = i+1, 4

                    Check: if ( dabs(vert_el - Omega_vert(j)) < zero_prec ) then

                        same_count = same_count + 1

                        if ( first_time_en ) then
                            first_time_en = .false.
                            typeOfequal = typeOfequal + 1
                        end if 

                        same_pairs(:, same_count) = (/i, j/)

                        check_vert(j) = .false.

                    end if Check

                end do loop2

            end if CheckToEnter

        end do vert_loop

        EnterCondition: if ( same_count == 0 ) then

            !d_is = CoeffAtVerticesTypeD_1984( Omega_vert, Omega_q0 )
            d_is = CoeffAtVerticesTypeD_2004( Omega_vert, Omega_q0 )    !*!

        else if ( same_count == 1 ) then EnterCondition

            SameIndx(:) = same_pairs(:, 1)
            d_is = CoeffAtVerticesType0_1984( Omega_vert, Omega_q0, SameIndx )    !*!
            !d_is = CoeffAtVerticesType0_2004( Omega_vert, Omega_q0, SameIndx )

        else if ( same_count == 2 ) then EnterCondition

            SameIndx(:) = same_pairs(1, 1:2) !Just a temporary array to check

            TypeOfEquals: if ( typeOfequal == 1 ) then

                InnerIf: if ( all((SameIndx - 1) == 0) ) then

                    !d_is = CoeffAtVerticesType1_1984( Omega_vert, Omega_q0 )
                    d_is = CoeffAtVerticesType1_2004( Omega_vert, Omega_q0 )    !*!

                else if ( all((SameIndx - 2) == 0) ) then InnerIf

                    !d_is = CoeffAtVerticesType3_1984( Omega_vert, Omega_q0 )
                    d_is = CoeffAtVerticesType3_2004( Omega_vert, Omega_q0 )    !*!

                else InnerIf

                    write(*, 20) SameIndx
                    20 FORMAT( "ERROR (in GFcoeffv_RealPart): SameIndx can not other than 1 or 2 (", I2, ", ", I2, " )" )

                end if InnerIf

            else if ( typeOfequal == 2 ) then TypeOfEquals
                
                !d_is = CoeffAtVerticesType2_1984( Omega_vert, Omega_q0 )
                d_is = CoeffAtVerticesType2_2004( Omega_vert, Omega_q0 )    !*!

            else TypeOfEquals

                write(*, 40) typeOfequal
                40 FORMAT( "ERROR (in GFcoeffv_RealPart):  typeOfequal can not be > 2 ", I3 )

            end if TypeOfEquals

        else if ( same_count == 3 ) then EnterCondition

            d_is(:) = 0.25_dp / ( Omega_q0 - Omega_vert(1) )

        else EnterCondition

            write(*, 80) same_count
            80 FORMAT( "ERROR (in GFcoeffv_RealPart): same_count can not be > 3 ", I3 ) 

        end if EnterCondition

    end subroutine GFcoeffv_RealPart


    Pure Function CoeffAtVerticesType0_1984( Wv, W0, SameIndx ) Result( d_is )

        implicit none

        real(dp), dimension(4), intent(in)          :: Wv
        real(dp), intent(in)                        :: W0
        integer, dimension(2), intent(in)           :: SameIndx

        real(dp), dimension(4)                      :: d_is !Result

        !! ===================================== Local variables =================================== !!

        real(dp)                        :: d_l_m, T_jm, T_km

        integer                         :: v, j, k, l, m, c, i
        integer, dimension(2)           :: DiffIndx

        !! ===================================== Local variables =================================== !!

        d_is(:) = 0.0_dp

        c = 0

        l = SameIndx(1)
        m = SameIndx(2)

        ! j and k (j /= k) are two sites different from l and m
        do v = 1, 4
            if ( (v /= l) .and. (v /= m) ) then
                c = c + 1
                DiffIndx(c) = v
            end if
        end do
        j = DiffIndx(1)
        k = DiffIndx(2)

        T_jm = (W0 - Wv(j)) / (Wv(m)-Wv(j))
        T_km = (W0 - Wv(k)) / (Wv(m)-Wv(k))

        d_l_m = ( ((W0 - Wv(j))**3) / ((Wv(j) - Wv(k)) * ((Wv(j) - Wv(m))**3)) ) * dlog( dabs(W0 - Wv(j)) ) + &
              & ( ((W0 - Wv(k))**3) / ((Wv(k) - Wv(j)) * ((Wv(k) - Wv(m))**3)) ) * dlog( dabs(W0 - Wv(k)) ) + &
              &  ( (W0 - Wv(m))     / ((Wv(m) - Wv(j)) *  (Wv(m) - Wv(k))) ) * &
              &  ( 0.5_dp + T_jm + T_km + ( T_jm**2 + T_km**2 + T_jm*T_km ) *  dlog( dabs(W0 - Wv(m)) ) )

        d_is(l) = d_l_m
        d_is(m) = d_l_m

        do v = 1, 2

            i = DiffIndx(v)
            c = MOD(v, 2) + 1 !! MODULO(v, 2) + 1 !R such that v=Q*2+R, where Q is an integer and R is between 0 (inclusive) and 2 (exclusive)
            l = DiffIndx(c) ! Reuse l

            d_is(i) = ( ((W0 - Wv(i))**2) / (((Wv(i)-Wv(m))**2) * (Wv(l) - Wv(i))) ) * ( 1.0_dp + &
                    & ( (2.0_dp*(W0-Wv(m))/(Wv(i)-Wv(m))) + ((W0-Wv(l))/(Wv(i)-Wv(l))) ) * dlog( dabs(W0 - Wv(i)) ) ) + &
                    & ( ((W0 - Wv(m))**2) / (((Wv(m)-Wv(i))**2) * (Wv(l) - Wv(m))) ) * ( 1.0_dp + &
                    & ( (2.0_dp*(W0-Wv(i))/(Wv(m)-Wv(i))) + ((W0-Wv(l))/(Wv(m)-Wv(l))) ) * dlog( dabs(W0 - Wv(m)) ) ) + &
                    & ( ((W0-Wv(l))**3) / (((Wv(l)-Wv(i))**2) * ((Wv(l)-Wv(m))**2)) ) * dlog( dabs(W0 - Wv(l)) ) 

        end do 

    end Function CoeffAtVerticesType0_1984


    Pure Function CoeffAtVerticesType1_2004( Wv, W0 ) Result( d_is )

        implicit none

        real(dp), dimension(4), intent(in)          :: Wv
        real(dp), intent(in)                        :: W0

        real(dp), dimension(4)                      :: d_is !Result

        !! ===================================== Local variables =================================== !!

        real(dp)                        :: q_123, g1, g2, g3, g4, &
                                         & g4_1, g1_4

        !! ===================================== Local variables =================================== !!

        d_is(:) = 0.0_dp

        g1 = W0 - Wv(1)
        g2 = W0 - Wv(2)
        g3 = W0 - Wv(3)
        g4 = W0 - Wv(4)

        g4_1 = g4 / g1
        g1_4 = g1 / g4

        q_123 = 11.0_dp / (6.0_dp * ( 1.0_dp - g4_1)) - 5.0_dp / ( 2.0_dp * ((1.0_dp - g4_1)**2) ) + &
              & 1.0_dp / ( (1.0_dp - g4_1)**3 ) + dlog( dabs(g4_1) ) / ( ((1.0_dp - g1_4)**4) * g4_1 )
        d_is(1:3) = q_123 / g1

        d_is(4) = ( 3.0_dp / ((1.0_dp - g4_1)**2) - 6.0_dp / ((1.0_dp - g4_1)**3) + 3.0_dp / ((1.0_dp - g4_1)**4) ) * &
                & ( dlog( dabs(g1_4) ) / g1_4 ) + 1.0_dp / ((1.0_dp - g1_4)**3) + 2.5 * g4_1 / ((1.0_dp - g4_1)**2) - &
                & 2.0_dp * g4_1 / ((1.0_dp - g4_1)**3)
        d_is(4) = d_is(4) / g4

    end Function CoeffAtVerticesType1_2004


    Pure Function CoeffAtVerticesType3_2004( Wv, W0 ) Result( d_is )

        implicit none

        real(dp), dimension(4), intent(in)          :: Wv
        real(dp), intent(in)                        :: W0

        real(dp), dimension(4)                      :: d_is !Result

        !! ===================================== Local variables =================================== !!

        real(dp)                        :: q_234, g1, g2, g3, g4, &
                                         & g1_2, g2_1

        !! ===================================== Local variables =================================== !!

        d_is(:) = 0.0_dp

        g1 = W0 - Wv(1)
        g2 = W0 - Wv(2)
        g3 = W0 - Wv(3)
        g4 = W0 - Wv(4)

        g1_2 = g1 / g2
        g2_1 = g2 / g1

        q_234 = 11.0_dp / (6.0_dp * ( 1.0_dp - g1_2)) - 5.0_dp / ( 2.0_dp * ((1.0_dp - g1_2)**2) ) + &
              & 1.0_dp / ( (1.0_dp - g1_2)**3 ) + dlog( dabs(g1_2) ) / ( ((1.0_dp - g2_1)**4) * g1_2 )
        d_is(2:4) = q_234 / g2

        d_is(1) = ( 3.0_dp / ((1.0_dp - g1_2)**2) - 6.0_dp / ((1.0_dp - g1_2)**3) + 3.0_dp / ((1.0_dp - g1_2)**4) ) * &
                & ( dlog( dabs(g2_1) ) / g2_1 ) + 1.0_dp / ((1.0_dp - g2_1)**3) + 2.5 * g1_2 / ((1.0_dp - g1_2)**2) - &
                & 2.0_dp * g1_2 / ((1.0_dp - g1_2)**3)
        d_is(1) = d_is(1) / g1

    end Function CoeffAtVerticesType3_2004


    Pure Function CoeffAtVerticesType2_2004( Wv, W0 ) Result( d_is )

        implicit none

        real(dp), dimension(4), intent(in)          :: Wv
        real(dp), intent(in)                        :: W0

        real(dp), dimension(4)                      :: d_is !Result

        !! ===================================== Local variables =================================== !!

        real(dp)                        :: q_12, q_34, g1, g2, g3, g4, &
                                         & g3_1, g1_3

        !! ===================================== Local variables =================================== !!

        d_is(:) = 0.0_dp

        g1 = W0 - Wv(1)
        g2 = W0 - Wv(2)
        g3 = W0 - Wv(3)
        g4 = W0 - Wv(4)

        g3_1 = g3 / g1
        g1_3 = g1 / g3

        q_12 = 2.5_dp / ( (1.0_dp-g3_1)**2 ) - 2.0_dp / ( (1.0_dp-g3_1)**3 ) + g1_3 / ( (1.0_dp-g1_3)**3 ) + &
             & 3.0_dp * ( 1.0_dp / ((1.0_dp-g1_3)**3) - 1.0_dp / ((1.0_dp-g1_3)**4) ) * ( dlog( dabs(g3_1) ) / g3_1 )

        q_34 = 2.5_dp / ( (1.0_dp-g1_3)**2 ) - 2.0_dp / ( (1.0_dp-g1_3)**3 ) + g3_1 / ( (1.0_dp-g3_1)**3 ) + &
             & 3.0_dp * ( 1.0_dp / ((1.0_dp-g3_1)**3) - 1.0_dp / ((1.0_dp-g3_1)**4) ) * ( dlog( dabs(g1_3) ) / g1_3 )

        d_is(1:2) = q_12 / g1
        d_is(3:4) = q_34 / g3

    end Function CoeffAtVerticesType2_2004


    Pure Function CoeffAtVerticesTypeD_2004( Wv, W0 ) Result( d_is )

        implicit none

        real(dp), dimension(4), intent(in)          :: Wv
        real(dp), intent(in)                        :: W0

        real(dp), dimension(4)                      :: d_is !Result

        !! ===================================== Local variables =================================== !!

        real(dp)                        :: g1, g2, g3, g4, q
                                         
        integer, dimension(4)           :: cycl_v
        integer                         :: v, vc, next_indx

        !! ===================================== Local variables =================================== !!

        d_is(:) = 0.0_dp

        do v = 1, 4

            cycl_v(1) = v
            do vc = 2, 4
                next_indx = v + (vc - 1) - 1
                next_indx = MOD(next_indx, 4) + 1 !! MODULO(next_indx, 4) + 1 !R such that next_indx=Q*2+R, where Q is an integer and R is
                                                  !! between 0 (inclusive) and 4 (exclusive)
                cycl_v(vc) = next_indx
            end do

            g1 = W0 - Wv( cycl_v(1) )
            g2 = W0 - Wv( cycl_v(2) )
            g3 = W0 - Wv( cycl_v(3) )
            g4 = W0 - Wv( cycl_v(4) )

            q = 1.0_dp / ( (1.0_dp - g2/g1) * (1.0_dp - g3/g1) * (1.0_dp - g4/g1) ) + &
            & ( 1.0_dp / ( ((1.0_dp - g1/g2)**2) * (1.0_dp - g3/g2) * (1.0_dp - g4/g2) ) ) * ( dlog(dabs(g2/g1)) / (g2/g1) ) + &
            & ( 1.0_dp / ( ((1.0_dp - g1/g3)**2) * (1.0_dp - g2/g3) * (1.0_dp - g4/g3) ) ) * ( dlog(dabs(g3/g1)) / (g3/g1) ) + &
            & ( 1.0_dp / ( ((1.0_dp - g1/g4)**2) * (1.0_dp - g2/g4) * (1.0_dp - g3/g4) ) ) * ( dlog(dabs(g4/g1)) / (g4/g1) )

            d_is(v) = q / g1

            !*!write(*, 67) cycl_v
            !*!67 FORMAT( "cycl_v : (", 3(I2, ' '), I2, " )" )
        end do

    end Function CoeffAtVerticesTypeD_2004


    Function DiracDeltaCoeffVertices_1984( Omega_vert, E ) Result( c_is )

        implicit none

        real(dp), dimension(4), intent(in)          :: Omega_vert !It must be sorted
        real(dp), intent(in)                        :: E

        real(dp), dimension(4)                      :: c_is !Result

        !! ===================================== Local variables =================================== !!

        real(dp)                :: E1, E2, E3, E4
        logical                 :: enterCase

        !! ===================================== Local variables =================================== !!

        E1 = Omega_vert(1)
        E2 = Omega_vert(2)
        E3 = Omega_vert(3)
        E4 = Omega_vert(4)

        enterCase = .false.

        ConditionOnE: if ( ( (E1 - E) <= cmp_prec ) .and. ( ( E - E2 ) <= cmp_prec ) ) then

            c_is(1) = ( (E2-E)/(E2-E1) + (E3-E)/(E3-E1) + (E4-E)/(E4-E1) ) * ( (E-E1)**2 / ((E4-E1)*(E3-E1)*(E2-E1)) )
            c_is(2) = (E-E1)**3 / ( (E2-E1)**2 * (E3-E1) * (E4-E1) )
            c_is(3) = (E-E1)**3 / ( (E2-E1) * (E3-E1)**2 * (E4-E1) )
            c_is(4) = (E-E1)**3 / ( (E2-E1) * (E3-E1) * (E4-E1)**2 )

            enterCase = .true.

        else if ( ( (E2 - E) <= cmp_prec ) .and. ( (E - E3) <= cmp_prec ) ) then ConditionOnE

            c_is(1) = ((E3 - E) / (E3 - E1)**2) * ( &
          & ( ((E3-E)*(E-E2)) / ((E4-E2)*(E3-E2)) ) + ( ((E4-E)*(E-E1)) / ((E4-E1)*(E4-E2)) ) + &
          & ( ((E3-E)*(E-E1)) / ((E3-E2)*(E4-E1)) ) ) + ((E4 - E) / (E4 - E1)**2) * ( &
          & ( ((E4-E)*(E-E1)) / ((E4-E2)*(E3-E1)) ) + ( ((E4-E)*(E-E2)) / ((E4-E2)*(E3-E2)) ) + &
          & ( ((E3-E)*(E-E1)) / ((E3-E1)*(E3-E2)) ) )

            c_is(2) = ((E3 - E) / (E3 - E2)**2) * ( &
          & ( ((E3-E)*(E-E2)) / ((E4-E2)*(E3-E1)) ) + ( ((E4-E)*(E-E2)) / ((E4-E2)*(E4-E1)) ) + &
          & ( ((E3-E)*(E-E1)) / ((E3-E1)*(E4-E1)) ) ) + ((E4 - E) / (E4 - E2)**2) * ( &
          & ( ((E3-E)*(E-E2)) / ((E3-E2)*(E3-E1)) ) + ( ((E4-E)*(E-E1)) / ((E4-E1)*(E3-E1)) ) + &
          & ( ((E4-E)*(E-E2)) / ((E3-E2)*(E4-E1)) ) )

            c_is(3) = ((E - E2) / (E3 - E2)**2) * ( &
          & ( ((E3-E)*(E-E2)) / ((E4-E2)*(E3-E1)) ) + ( ((E4-E)*(E-E2)) / ((E4-E2)*(E4-E1)) ) + &
          & ( ((E3-E)*(E-E1)) / ((E3-E1)*(E4-E1)) ) ) + ((E - E1) / (E3 - E1)**2) * ( &
          & ( ((E3-E)*(E-E2)) / ((E4-E2)*(E3-E2)) ) + ( ((E4-E)*(E-E1)) / ((E4-E1)*(E4-E2)) ) + &
          & ( ((E3-E)*(E-E1)) / ((E3-E2)*(E4-E1)) ) )
            
            c_is(4) = ((E - E2) / (E4 - E2)**2) * ( &
          & ( ((E3-E)*(E-E2)) / ((E3-E2)*(E3-E1)) ) + ( ((E4-E)*(E-E1)) / ((E4-E1)*(E3-E1)) ) + &
          & ( ((E4-E)*(E-E2)) / ((E3-E2)*(E4-E1)) ) ) + ((E - E1) / (E4 - E1)**2) * ( &
          & ( ((E4-E)*(E-E1)) / ((E4-E2)*(E3-E1)) ) + ( ((E4-E)*(E-E2)) / ((E4-E2)*(E3-E2)) ) + &
          & ( ((E3-E)*(E-E1)) / ((E3-E1)*(E3-E2)) ) )

            c_is = 0.5_dp * c_is
            enterCase = .true.

        else if ( ( (E3 - E) <= cmp_prec ) .and. ( (E - E4) <= cmp_prec ) ) then ConditionOnE

            c_is(1) = (E4 - E)**3 / ( (E4 - E1)**2 * (E4 - E2) * (E4 - E3) )
            c_is(2) = (E4 - E)**3 / ( (E4 - E1) * (E4 - E2)**2 * (E4 - E3) )
            c_is(3) = (E4 - E)**3 / ( (E4 - E1) * (E4 - E2) * (E4 - E3)**2 )
            c_is(4) = ( (E-E3)/(E4-E3) + (E-E2)/(E4-E2) + (E-E1)/(E4-E1) ) * &
                  & ( (E4 - E)**2 / ( (E4 - E1) * (E4 - E2) * (E4 - E3) ) )

            enterCase = .true.

        else ConditionOnE

            does_not_enters: if ( (E >= E1) .and. (E <= E4) .and. (.not. enterCase) ) then
                
                write(*, 55) E
                write(*, 65) Omega_vert

                55 FORMAT("WARNING: E0 = ", F12.5, "THz. It does not enters any cases.")
                65 FORMAT("         The values of E at vertices: [", 3(F12.5, ', '), F12.5, "] THz")

            end if does_not_enters

            c_is = 0.0_dp
            enterCase = .false.

        end if ConditionOnE

    end Function DiracDeltaCoeffVertices_1984

    ! **** Following Procedures were created just to check the correctness **** !

    !-! Pure Function CoeffAtVerticesType0_2004( Wv, W0, SameIndx ) Result( d_is )

    !-!     implicit none

    !-!     real(dp), dimension(4), intent(in)          :: Wv
    !-!     real(dp), intent(in)                        :: W0
    !-!     integer, dimension(2), intent(in)           :: SameIndx

    !-!     real(dp), dimension(4)                      :: d_is !Result

    !-!     !! ===================================== Local variables =================================== !!

    !-!     real(dp), dimension(4)          :: g
    !-!     real(dp)                        :: q_l_m

    !-!     integer                         :: v, j, k, l, m, c, i
    !-!     integer, dimension(2)           :: DiffIndx

    !-!     !! ===================================== Local variables =================================== !!

    !-!     d_is(:) = 0.0_dp
    !-!     c = 0
    !-!     g(1) = W0 - Wv(1)
    !-!     g(2) = W0 - Wv(2)
    !-!     g(3) = W0 - Wv(3)
    !-!     g(4) = W0 - Wv(4)

    !-!     l = SameIndx(1)
    !-!     m = SameIndx(2)
    !-!     ! j and k (j /= k) are two sites different from l and m
    !-!     do v = 1, 4
    !-!         if ( (v /= l) .and. (v /= m) ) then
    !-!             c = c + 1
    !-!             DiffIndx(c) = v
    !-!         end if
    !-!     end do
    !-!     j = DiffIndx(1)
    !-!     k = DiffIndx(2)

    !-!     q_l_m = 5.0_dp / ( 2.0_dp * (1.0_dp - (g(j)/g(m))) * (1.0_dp - (g(k)/g(m))) ) - &
    !-!           & 1.0_dp / ( ((1.0_dp - (g(j)/g(m)))**2) * (1.0_dp - (g(k)/g(m))) )  - &
    !-!           & 1.0_dp / ( (1.0_dp - (g(j)/g(m))) * ((1.0_dp - (g(k)/g(m)))**2) ) + &
    !-!           & (1.0_dp / ( ((1.0_dp - (g(m)/g(j)))**3) * (1.0_dp - (g(k)/g(j))) )) * &
    !-!           & ( dlog(dabs(g(j)/g(m))) / (g(j) / g(m)) ) + &
    !-!           & (1.0_dp / ( ((1.0_dp - (g(m)/g(k)))**3) * (1.0_dp - (g(j)/g(k))) )) * &
    !-!           & ( dlog(dabs(g(k)/g(m))) / (g(k) / g(m)) )

    !-!     d_is(l) = q_l_m / g(l)
    !-!     d_is(m) = q_l_m / g(m)

    !-!     do v = 1, 2

    !-!         i = DiffIndx(v)
    !-!         c = MOD(v, 2) + 1 !! MODULO(v, 2) + 1 !R such that v=Q*2+R, where Q is an integer and R is between 0 (inclusive) and 2 (exclusive)
    !-!         l = DiffIndx(c) ! Reuse l

    !-!         d_is(i) = 1.0_dp / ( ((1.0_dp - (g(m)/g(i)))**2) * (1.0_dp - (g(l)/g(i))) ) + &
    !-!            & (g(i)/g(m)) / ( ((1.0_dp - (g(i)/g(m)))**2) * (1.0_dp - (g(l)/g(m))) ) + &
    !-!            & ( dlog(dabs(g(l)/g(i))) / (g(l) / g(i)) ) / &
    !-!                          & ( ((1.0_dp - (g(m)/g(l)))**2) * ((1.0_dp - (g(i)/g(l)))**2 ) ) + ( &
    !-!            & ( 3.0_dp / ( ((1.0_dp - (g(i)/g(m)))**2) * (1.0_dp - (g(l)/g(m))) ) ) - &
    !-!            & ( 1.0_dp / ( ((1.0_dp - (g(i)/g(m)))**2) * ((1.0_dp - (g(l)/g(m)))**2) ) ) - &
    !-!            & ( 2.0_dp / ( ((1.0_dp - (g(i)/g(m)))**3) * (1.0_dp - (g(l)/g(m))) ) ) ) * &
    !-!            & ( dlog(dabs(g(m)/g(i))) / (g(m)/g(i)) ) 

    !-!         d_is(i) = d_is(i) / g(i)

    !-!     end do

    !-! end Function CoeffAtVerticesType0_2004


    !-! Pure Function CoeffAtVerticesType1_1984( Wv, W0 ) Result( d_is )

    !-!     implicit none

    !-!     real(dp), dimension(4), intent(in)          :: Wv
    !-!     real(dp), intent(in)                        :: W0

    !-!     real(dp), dimension(4)                      :: d_is !Result

    !-!     !! ===================================== Local variables =================================== !!

    !-!     real(dp)                        :: d_123, EE4, EE3, E4E3

    !-!     !! ===================================== Local variables =================================== !!

    !-!     d_is(:) = 0.0_dp

    !-!     EE4 = W0 - Wv(4)
    !-!     EE3 = W0-Wv(3)
    !-!     E4E3 = Wv(4) - Wv(3)
    !-! 
    !-!     d_123 = ( (EE4**3) / (E4E3**4) ) * dlog( dabs(EE4 / EE3) ) + &
    !-!           & ( 6.0_dp * (EE4**2) - 3.0_dp * E4E3 * EE4 + 2.0_dp * (E4E3**2) ) / ( 6.0_dp * (E4E3**3) )

    !-!     d_is(1:3) = d_123

    !-!     d_is(4) = ( 3.0_dp * (EE4**2) * EE3 / (E4E3**4) ) * dlog( dabs(EE3 / EE4) ) - &
    !-!             & 1.5_dp * EE3 * ( (2.0_dp * EE4 - E4E3) / (E4E3**3) ) - (1.0_dp / E4E3)

    !-! end Function CoeffAtVerticesType1_1984


    !-! Pure Function CoeffAtVerticesType3_1984( Wv, W0 ) Result( d_is )

    !-!     implicit none

    !-!     real(dp), dimension(4), intent(in)          :: Wv
    !-!     real(dp), intent(in)                        :: W0

    !-!     real(dp), dimension(4)                      :: d_is !Result

    !-!     !! ===================================== Local variables =================================== !!

    !-!     real(dp)                        :: d_234, EE1, EE2, E2E1

    !-!     !! ===================================== Local variables =================================== !!

    !-!     d_is(:) = 0.0_dp

    !-!     EE1 = W0 - Wv(1)
    !-!     EE2 = W0-Wv(2)
    !-!     E2E1 = Wv(2) - Wv(1)
    !-! 
    !-!     d_234 = ( (EE1**3) / (E2E1**4) ) * dlog( dabs(EE1 / EE2) ) - &
    !-!           & ( 6.0_dp * (EE1**2) + 3.0_dp * EE1 * E2E1 + 2.0_dp * (E2E1**2) ) / ( 6.0_dp * (E2E1**3) )

    !-!     d_is(2:4) = d_234

    !-!     d_is(1) = ( 3.0_dp * (EE1**2) * EE2 / (E2E1**4) ) * dlog( dabs(EE2 / EE1) ) + &
    !-!             & 1.5_dp * EE2 * ( (2.0_dp * EE1 + E2E1) / (E2E1**3) ) + (1.0_dp / E2E1)

    !-! end Function CoeffAtVerticesType3_1984


    !-! Pure Function CoeffAtVerticesType2_1984( Wv, W0 ) Result( d_is )

    !-!     implicit none

    !-!     real(dp), dimension(4), intent(in)          :: Wv
    !-!     real(dp), intent(in)                        :: W0

    !-!     real(dp), dimension(4)                      :: d_is !Result

    !-!     !! ===================================== Local variables =================================== !!

    !-!     real(dp)                        :: d_12, d_34, EE3, EE2, E3E2

    !-!     !! ===================================== Local variables =================================== !!

    !-!     d_is(:) = 0.0_dp

    !-!     EE3 = W0 - Wv(3)
    !-!     EE2 = W0 - Wv(2)
    !-!     E3E2 = Wv(3) - Wv(2)

    !-!     d_12 = 3.0_dp * ( (EE3**2) * EE2 / (E3E2**4) ) * dlog( dabs(EE2 / EE3) ) - &
    !-!          & 1.5_dp * EE2 * ( (2.0_dp * EE3 - E3E2) / (E3E2**3) ) - (1.0_dp / E3E2)

    !-!     d_34 = 3.0_dp * ( (EE2**2) * EE3 / (E3E2**4) ) * dlog( dabs(EE3 / EE2) ) + &
    !-!          & 1.5_dp * EE3 * ( (2.0_dp * EE2 + E3E2) / (E3E2**3) ) + (1.0_dp / E3E2)

    !-!     d_is(1:2) = d_12
    !-!     d_is(3:4) = d_34

    !-! end Function CoeffAtVerticesType2_1984


    !-! Pure Function CoeffAtVerticesTypeD_1984( Wv, W0 ) Result( d_is )

    !-!     implicit none

    !-!     real(dp), dimension(4), intent(in)          :: Wv
    !-!     real(dp), intent(in)                        :: W0

    !-!     real(dp), dimension(4)                      :: d_is !Result

    !-!     !! ===================================== Local variables =================================== !!

    !-!     real(dp)                        :: denom1, denom2, term1, term2, d

    !-!     integer                         :: i, j, k

    !-!     !! ===================================== Local variables =================================== !!

    !-!     d_is(:) = 0.0_dp

    !-!     VerticesLoop: do i = 1, 4

    !-!         denom1 = 1.0_dp
    !-!         do k = 1, 4
    !-!             if ( k /= i ) denom1 = denom1 * (Wv(k) - Wv(i)) 
    !-!         end do

    !-!         term1 = 0.0_dp
    !-!         do j = 1, 4
    !-!             if ( j /= i ) term1 = term1 + ( (W0 - Wv(j)) / (Wv(i) - Wv(j)) ) * dlog( dabs(W0 - Wv(i)) )
    !-!         end do

    !-!         term2 = 0.0_dp
    !-!         do j = 1, 4

    !-!             if ( j /= i ) then

    !-!                 denom2 = 1.0_dp
    !-!                 do k = 1, 4
    !-!                     if ( k /= j ) denom2 = denom2 * (Wv(k) - Wv(j))
    !-!                 end do

    !-!                 term2 = term2 + ( ((W0 - Wv(j))**3) / denom2 ) * dlog( dabs(W0 - Wv(j)) ) / (Wv(i) - Wv(j))

    !-!             end if

    !-!         end do

    !-!         d = ( ((W0 - Wv(i))**2) / denom1 ) * ( 1.0_dp + term1 ) + term2

    !-!         d_is(i) = d

    !-!     end do VerticesLoop

    !-! end Function CoeffAtVerticesTypeD_1984


    !-! Function DiracDeltaCoeffVertices_1979( E, E0 ) Result( c_is )

    !-!     implicit none

    !-!     real(dp), dimension(4), intent(in)          :: E !It must be sorted
    !-!     real(dp), intent(in)                        :: E0

    !-!     real(dp), dimension(4)                      :: c_is !Result
    !-!     
    !-!     ! ====================================== Local Variables ====================================== !

    !-!     real(dp)                                :: g
    !-!     real(dp), dimension(4)                  :: Iver

    !-!     integer                                 :: vv
    !-!     logical                                 :: enterCase

    !-!     ! ====================================== Local Variables ====================================== !

    !-!     enterCase = .false.

    !-!     tetra_coef_case: if ( ( (E(1) - E0) <= cmp_prec ) .and. ( (E0 - E(2)) <= cmp_prec ) ) then

    !-!         ! ------------------- value of g and values of I at vertices ------------------- !
    !-!         g = ( 3.0_dp * ( (E0 - E(1)) ** 2 ) ) / &
    !-!           & ( (E(2) - E(1)) * & 
    !-!           &   (E(3) - E(1)) * &
    !-!           &   (E(4) - E(1)) )

    !-!         Iver(1) = ( ( (E0 - E(2)) / (E(1) - E(2)) ) + &
    !-!                   & ( (E0 - E(3)) / (E(1) - E(3)) ) + &
    !-!                   & ( (E0 - E(4)) / (E(1) - E(4)) ) ) / 3.0_dp

    !-!         do vv = 2, 4
    !-!             Iver(vv) = ( (E0 - E(1)) / (E(vv) - E(1)) ) / 3.0_dp
    !-!         end do
    !-!         ! ------------------- value of g and values of I at vertices ------------------- !

    !-!         enterCase = .true.
    !-!         c_is = g * Iver


    !-!     else if ( ( (E(2) - E0) <= cmp_prec ) .and. ( (E0 - E(3)) <= cmp_prec ) ) then tetra_coef_case

    !-!         ! ------------------- value of g and values of I at vertices ------------------- !
    !-!         g = ( 3.0_dp / ( (E(2) - E(3)) * &
    !-!           &              (E(4) - E(1)) ) ) * &
    !-!           & ( ( ( (E0 - E(1)) * (E0 - E(3)) ) / &
    !-!           &       (E(3) - E(1)) ) - &
    !-!           &   ( ( (E0 - E(2)) * (E0 - E(4)) ) / &
    !-!           &       (E(2) - E(4)) ) )

    !-!         Iver(1) = ( ( (E0 - E(4)) / (E(1) - E(4)) ) / 3.0_dp ) + &
    !-!                 & ( ( ( (E0 - E(3)) / (E(1) - E(3)) ) * &
    !-!                 &     ( (E0 - E(1)) / (E(3) - E(1)) ) * &
    !-!                 &     ( (E0 - E(3)) / (E(2) - E(3)) ) ) / &
    !-!                 &   ( g * (E(4) - E(1)) ) )

    !-!         Iver(2) = ( ( (E0 - E(3)) / (E(2) - E(3)) ) / 3.0_dp ) + &
    !-!                 & ( ( ( ((E0 - E(4)) / (E(2) - E(4))) ** 2 ) * &
    !-!                 &     ( (E0 - E(2)) / (E(3) - E(2)) ) ) / &
    !-!                 &   ( g * (E(4) - E(1)) ) )

    !-!         Iver(3) = ( ( (E0 - E(2)) / (E(3) - E(2)) ) / 3.0_dp ) + &
    !-!                 & ( ( ( ((E0 - E(1)) / (E(3) - E(1))) ** 2 ) * &
    !-!                 &     ( (E0 - E(3)) / (E(2) - E(3)) ) ) / &
    !-!                 &   ( g * (E(4) - E(1)) ) )
    !-!     
    !-!         Iver(4) = ( ( (E0 - E(1)) / (E(4) - E(1)) ) / 3.0_dp ) + &
    !-!                 & ( ( ( (E0 - E(2)) / (E(4) - E(2)) ) * &
    !-!                 &     ( (E0 - E(4)) / (E(2) - E(4)) ) * &
    !-!                 &     ( (E0 - E(2)) / (E(3) - E(2)) ) ) / &
    !-!                 &   ( g * (E(4) - E(1)) ) )
    !-!         ! ------------------- value of g and values of I at vertices ------------------- !

    !-!         enterCase = .true.
    !-!         c_is = g * Iver


    !-!     else if ( ( (E(3) - E0) <= cmp_prec ) .and. ( (E0 - E(4)) <= cmp_prec ) ) then tetra_coef_case

    !-!         ! ------------------- value of g and values of I at vertices ------------------- !
    !-!         g = ( 3.0_dp * ((E0 - E(4)) ** 2) ) / &
    !-!           & ( ( E(4) - E(1) ) * ( E(4) - E(2) ) * &
    !-!           &   ( E(4) - E(3) ) )

    !-!         do vv = 1, 3
    !-!             Iver(vv) = ( (E0 - E(4)) / (E(vv) - E(4)) ) / 3.0_dp
    !-!         end do

    !-!         Iver(4) = ( ( (E0 - E(1)) / (E(4) - E(1)) ) + &
    !-!                 &   ( (E0 - E(2)) / (E(4) - E(2)) ) + &
    !-!                 &   ( (E0 - E(3)) / (E(4) - E(3)) ) ) / 3.0_dp
    !-!         ! ------------------- value of g and values of I at vertices ------------------- !

    !-!         enterCase = .true.
    !-!         c_is = g * Iver

    !-!     else tetra_coef_case

    !-!         does_not_enters: if ( (E0 >= E(1)) .and. &
    !-!                             & (E0 <= E(4)) .and. (.not. enterCase) ) then
    !-!             
    !-!             write(*, 55) E0
    !-!             write(*, 65) E

    !-!             55 FORMAT("WARNING: E0 = ", F12.5, "THz. It does not enters any cases.")
    !-!             65 FORMAT("         The values of E at vertices: [", 3(F12.5, ', '), F12.5, "] THz")

    !-!         end if does_not_enters

    !-!         g = 0.0_dp
    !-!         Iver(:) = 0.0_dp
    !-!         c_is = g * Iver
    !-!         enterCase = .false.

    !-!     end if tetra_coef_case

    !-!     chk_enterCase: if ( enterCase ) then

    !-!         chk_sumeq1: if ( dabs(sum(Iver) - 1.0_dp) > zero_prec ) then

    !-!             write(*, 45) Iver
    !-!             45 FORMAT("WARNING: sum( Iver(:) ) .ne. 1. Values of Iver(1:4): [", 3(F12.5, ', '), F12.5, "]")

    !-!         end if chk_sumeq1

    !-!     end if chk_enterCase

    !-! end Function DiracDeltaCoeffVertices_1979

    ! **** Above Procedures were created just to check the correctness **** !

end module GFIntegrationCoeff

