
module TetrahedronIntegration

    use kinds,          only : dp
    use constants,      only : zero_prec, OnebySix
    use phonon_m,       only : Phon

    implicit none

    EXTERNAL        :: dlasrt

    private
    public          :: FindFandOmegaAtVertices, FindgAndIver, FindIverxFver

contains

    subroutine my_sort( Omega_vert, F_ver, tetra_ver_q1, tetra_ver_q2 )

        implicit none

        real(dp), dimension(4), intent(inout)                       :: Omega_vert
        real(dp), dimension(4), intent(inout)                       :: F_ver
        integer, dimension(4), intent(inout)                        :: tetra_ver_q1
        integer, dimension(4), intent(inout)                        :: tetra_ver_q2

        !! ===================================== Local variables =================================== !!

        real(dp), dimension(4)                      :: Omega_vertc, F_verc
        integer, dimension(4)                       :: tetra_ver_q1c, tetra_ver_q2c

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
        F_verc = F_ver

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
            F_ver(ii) = F_verc( ActualIndx )

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


    subroutine FindgAndIver( Mqs, g, Iver, Mq1q2s, F_ver_srt, ver_q1_srt, ver_q2_srt, enterCase )

        implicit none

        real(dp), intent(in)                                        :: Mqs
        real(dp), intent(out)                                       :: g
        real(dp), dimension(4), intent(out)                         :: Iver
        real(dp), dimension(4), intent(inout)                       :: Mq1q2s, F_ver_srt
        integer, dimension(4), intent(inout)                        :: ver_q1_srt, ver_q2_srt
        logical, intent(out)                                        :: enterCase
        
        ! ====================================== Local Variables ====================================== !

        integer                                 :: vv

        ! ====================================== Local Variables ====================================== !

        call my_sort(Mq1q2s, F_ver_srt, ver_q1_srt, ver_q2_srt)

        enterCase = .false.

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
                    &     ( (Mqs - Mq1q2s(1)) / (Mq1q2s(3) - Mq1q2s(1)) ) * &
                    &     ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) ) / &
                    &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )

            Iver(2) = ( ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) / 3.0_dp ) + &
                    & ( ( ( ((Mqs - Mq1q2s(4)) / (Mq1q2s(2) - Mq1q2s(4))) ** 2 ) * &
                    &     ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) ) / &
                    &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )

            Iver(3) = ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) / 3.0_dp ) + &
                    & ( ( ( ((Mqs - Mq1q2s(1)) / (Mq1q2s(3) - Mq1q2s(1))) ** 2 ) * &
                    &     ( (Mqs - Mq1q2s(3)) / (Mq1q2s(2) - Mq1q2s(3)) ) ) / &
                    &   ( g * (Mq1q2s(4) - Mq1q2s(1)) ) )
        
            Iver(4) = ( ( (Mqs - Mq1q2s(1)) / (Mq1q2s(4) - Mq1q2s(1)) ) / 3.0_dp ) + &
                    & ( ( ( (Mqs - Mq1q2s(2)) / (Mq1q2s(4) - Mq1q2s(2)) ) * &
                    &     ( (Mqs - Mq1q2s(4)) / (Mq1q2s(2) - Mq1q2s(4)) ) * &
                    &     ( (Mqs - Mq1q2s(2)) / (Mq1q2s(3) - Mq1q2s(2)) ) ) / &
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

            g = 0.0_dp
            Iver(:) = 0.0_dp
            enterCase = .false.

        end if tetra_coef_case

        chk_enterCase: if ( enterCase ) then

            chk_sumeq1: if ( dabs(sum(Iver) - 1.0_dp) > zero_prec ) then

                write(*, 45) Iver
                45 FORMAT("WARNING: sum( Iver(:) ) .ne. 1. Values of Iver(1:4): [", 3(F12.5, ', '), F12.5, "]")

            end if chk_sumeq1

        end if chk_enterCase

    end subroutine FindgAndIver

    !** Debug Purpose **!
    subroutine FindIverxFver( Iver, F_ver_srt, enterCase, IxF )

        implicit none

        real(dp), dimension(4), intent(in)                  :: Iver
        real(dp), dimension(4), intent(in)                  :: F_ver_srt
        logical, intent(in)                                 :: enterCase

        real(dp), intent(out)                               :: IxF

        ! =================================== Local Variables =================================== !

        integer                                             :: vv

        ! =================================== Local Variables =================================== !

        chk_entercase: if ( enterCase ) then

            IxF = 0.0_dp

            vert_loop: do vv = 1, 4

                IxF = IxF + Iver(vv) * F_ver_srt(vv)

            end do vert_loop

        else chk_entercase

            IxF = 0.0_dp

        end if chk_entercase

        IxF = OnebySix * IxF

    end subroutine FindIverxFver
    !** Debug Purpose **!

    subroutine FindFandOmegaAtVertices( Ndof, Nq, s0, s1, s2, &
                                      & q1IndxVer_unsrt, ph_q1, W_s0s1s2, &
                                      & omega_q2All, nBE_q2All, &
                                      & F_ver1, F_ver2, M_ver1, M_ver2, M_ver3)

        implicit none

        integer, intent(in)                                     :: Ndof, Nq, s0, s1, s2
        integer, dimension(4), intent(in)                       :: q1IndxVer_unsrt
        type(Phon), intent(in)                                  :: ph_q1
        real(dp), dimension(Ndof,Ndof,Ndof,Nq), intent(in)      :: W_s0s1s2
        real(dp), dimension(Ndof, Nq), intent(in)               :: omega_q2ALl, nBE_q2All

        real(dp), dimension(4), intent(out)                     :: F_ver1, F_ver2, & 
                                                                 & M_ver1, M_ver2, M_ver3

        ! ==================================== Local Variables ==================================== !

        integer                             :: vv

        ! ==================================== Local Variables ==================================== !

        vert_loop: do vv = 1, 4

            F_ver1(vv) = W_s0s1s2( s2, s1, s0, q1IndxVer_unsrt(vv) ) * &
                       & ( ph_q1%nBE(s1, q1IndxVer_unsrt(vv)) + nBE_q2ALl(s2, q1IndxVer_unsrt(vv)) + 1.0_dp )

            F_ver2(vv) = W_s0s1s2( s2, s1, s0, q1IndxVer_unsrt(vv) ) * &
                       & ( ph_q1%nBE(s1, q1IndxVer_unsrt(vv)) - nBE_q2All(s2, q1IndxVer_unsrt(vv)) )

            M_ver1(vv) = ph_q1%omega(s1, q1IndxVer_unsrt(vv)) + omega_q2All(s2, q1IndxVer_unsrt(vv))
            M_ver2(vv) = omega_q2All(s2, q1IndxVer_unsrt(vv)) - ph_q1%omega(s1, q1IndxVer_unsrt(vv))
            M_ver3(vv) = ph_q1%omega(s1, q1IndxVer_unsrt(vv)) - omega_q2All(s2, q1IndxVer_unsrt(vv))

        end do vert_loop

    end subroutine FindFandOmegaAtVertices

end module TetrahedronIntegration

