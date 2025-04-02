
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

module Helper

    use kinds,              only : dp
    use constants,          only : FunitE_11m_K, zero_prec, kW_mK
    use Irr_q_point,        only : q_points_data
    use phonon_m,           only : Phon

    implicit none
    private

    public                  :: ReadkappaInfo, CreateFzero, multiplyTaudFnew, kappa

contains

    subroutine ReadkappaInfo(Nbasis, filename, mesh, isotope, g2iso, eps_k, &
                           & time_limit, LongEW, ReadRestart) 

        implicit none

        integer, intent(in)                         :: Nbasis
        character(len=*), intent(in)                :: filename

        integer, dimension(3), intent(out)          :: mesh
        logical, intent(out)                        :: isotope
        real(dp), dimension(Nbasis), intent(out)    :: g2iso
        real(dp), intent(out)                       :: eps_k, time_limit
        logical, intent(out)                        :: LongEw, ReadRestart

        ! ================================= Local variable =================================== !

        integer, dimension(3)                       :: Recpmesh
        logical                                     :: IsoScatter
        real(dp), dimension(Nbasis)                 :: massVar
        logical                                     :: Ew_correc, Restart
        real(dp)                                    :: kappa_accu, Maxtime

        namelist        /KappaInfo/         Recpmesh, IsoScatter, massVar, Ew_correc, Restart, kappa_accu, Maxtime

        integer                                     :: err
        character(len=512)                          :: err_msg

        ! ================================= Local variable =================================== !

        Maxtime = 5.0_dp * 24.0_dp * 60.0_dp !minutes
        Restart = .false.

        open(unit=5, file=filename, status='OLD', iostat=err, iomsg=err_msg, &
           & action='READ', delim='APOSTROPHE')

            open_chk1: if ( err /= 0 ) then
                write(*, *) 'Input file OPEN failed: iostat = ', err
                write(*, *) 'Error message : ', err_msg
            end if open_chk1

            read(unit=5, nml=KappaInfo, iostat=err, iomsg=err_msg)

            ReadChk: if ( err /= 0 ) then
                write(*, 134) filename
                write(*, *) 'Error message : ', err_msg

                134 FORMAT('Reading Namelist file: ', A, ' failed')
            end if ReadChk

        close(unit=5, status='KEEP', iostat=err)

        mesh = Recpmesh
        isotope = IsoScatter
        g2iso = massVar
        LongEw = Ew_correc
        eps_k = kappa_accu

        time_limit = Maxtime
        ReadRestart = Restart

        img1_chk: if ( this_image() == 1 ) then

            write(*, *)
            write(*, 10) mesh
            write(*, 20) isotope
            !if ( isotope) write(*, 30) g2iso
            write(*, 40) LongEW
            write(*, 50) eps_k
            write(*, 60) time_limit
            write(*, 70) ReadRestart

            10 FORMAT("q-mesh for kappa is set to : ", 2(I4, ' x '), I4)
            20 FORMAT("Calculate isotopic scattering: ", L4)
            40 FORMAT("Ewald Correction in phonon calcualtion: ", L4)
            50 FORMAT("Accuracy for iterative solution of kappa: ", E14.6)
            60 FORMAT("Maximum execution time before writing restart file: ", F7.2, " min.")
            70 FORMAT("Starting from restart file ? : ", L4)

        end if img1_chk

    end subroutine ReadkappaInfo


    subroutine CreateFzero( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, Ndof, num_irr_q, &
                          & Temp, Inv_tauqs, FZero)

        implicit none

        type(q_points_data), intent(in)                             :: Qpoints
        type(Phon), intent(in)                                      :: phonon_dat

        integer, dimension(3), intent(in)                           :: mesh
        integer, intent(in)                                         :: my_Qsize, my_offset, &
                                                                     & Ndof, num_irr_q

        real(dp), intent(in)                                        :: Temp
        real(dp), dimension(Ndof, num_irr_q), intent(in)            :: Inv_tauqs

        real(dp), dimension(3, Ndof, num_irr_q), intent(inout)      :: FZero

        ! ================================= Local variable =================================== !
        
        real(dp)                                    :: tauqs_inv

        integer                                     :: qindx, qindx2, strt
        integer, dimension(3)                       :: q_int, mesh_mul
        integer                                     :: qi, s

        logical                                     :: q0chk

        ! ================================= Local variable =================================== !

        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)

        irrQloop: do qi = 1, my_Qsize

            qindx = my_offset + qi
            q_int = Qpoints%irr_q_int(1:3, qindx)
            qindx2 = Qpoints%indx_map( dot_product( modulo( q_int, mesh ), mesh_mul ) + 1 )

            q0chk = all( q_int == 0 )

            avoid_0div: if (q0chk ) then
                strt = 4
            else avoid_0div
                strt = 1
            end if avoid_0div

            s_loop: do s = strt, Ndof

                tauqs_inv = Inv_tauqs(s, qindx)

                tau0chk: if ( dabs(tauqs_inv) > zero_prec ) then

                    FZero(:, s, qindx) = FunitE_11m_K * ( phonon_dat%omega(s, qindx2) ) * &
                                       & ( phonon_dat%grp_vel(:, s, qindx2) ) / &
                                       & ( tauqs_inv * (Temp ** 2) )

                else tau0chk

                    write(*, 100) qindx, s, tauqs_inv
                    100 FORMAT("WARNING(in CreateFzero): tau inverse is zero, for qindx = ", I4, &
                             & ", s = ", I4, ". tau^-1 = ", E12.5)

                end if tau0chk

            end do s_loop

        end do irrQloop

    end subroutine CreateFzero


    subroutine multiplyTaudFnew( Qpoints, my_Qsize, my_offset, Ndof, num_irr_q, Inv_tauqs, dFnew )

        implicit none

        type(q_points_data), intent(in)                             :: Qpoints

        integer, intent(in)                                         :: my_Qsize, my_offset, &
                                                                     & Ndof, num_irr_q

        real(dp), dimension(Ndof, num_irr_q), intent(in)            :: Inv_tauqs !Explicit-shape dummy Array

        real(dp), dimension(3, Ndof, num_irr_q), intent(inout)      :: dFnew !Explicit-shape dummy Array

        ! ================================= Local variable =================================== !
        
        real(dp)                                    :: tauqs_inv

        integer                                     :: qindx, strt
        integer, dimension(3)                       :: q_int
        integer                                     :: qi, s

        logical                                     :: q0chk

        ! ================================= Local variable =================================== !

        irrQloop: do qi = 1, my_Qsize

            qindx = my_offset + qi
            q_int = Qpoints%irr_q_int(1:3, qindx)

            q0chk = all( q_int == 0 )

            avoid_0div: if (q0chk ) then
                strt = 4
            else avoid_0div
                strt = 1
            end if avoid_0div

            s_loop: do s = strt, Ndof

                tauqs_inv = Inv_tauqs(s, qindx)

                tau0chk: if ( dabs(tauqs_inv) > zero_prec ) then

                    dFnew(:, s, qindx) = dFnew(:, s, qindx) / tauqs_inv

                else tau0chk

                    write(*, 100) qindx, s, tauqs_inv
                    100 FORMAT("WARNING(in multiplyTau_dF): tau inverse is zero, for qindx = ", I4, &
                             & ", s = ", I4, ". tau^-1 = ", E12.5)

                end if tau0chk

            end do s_loop

        end do irrQloop

    end subroutine multiplyTaudFnew


    Function kron_prod(invec1, invec2) Result(outvec)

        implicit none

        real(dp), dimension(3), intent(in)                  :: invec1, invec2
        real(dp), dimension(9)                              :: outvec !Result

        ! ================================== Local Variables ================================== !

        integer                                             :: vi, strtp, endp

        ! ================================== Local Variables ================================== !

        strtp = 1
        endp = 3

        invec1_loop: do vi = 1, 3

            outvec(strtp : endp) = invec1(vi) * invec2(:)

            strtp = endp + 1
            endp = endp + 3

        end do invec1_loop

    end Function kron_prod


    Function kappa( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                  & Ndof, num_irr_q, Nqpoints, F_new, vol ) Result(k_tensor)

        implicit none

        type(q_points_data), intent(in)                             :: Qpoints
        type(Phon), intent(in)                                      :: phonon_dat

        integer, dimension(3), intent(in)                           :: mesh
        integer, intent(in)                                         :: my_Qsize, my_offset, Ndof, &
                                                                     & num_irr_q, Nqpoints

        real(dp), dimension(3, Ndof, num_irr_q), intent(in)         :: F_new !Explicit-shape dummy Array
        real(dp), intent(in)                                        :: vol

        real(dp), dimension(9)                                      :: k_tensor !Result

        ! ================================= Local variable =================================== !

        real(dp), dimension(3)                      :: grp_vel, F_vec
        real(dp), dimension(9)                      :: k_tns
        real(dp)                                    :: Mqs, nqs, multiply_fac
        
        integer                                     :: qindx, qPosindx, strt, deg
        integer, dimension(3)                       :: q_int, mesh_mul
        integer                                     :: qi, s
        logical                                     :: q0chk

        ! ================================= Local variable =================================== !

        multiply_fac = 1.0_dp / (dble(Nqpoints) * vol)

        mesh_mul = (/1, mesh(1), mesh(1)*mesh(2)/)

        k_tensor = 0.0_dp

        irrQloop: do qi = 1, my_Qsize

            qindx = my_offset + qi
            q_int = Qpoints%irr_q_int(1:3, qindx)
            deg = Qpoints%irr_q_int(4, qindx)

            qPosindx = Qpoints%indx_map( dot_product( modulo( q_int, mesh ), mesh_mul ) + 1 )

            q0chk = all( q_int == 0 )

            avoid_0div: if (q0chk ) then
                strt = 4
            else avoid_0div
                strt = 1
            end if avoid_0div

            s_loop: do s = strt, Ndof

                Mqs = phonon_dat%omega(s, qPosindx)
                nqs = phonon_dat%nBE(s, qPosindx)
                grp_vel = phonon_dat%grp_vel(:, s, qPosindx)

                F_vec = F_new(:, s, qindx)

                k_tns = kron_prod(grp_vel, F_vec)

                k_tensor = k_tensor + ( (dble(deg) * Mqs * nqs * (nqs + 1.0_dp)) * k_tns )

            end do s_loop

        end do irrQloop
        
        k_tensor = ( kW_mK * multiply_fac * k_tensor )

    end Function kappa

end module Helper

