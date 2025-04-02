
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


!** Why I am doing this in my life ??? **!
program main

    use kinds,                      only : dp
    use parse_cmd_line,             only : get_command_line, ShowWelcomeBanner
    use unit_cell,                  only : cell
    use timer_class,                only : timer
    use FC2_mod,                    only : FC2type
    use FC3_mod,                    only : FC3type
    use EwaldMod,                   only : EwaldParam
    use CreateMesh,                 only : mesh_points, q_points_highsymm
    use DistributeQs,               only : DistributeQpoints, ShowProcessorDistribution
    use phonon_m,                   only : Phon
    use constants,                  only : MAX_ITER

    use Irr_q_point,                only : q_points_data
    use ConserveQuasiMomentum,      only : AllQ
    use IsoScatter,                 only : CalculateIsoW, TauIsoTetra, TetrahedronIterIso
    use Helper,                     only : ReadkappaInfo, CreateFzero, multiplyTaudFnew, kappa

#ifndef GNU
    use ifport,                     only : SYSTEM
#endif

    implicit none

    type(cell)                                              :: sys
    type(timer)                                             :: time_passed
    type(FC2type)                                           :: FC2
    type(FC3type)                                           :: FC3

    type(Phon)                                              :: phonon_dat
    type(EwaldParam)                                        :: EwaldConst

    type(AllQ)                                              :: my_data

    type(q_points_data)                                     :: Qpoints

    real(dp)                                                :: vol
    real(dp)                                                :: T, Temp_k
    real(dp)                                                :: eps_k, k_accu, time_limit

    ! .............:::::::::::::: Coindexed Object ::::::::::::::............. !
    real(dp), allocatable, dimension(:,:), codimension[:]   :: Inv_tauqs
    real(dp), allocatable, dimension(:,:,:), codimension[:] :: FZero
    real(dp), allocatable, dimension(:,:,:), codimension[:] :: dF_new

    real(dp), dimension(9), codimension[*]                  :: k_tensor
    ! .............:::::::::::::: Coindexed Object ::::::::::::::............. !

    real(dp), dimension(9)                                  :: k_tensor_old
    real(dp), allocatable, dimension(:,:,:)                 :: F_vecAll

    real(dp), dimension(:,:,:,:), allocatable               :: W_iso
    real(dp), dimension(:), allocatable                     :: g2iso

    character(len=256)                                      :: filename, RestartDir, msg
    character(len=512)                                      :: command
    character(len=24)                                       :: Temp_char
    character(len=8)                                        :: frmt
    character(len=64)                                       :: fc2file, fc3file

    integer, dimension(3)                                   :: mesh, shift
    integer                                                 :: time_rev

    integer                                                 :: Nqpoints, num_irr_q, Nbasis, Ndof
    integer                                                 :: my_Qsize, my_offset
    integer, dimension(2)                                   :: my_edge, progress
    integer                                                 :: iterNum, err

    logical                                                 :: isotope, LongEW, ReadRestart, &
                                                             & STOP_FLAG

    if ( this_image() == 1 ) call ShowWelcomeBanner()

    call time_passed%start_timer()

    call get_command_line(filename, T, Temp_k)
    call sys%init_data(filename)

    Nbasis = sys%natm
    vol = sys%vol
    Ndof = 3 * Nbasis

    frmt = '(F6.1)'
    write(Temp_char, frmt) T

    fc2file = 'FC2nd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
    call FC2%set_FC2(fc2file)

    fc3file = 'FC3rd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
    call FC3%set_FC3(sys, fc3file)

    allocate( g2iso(Nbasis) )
    g2iso = 0.0_dp

    call ReadkappaInfo(Nbasis, filename, mesh, isotope, g2iso, eps_k, &
                     & time_limit, LongEW, ReadRestart)

    Nqpoints = product( mesh )

    !* Restart *!
    write(Temp_char, frmt) Temp_k
    RestartDir = 'Restart_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K'
    command = 'mkdir -p -v ./'//trim(adjustl(adjustr(RestartDir)))

    CreateDir: if ( this_image() == 1 ) then

#ifdef GNU
        call SYSTEM(command, status=err)
#else
        err = SYSTEM(command)
#endif
        dir_create_err: if ( err /= 0 ) then
            write(*, *) 'Directory creating failed: iostat = ', err
        end if dir_create_err

    end if CreateDir
    !* Restart *!

    shift(:) = 0
    time_rev = 0

    call Qpoints%make_tetrahedron( sys, mesh, shift )
    num_irr_q = size( Qpoints%irr_q_int, 2 )

    ! ======================= find phonon freq., velocity and eigenvectors ====================== !

    Ew_chk: if ( LongEW ) then

        call EwaldConst%set_EwaldParam(sys, filename)

        img1_chk1: if ( this_image() == 1 ) then
            write(*, *)
            call execute_command_line(' ')
            write(*, 100) EwaldConst%Lmb
            call execute_command_line(' ')

            write(*, 200) EwaldConst%Gmesh
            call execute_command_line(' ')
            write(*, *)
            call execute_command_line(' ')

            100 FORMAT('The alpha parameter for Ewald Sum: ', ES10.3)
            200 FORMAT('The NKcut for Ewald Sum:', 3I3)
        end if img1_chk1

        SYNC ALL

    end if Ew_chk

    call phonon_dat%set_Phon_mp(sys, FC2, EwaldConst, Qpoints%q_pnt, Temp_k, Ndof, LongEW)

    img1_chk: if ( this_image() == 1 ) then
        write(*, *)
        call execute_command_line(' ')
        write(*, 12) num_irr_q
        call execute_command_line(' ')
        write(*, *)
        call execute_command_line(' ')

        12 FORMAT('Number of irreducible q-points: ', I5)
    end if img1_chk
    SYNC ALL

    ! ======================= find phonon freq., velocity and eigenvectors ====================== !

    ! ========================== distribute the q's through the images ========================== !
    ! For each images there are my_Qsize number of qs. So the shape of the allocatable array is   !
    ! my_Qsize with index from 1 to my_Qsize. To recover the q-point add my_offset with the index !

    !*  New Distribute  *!
    call DistributeQpoints(num_irr_q, my_Qsize, my_offset, my_edge)
    !*  New Distribute  *!

    call ShowProcessorDistribution( "q-points distribution", &
                                  & my_edge, my_Qsize, my_offset )
    SYNC ALL

    ! ========================== distribute the q's through the images ========================== !

    ! =============================== coindexed object allocation =============================== !

    allocate( Inv_tauqs(Ndof, num_irr_q)[*], stat=err, ERRMSG=msg )
    if ( err /= 0 ) write(*, "(' Memory allocation ERROR in image : ', I5, '. Error message: ', A128)") this_image(), msg
    SYNC ALL
    Inv_tauqs = 0.0_dp

    allocate( FZero(3, Ndof, num_irr_q)[*], stat=err, ERRMSG=msg )
    if ( err /= 0 ) write(*, "(' Memory allocation ERROR in image : ', I5, '. Error message: ', A128)") this_image(), msg
    SYNC ALL
    FZero = 0.0_dp

    allocate( dF_new(3, Ndof, num_irr_q)[*], stat=err, ERRMSG=msg )
    if ( err /= 0 ) write(*, "(' Memory allocation ERROR in image : ', I5, '. Error message: ', A128)") this_image(), msg
    SYNC ALL
    dF_new = 0.0_dp

    ! =============================== coindexed object allocation =============================== !

    allocate( F_vecAll(3, Ndof, num_irr_q), stat=err, ERRMSG=msg )
    if ( err /= 0 ) write(*, "(' Memory allocation ERROR in image : ', I5, '. Error message: ', A128)") this_image(), msg

    allocate( W_iso(Ndof, Nqpoints, Ndof, my_Qsize), stat=err, ERRMSG=msg )
    if ( err /= 0 ) write(*, "(' Memory allocation ERROR in image : ', I5, '. Error message: ', A128)") this_image(), msg
    W_iso = 0.0_dp

    call my_data%QconsrvMomentum(time_passed, time_limit, Qpoints, mesh, my_Qsize, my_offset, &
                               & sys, FC3, phonon_dat, RestartDir, ReadRestart, progress, STOP_FLAG)

    STOP_EXEC: if ( STOP_FLAG ) then

        SYNC ALL
        ERROR STOP

    end if STOP_EXEC

    if ( isotope ) then
        call CalculateIsoW( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, Nbasis, &
                          & Ndof, Nqpoints, g2iso, W_iso )
    end if 

    call my_data%TetrahedronWplus( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                                 & Ndof, num_irr_q, Inv_tauqs )

    call my_data%TetrahedronWminus( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                                  & Ndof, num_irr_q, Inv_tauqs )


    if ( isotope ) then
        call TauIsoTetra( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, Ndof, &
                        & Nqpoints, num_irr_q, W_iso, Inv_tauqs )
    end if

    SYNC ALL
    call co_sum( Inv_tauqs, stat=err, ERRMSG=msg )
    if ( err /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
    SYNC ALL

    call CreateFzero( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, Ndof, num_irr_q, &
                    & Temp_k, Inv_tauqs, FZero)

    SYNC ALL
    call co_sum( FZero, stat=err, ERRMSG=msg )
    if ( err /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
    SYNC ALL

    F_vecAll = FZero

    k_tensor = kappa( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                    & Ndof, num_irr_q, Nqpoints, F_vecAll, vol )

    SYNC ALL
    call co_sum( k_tensor, stat=err, ERRMSG=msg )
    if ( err /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
    SYNC ALL

    iterNum = 0
    img1_chk3: if ( this_image() == 1 ) then

        write(*, *)
        write(*, "(5X, '|', 40('-'), 'Iterative solution of BTE', 39('-'), '|')")
        write(*, "(5X, '|', 104X, '|')")

        write(*, 1001) k_tensor(1:3)
        write(*, 1002) iterNum, k_tensor(4:6)
        write(*, 1003) k_tensor(7:9)
        write(*, 1004)

    end if img1_chk3

    k_tensor_old = k_tensor

    sc_iter: do

        iterNum = iterNum + 1

        call my_data%TetrahedronIterPlus( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                                        & Ndof, num_irr_q, F_vecAll, dF_new)

        call my_data%TetrahedronIterMinus( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                                         & Ndof, num_irr_q, F_vecAll, dF_new )

        if ( isotope ) then

            call TetrahedronIterIso( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                                   & Ndof, Nqpoints, num_irr_q, W_iso, F_vecAll, dF_new )

        end if

        call multiplyTaudFnew( Qpoints, my_Qsize, my_offset, Ndof, num_irr_q, Inv_tauqs, dF_new )

        SYNC ALL
        call co_sum( dF_new, stat=err, ERRMSG=msg )
        if ( err /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
        SYNC ALL

        F_vecAll = FZero + dF_new ! This is new updated F vector

        k_tensor = kappa( Qpoints, phonon_dat, mesh, my_Qsize, my_offset, &
                        & Ndof, num_irr_q, Nqpoints, F_vecAll, vol )

        SYNC ALL
        call co_sum( k_tensor, stat=err, ERRMSG=msg )
        if ( err /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
        SYNC ALL

        img1_chk4: if ( this_image() == 1 ) then

            write(*, 1001) k_tensor(1:3)
            write(*, 1002) iterNum, k_tensor(4:6)
            write(*, 1003) k_tensor(7:9)
            write(*, 1004)

        end if img1_chk4

        k_accu = maxval( dabs(k_tensor - k_tensor_old) )

        ExitLoop: if ( k_accu < eps_k ) then

            if ( this_image() == 1 ) then
                write(*, 222) iterNum
                write(*, "(5X, '|', 104('-'), '|')")
            end if
            exit sc_iter

        else if ( iterNum > MAX_ITER ) then ExitLoop

            if ( this_image() == 1 ) then 
                write(*, 777) iterNum
                write(*, "(5X, '|', 104('-'), '|')")
            end if
            exit sc_iter

        end if ExitLoop

        k_tensor_old = k_tensor
        dF_new = 0.0_dp

    end do sc_iter

    1001 FORMAT( 5X, '|', 21X, '|kxx = ', G20.12, ', kxy = ', G20.12, ', kxz = ', G20.12, '|' )
    1002 FORMAT( 5X, '| iteration no: ', I5, 1X, '|kyx = ', G20.12, ', kyy = ', G20.12, ', kyz = ', G20.12, '|' )
    1003 FORMAT( 5X, '|', 21X, '|kzx = ', G20.12, ', kzy = ', G20.12, ', kzz = ', G20.12, '|' )
    1004 FORMAT( 5X, '|', 104('-'), '|')

    222 FORMAT(5X, '|**', 19X, 'Convergence in Self Consistence loop achieved at iter no: ', I4, 19X, '**|')
    777 FORMAT(5X, '|', 26X, 'Convergence is not achieved, MAX_ITER exceeded: ', I4, 26X, '|')

    SYNC ALL
    if ( this_image() == 1) write(*, '(/)')

end program main

