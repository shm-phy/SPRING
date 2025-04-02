
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

program main

    use kinds,                      only : dp
    use parse_cmd_line,             only : get_command_line, ShowWelcomeBanner
    use unit_cell,                  only : cell
    use timer_class,                only : timer
    use FC2_mod,                    only : FC2type
    use FC3_mod,                    only : FC3type
    use EwaldMod,                   only : EwaldParam
    use CreateMesh,                 only : q_points_highsymm
    use TetrahedronQ,               only : q_tetrahedron
    use Helper,                     only : DistributeQpoints, ShowProcessorDistribution, &
                                         & WriteAsciiOutFile
    use phonon_m,                   only : Phon

    use Linewidth,                  only : ScatterProbAllq1, CalculateLinewidth, &
                                         & CalculateLwTetrahedron

    implicit none

    type(cell)                                              :: sys
    type(timer)                                             :: time_passed
    type(FC2type)                                           :: FC2
    type(FC3type)                                           :: FC3

    type(q_tetrahedron)                                     :: tetraQZeroCentered, &
                                                             & tetraQNotZeroCentered

    type(EwaldParam)                                        :: EwaldConst
    type(Phon)                                              :: ph_q0, ph_ZeroCntrdq1, &
                                                             & ph_NotZeroCntrdq1

    character(len=256)                                      :: filename
    character(len=24)                                       :: Temp_char, Ndof_char
    character(len=8)                                        :: frmt
    character(len=64)                                       :: fc2file, fc3file
    character(len=64)                                       :: outfile_ascii
    character(len=256)                                      :: err_msg
    character(len=128)                                      :: outfile_frmt

    real(dp)                                                :: vol
    real(dp)                                                :: T, freq_cut
    real(dp)                                                :: time_2nd, time_LW
    real(dp), dimension(3)                                  :: q0_single

    real(dp), dimension(:, :), allocatable                  :: q0
    real(dp), dimension(:), allocatable                     :: dist_arr_q0, lw

    real(dp), allocatable                                   :: q_high_sym(:, :)
    namelist        /HighSymPathLW/     q_high_sym

    real(dp)                                                :: Sigma, OnebySigma, &
                                                             & OnebySigma2
    integer, dimension(3)                                   :: qmesh
    integer                                                 :: num_points
    namelist        /LWInfo/     Sigma, qmesh, num_points

    ! ..............:::::::::::::: Coindexed Object ::::::::::::::.............. !
    real(dp), allocatable, dimension(:,:,:,:), codimension[:]       :: Sctrq1
    real(dp), allocatable, dimension(:), codimension[:]             :: LW_q0

    real(dp), dimension(:, :), allocatable, codimension[:]          :: omegaq2
    real(dp), dimension(:, :), allocatable, codimension[:]          :: nBEq2
    real(dp), dimension(:, :, :), allocatable, codimension[:]       :: grp_velq2
    integer, dimension(:), allocatable, codimension[:]              :: Ifq2Exists
    ! ..............:::::::::::::: Coindexed Object ::::::::::::::.............. !

    integer                                                 :: Nbasis, Ndof, Nq1points, &
                                                             & Nq0points, Ntetrahedrons

    integer, dimension(3)                                   :: q1mesh, shift

    integer                                                 :: my_Qsize, my_offset
    integer, dimension(2)                                   :: my_edge

    integer                                                 :: my_QsizeT, my_offsetT
    integer, dimension(2)                                   :: my_edgeT

    integer                                                 :: i0

    integer                                                 :: err

    logical                                                 :: LongEW, ZroCntrd, &
                                                             & HeaderWrite
    logical, dimension(:), allocatable                      :: ZeroCenteredMesh


    if ( this_image() == 1 ) then
        call ShowWelcomeBanner()
        call time_passed%start_timer()
    end if

    call get_command_line(filename, T, freq_cut)
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

    outfile_ascii = 'linewidth_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.txt'

    frmt = '(I4)'
    write(Ndof_char, frmt) (2*Ndof-1)
    outfile_frmt = "(G20.12,'   ', "//trim(adjustl(adjustr(Ndof_char)))//"(G20.12, '  '), G20.12)"

    ! ========================== Read the LW input from namelist file =========================== !

    open(unit=5, file=filename, status='OLD', iostat=err, iomsg=err_msg, &
         action='READ', delim='APOSTROPHE')

        open_chk: if ( err /= 0 ) then
            write(*, 300) err
            300 FORMAT('Input file OPEN failed: iostat =', I3)
            write(*, 400) err_msg
            400 FORMAT('Error message = ', A)
        end if open_chk

        read(unit=5, nml=LWInfo, iostat=err, iomsg=err_msg)        ! ** !

        if ( (err /= 0) .and. (this_image() == 1) ) then
            write(*, 301) filename
            write(*, *) err_msg
            ERROR STOP
            301 FORMAT(" ERROR in reading LWInfo namelist in file: ", A128, "ERROR MESSAGE ==> ")
        end if

        allocate(q_high_sym(4, num_points))

        read(unit=5, nml=HighSymPathLW, iostat=err, iomsg=err_msg) ! ** !

        if ( (err /= 0) .and. (this_image() == 1) ) then
            write(*, 901) filename
            write(*, *) err_msg
            ERROR STOP
            901 FORMAT(" ERROR in reading HighSymPathLW namelist in file: ", A128, "ERROR MESSAGE ==> ")
        end if

    close(unit=5, status='KEEP', iostat=err, iomsg=err_msg)

    close_chk: if ( (err /= 0) ) then
        write(*, 102) filename
        write(*, *) err_msg
        ERROR STOP
        102 FORMAT(" ERROR in closing namelist file: ", A128, "ERROR MESSAGE ==> ")
    end if close_chk

    OnebySigma = 1.0_dp / Sigma
    OnebySigma2 = (OnebySigma ** 2)

    ! ========================== Read the LW input from namelist file =========================== !

    ! =============================== q0 from high symmetry path ================================ !

    call q_points_highsymm(q_high_sym, sys%G, q0, dist_arr_q0, ZeroCenteredMesh)

    Nq0points = size( q0, 2 )

    ! =============================== q0 from high symmetry path ================================ !

    ! ===================== q1 with 0-centered mesh, and not 0-centerd mesh ===================== !

    q1mesh = qmesh !(/5, 5, 5/)
    shift = (/0, 0, 0/)

    Nq1points = product( q1mesh )

    call tetraQZeroCentered%make_tetrahedron( sys, q1mesh, shift, .true. )

    call tetraQNotZeroCentered%make_tetrahedron( sys, q1mesh, shift, .false. )

    ! ===================== q1 with 0-centered mesh, and not 0-centerd mesh ===================== !

    ! ======================= find phonon freq., velocity and eigenvectors ====================== !

    LongEW = sys%Force_dd

    Ew_chk: if ( LongEW ) then

        call EwaldConst%set_EwaldParam(sys, filename)

        img1_chk0: if ( this_image() == 1 ) then
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
        end if img1_chk0

        SYNC ALL

    end if Ew_chk

    call ph_q0%set_Phon_mp( sys, FC2, EwaldConst, q0, T, freq_cut, Ndof, LongEW)

    call ph_ZeroCntrdq1%set_Phon_mp( sys, FC2, EwaldConst, tetraQZeroCentered%q_pnt, &
                                   & T, freq_cut, Ndof, LongEW )

    call ph_NotZeroCntrdq1%set_Phon_mp( sys, FC2, EwaldConst, tetraQNotZeroCentered%q_pnt, &
                                   & T, freq_cut, Ndof, LongEW )

    SYNC ALL

    ! ======================= find phonon freq., velocity and eigenvectors ====================== !

    ! ========================== distribute the q1's through the images ========================= !
    ! For each images there are my_Qsize number of qs. So the shape of the allocatable array is   !
    ! my_Qsize with index from 1 to my_Qsize. To recover the q-point add my_offset with the index !

    !*  New Distribute  *!
    call DistributeQpoints(Nq1points, my_Qsize, my_offset, my_edge)
    !*  New Distribute  *!

    call ShowProcessorDistribution( "q-points distribution", my_edge, my_Qsize, my_offset )

    ! ========================== distribute the q1's through the images ========================= !

    Ntetrahedrons = Nq1points * 6
    ! ====================== distribute the tetrahedrons through the images ===================== !
    ! For each images there are my_Qsize number of qs. So the shape of the allocatable array is   !
    ! my_Qsize with index from 1 to my_Qsize. To recover the q-point add my_offset with the index !

    !*  New Distribute  *!
    call DistributeQpoints(Ntetrahedrons, my_QsizeT, my_offsetT, my_edgeT)
    !*  New Distribute  *!

    call ShowProcessorDistribution( "q-points distribution (for tetrahedron integration)", &
                                  & my_edgeT, my_QsizeT, my_offsetT )

    ! ====================== distribute the tetrahedrons through the images ===================== !

    ! ...............:::::::::::::::: Coindexed Object allocation ::::::::::::::::............... !
    SYNC ALL
    allocate( LW_q0(Ndof)[*] )

    allocate( Sctrq1(Ndof, Ndof, Ndof, Nq1points)[*], STAT=err, ERRMSG=err_msg )
    if ( err /= 0 ) write(*, "(' Memory allocation ERROR in image : ', I5, '. Error message: ', A128)") this_image(), err_msg

    allocate( omegaq2(Ndof, Nq1points)[*], STAT=err, ERRMSG=err_msg )
    if ( err /= 0 ) write(*, "(' Memory allocation ERROR in image : ', I5, '. Error message: ', A128)") this_image(), err_msg

    allocate( nBEq2(Ndof, Nq1points)[*], STAT=err, ERRMSG=err_msg )
    if ( err /= 0 ) write(*, "(' Memory allocation ERROR in image : ', I5, '. Error message: ', A128)") this_image(), err_msg

    allocate( grp_velq2(3, Ndof, Nq1points)[*], STAT=err, ERRMSG=err_msg )
    if ( err /= 0 ) write(*, "(' Memory allocation ERROR in image : ', I5, '. Error message: ', A128)") this_image(), err_msg

    allocate( Ifq2Exists(Nq1points)[*] )
    SYNC ALL
    ! ...............:::::::::::::::: Coindexed Object allocation ::::::::::::::::............... !

    LW_q0 = 0.0_dp
    ! ===================================== Phonon Linewidth ==================================== !

    allocate( lw(Ndof) ) 

    HeaderWrite = .true.
    q0Loop: do i0 = 1, Nq0points

        Sctrq1 = 0.0_dp

        omegaq2 = 0.0_dp
        nBEq2 = 0.0_dp
        grp_velq2 = 0.0_dp

        Ifq2Exists = 0

        q0_single = q0(1:3, i0)
        ZroCntrd = ZeroCenteredMesh(i0)

        if ( ZroCntrd ) then

            call ScatterProbAllq1( sys, FC2, FC3, EwaldConst, ph_q0, ph_ZeroCntrdq1, &
                                 & tetraQZeroCentered, T, q0_single, q1mesh, &
                                 & i0, Nq1points, Nbasis, Ndof, my_Qsize, my_offset, ZroCntrd, LongEw, &
                                 & omegaq2, nBEq2, grp_velq2, Sctrq1, Ifq2Exists )

        else

            call ScatterProbAllq1( sys, FC2, FC3, EwaldConst, ph_q0, ph_NotZeroCntrdq1, &
                                 & tetraQNotZeroCentered, T, q0_single, q1mesh, &
                                 & i0, Nq1points, Nbasis, Ndof, my_Qsize, my_offset, ZroCntrd, LongEw, &
                                 & omegaq2, nBEq2, grp_velq2, Sctrq1, Ifq2Exists )

        end if

        SYNC ALL
        call co_sum( Sctrq1, stat=err, ERRMSG=err_msg )
        if ( err /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") err_msg

        SYNC ALL
        call co_sum( omegaq2, stat=err, ERRMSG=err_msg )
        if ( err /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") err_msg

        SYNC ALL
        call co_sum( nBEq2, stat=err, ERRMSG=err_msg )
        if ( err /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") err_msg

        SYNC ALL
        call co_sum( Ifq2Exists, stat=err, ERRMSG=err_msg )
        if ( err /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") err_msg

        SYNC ALL

        !-Off-! if ( ZroCntrd ) then

        !-Off-!     call CalculateLinewidth( Ndof, Nq1points, my_QsizeT, my_offsetT, &
        !-Off-!                            & tetraQZeroCentered, ph_q0, ph_ZeroCntrdq1, omegaq2, nBEq2, &
        !-Off-!                            & Sctrq1, OnebySigma, OnebySigma2, i0, q1mesh, Ifq2Exists, lw )
        !-Off-!     
        !-Off-! else

        !-Off-!     call CalculateLinewidth( Ndof, Nq1points, my_QsizeT, my_offsetT, &
        !-Off-!                            & tetraQNotZeroCentered, ph_q0, ph_NotZeroCntrdq1, omegaq2, nBEq2, &
        !-Off-!                            & Sctrq1, OnebySigma, OnebySigma2, i0, q1mesh, Ifq2Exists, lw )
        !-Off-! end if

        if ( ZroCntrd ) then

            call CalculateLwTetrahedron( Ndof, Nq1points, my_QsizeT, my_offsetT, tetraQZeroCentered, &
                                       & ph_q0, ph_ZeroCntrdq1, omegaq2, nBEq2, Sctrq1, i0, lw )

        else

            call CalculateLwTetrahedron( Ndof, Nq1points, my_QsizeT, my_offsetT, tetraQNotZeroCentered, &
                                       & ph_q0, ph_NotZeroCntrdq1, omegaq2, nBEq2, Sctrq1, i0, lw )

        end if

        LW_q0 = lw

        SYNC ALL
        call co_sum( LW_q0, stat=err, ERRMSG=err_msg )
        if ( err /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") err_msg
        SYNC ALL

        img1_chk3: if ( this_image() == 1 ) then

            call WriteAsciiOutFile(outfile_ascii, outfile_frmt, dist_arr_q0(i0), &
                                 & ph_q0%omega(:, i0), LW_q0, Ndof, HeaderWrite)

        end if img1_chk3

    end do q0Loop

    ! ===================================== Phonon Linewidth ==================================== !

    SYNC ALL

    img1_chk10: if ( this_image() == 1 ) then

        time_LW = time_passed%elapsed_time()

        write(*, *)
        write(*, 111)
        write(*, 222)
        write(*, 666) time_2nd / 60.0_dp
        write(*, 333) time_LW / 60.0_dp
        write(*, 222)
        write(*, 111)
        write(*, *)

    end if img1_chk10

    if ( this_image() == 1) write(*, '(/)')

    666 FORMAT(10X, '|', 23X, '    Phonon calculation: ', F10.4, ' min', 24X, '|')
    333 FORMAT(10X, '|', 23X, '      Phonon Linewidth: ', F10.4, ' min', 24X, '|')

    111 FORMAT(10X,'|', 23('='), ' Time taken for different subroutines ', 24('='), '|')
    222 FORMAT(10X,'|',85X,'|')

end program main

