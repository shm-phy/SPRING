
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
    use IndxChng,                   only : indx2num
    use unit_cell,                  only : cell
    use timer_class,                only : timer
    use Symmetry,                   only : SpaceGr_typ
    use hdf5_wrap,                  only : r_RestartLW, w_RestartLW, &
                                         & w1dh5_no_grp, w2dh5_no_grp, w2dh5_no_grp_type2
    use FC2_mod,                    only : FC2type
    use FC3_mod,                    only : FC3type
    use FC4_mod,                    only : FC4type
    use EwaldMod,                   only : EwaldParam
    use TetrahedronQ,               only : q_tetrahedron
    use Helper,                     only : DistributeQpoints, WriteAsciiOutFile, &
                                         & ReadNmlKappa, ReadPointDefVar, &
                                         & ShowProcessorDistribution, PolarizationStrtEnd
    use phonon_m,                   only : Phon

    use MomentumConsrv3rd_mod,      only : MomentumConsrv3rd
    use MomentumConsrv4th_mod,      only : MomentumConsrv4th
    use ThirdOrder,                 only : ScatterProbAllq1, CalculateLinewidth, &
                                         & CalculateLwTetrahedron
    use IsoScatter,                 only : CalculateIsoW, CalculateIsoLwTetra
    use FourthOrder,                only : ScatterProbq1q2_per
    use Tetrahedron4th,             only : TetrahedronIntegration4th
    use SymmVelKronProd,            only : GrVelKronProdSymm
    use PointDefect,                only : ScatterRatePointDefect, PerturbationMatrix
    use kappa_phonon,               only : kappa_noIter, kappa_noIter_symm

#ifndef GNU
    use ifport,                     only : SYSTEM
#endif

    implicit none

    type(cell)                                              :: sys
    type(timer)                                             :: time_passed, restrt_timer
    type(SpaceGr_typ)                                       :: SpGr
    type(FC2type)                                           :: FC2, &
                                                             & FC2_pd !** Variable for point-defect calculation **!
    type(FC3type)                                           :: FC3
    type(FC4type)                                           :: FC4

    type(q_tetrahedron)                                     :: tetraQ, &
                                                             & qPnts_pd !** Variable for point-defect calculation **!

    type(EwaldParam)                                        :: EwaldConst
    type(Phon)                                              :: ph_q0, ph_q1, &
                                                             & ph_pd !** Variable for point-defect calculation **!

    character(len=256)                                      :: filename, msg, LW_Restartfile, &
                                                             & kappa_outFile
    character(len=512)                                      :: SctrEl4thRDir
    character(len=24)                                       :: Temp_char, Ndof_char
    character(len=8)                                        :: frmt
    character(len=64)                                       :: fc2file, fc3file, fc4file
    character(len=64)                                       :: outfile_ascii
    character(len=128)                                      :: outfile_frmt, RestartDir, command, &
                                                             & IFC2_file_pd, outfile_pdScatter, & !** Variable for point-defect calculation **!
                                                             & SctrEl4thWDir
    character(len=1), dimension(3)                          :: SelectionRule3
    character(len=1), dimension(4)                          :: SelectionRule4

    real(dp)                                                :: vol, timelimit
    real(dp)                                                :: T, freq_cut, &
                                                             & rho_def !** Variable for point-defect calculation **!
    real(dp)                                                :: time_2nd, time_LW

    real(dp), dimension(:), allocatable                     :: lw, g2iso, P3_q0, &
                                                             & Mass_def, tau_q0pd_inv !** Variable for point-defect calculation **!
    real(dp), dimension(:, :), allocatable                  :: V_pertbFC, V_pertbM, & !** Variable for point-defect calculation **!
                                                               ScatterAllq0_pd !** Variable for point-defect calculation **!

    real(dp)                                                :: Sigma, OnebySigma, &
                                                             & OnebySigma2, symprec, P3

    ! ..............:::::::::::::: Coindexed Objects ::::::::::::::.............. !

    real(dp), allocatable, dimension(:,:,:,:,:), codimension[:]     :: Sctrq1q2
    real(dp), allocatable, dimension(:,:,:,:), codimension[:]       :: Sctrq1
    real(dp), allocatable, dimension(:,:,:), codimension[:]         :: W_iso, KrnGrpVelSym
    real(dp), allocatable, dimension(:), codimension[:]             :: LW_q0, P_q0
    real(dp), allocatable, dimension(:,:), codimension[:]           :: LW_Allq0, P_Allq0, Spectral_k, Mode_k
    real(dp), dimension(9), codimension[*]                          :: kappa
    real(dp), dimension(1), codimension[*]                          :: PhaseSpace3 !gfrotran breaks when -finit-real=snan, 
                                                                                   !so define a rank=1, dim=1 instead of scalar
    integer, dimension(2), codimension[*]                           :: progress

    ! ..............:::::::::::::: Coindexed Objects ::::::::::::::.............. !

    integer, dimension(3)                                   :: qmesh, shift, &
                                                             & qMesh_pd !** Variable for point-defect calculation **!
    integer                                                 :: time_rev

    integer                                                 :: Nbasis, Ndof, Nq1points, &
                                                             & Nq1Perm, Nq0points, &
                                                             & Ntetrahedrons, Nq1q2Perm, &
                                                             & NTtrhdrn_pd, NQpnts_pd !** Variable for point-defect calculation **!

    integer, dimension(3)                                   :: q1mesh

    integer, dimension(2, 3)                                :: Polarization3rd
    integer, dimension(2, 4)                                :: Polarization4th

    integer                                                 :: my_Qsizeq0, my_offsetq0
    integer, dimension(2)                                   :: my_edgeq0

    integer                                                 :: my_QsizeAllq, my_offsetAllq
    integer, dimension(2)                                   :: my_edgeAllq

    integer                                                 :: my_QsizePermq, my_offsetPermq
    integer, dimension(2)                                   :: my_edgePermq

    integer                                                 :: my_QsizeT, my_offsetT
    integer, dimension(2)                                   :: my_edgeT

    integer                                                 :: my_Qsize4th, my_offset4th
    integer, dimension(2)                                   :: my_edge4th

    integer                                                 :: my_Qsize_pd, my_offsetq_pd !** Variable for point-defect calculation **!
    integer, dimension(2)                                   :: my_edgeq_pd                !** Variable for point-defect calculation **!

    integer                                                 :: qi0, s, current_indx
    integer, dimension(2)                                   :: bound, current_address

    integer                                                 :: i0, start, istat

    integer                                                 :: QPerProcsMax_4th, &
                                                             & NumDef_sup !** Variable for point-defect calculation **!
    integer, dimension(:, :), allocatable                   :: Allq1q2, q1q2_permuted, &
                                                             & q1q2q3_permuted, QinWhichProcs_4th, &
                                                             & cell_bsatm !** Variable for point-defect calculation **!

    integer, dimension(:,:,:), allocatable                  :: Allq1q2q3

    logical                                                 :: LongEW, HeaderWrite, &
                                                             & HeaderWrite_pd, & !** Variable for point-defect calculation **!
                                                             & isotope, FourPhonon, WriteRestart, &
                                                             & ReadRestart, StartFromBegin4th, TimeUp, &
                                                             & SpectralResolve, ModeResolve, Symm, &
                                                             & PointDef, MatInvFull, OptcTherm, & !** Variable for point-defect calculation **!
                                                             & WriteSctrEl4th, ReadSctrEl4th, CalSctrEl4th
    logical                                                 :: nd_Select_3rd, nd_Select_4th

    if ( this_image() == 1 ) then
        call ShowWelcomeBanner()
        call time_passed%start_timer()
    end if

    TimeUp = .false.
    call restrt_timer%start_timer()

    call get_command_line(filename, T, freq_cut, SpectralResolve, ModeResolve, Symm)
    call sys%init_data(filename)

    Nbasis = sys%natm
    vol = sys%vol !* unit cell volume *!
    Ndof = 3 * Nbasis

    frmt = '(F6.1)'
    write(Temp_char, frmt) T

    !* Restart *!
    RestartDir = 'Restart_'//trim(sys%prefix)
    LW_Restartfile = trim(adjustl(adjustr(RestartDir)))//'/'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_LW_Restart.h5'

    SctrEl4thWDir = 'ScatterMatEl4th_'//trim(sys%prefix)
    !* Restart *!

    fc2file = 'FC2nd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
    call FC2%set_FC2(fc2file)

    fc3file = 'FC3rd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
    call FC3%set_FC3(sys, fc3file)

    kappa_outFile = 'kappa_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

    outfile_ascii = 'linewidth_IBZ_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.txt'

    frmt = '(I4)'
    write(Ndof_char, frmt) (2*Ndof-1)
    outfile_frmt = "("//trim(adjustl(adjustr(Ndof_char)))//"(G20.12, '  '), G20.12)"

    ! ========================== Read the LW input from namelist file =========================== !

    allocate( g2iso(Nbasis) )

    call ReadNmlKappa( filename, sys%prefix, SelectionRule3, SelectionRule4, SctrEl4thRDir, &
                     & Sigma, symprec, timelimit, g2iso, qmesh, shift, time_rev, &
                     & isotope, FourPhonon, WriteRestart, ReadRestart, WriteSctrEl4th, &
                     & ReadSctrEl4th )
    CalSctrEl4th = .not. ReadSctrEl4th

    call PolarizationStrtEnd( 3, Ndof, SelectionRule3, Polarization3rd, nd_Select_3rd )

    if ( FourPhonon ) call PolarizationStrtEnd( 4, Ndof, SelectionRule4, Polarization4th, nd_Select_4th )

    OnebySigma = 1.0_dp / Sigma
    OnebySigma2 = (OnebySigma ** 2)

    ! ========================== Read the LW input from namelist file =========================== !

    ! =================== Read Point-Defect realted input from namelist file ==================== !

    allocate( Mass_def(Nbasis) )

    call ReadPointDefVar( filename, Temp_char, Nbasis, sys, IFC2_file_pd, Mass_def, rho_def, &
                        & qMesh_pd, NumDef_sup, cell_bsatm, MatInvFull, OptcTherm, PointDef )

    ! =================== Read Point-Defect realted input from namelist file ==================== !

    FourPhonon0: if ( FourPhonon ) then
        fc4file = 'FC4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
        call FC4%set_FC4(sys, fc4file)
    end if FourPhonon0

    ! ========================= Find space group operations from spglib ========================= !
    call SpGr%set_SpaceGrSPGlib(sys, symprec)
    ! ========================= Find space group operations from spglib ========================= !

    !* Restart *!
    command = 'mkdir -p -v ./'//trim(adjustl(adjustr(RestartDir)))

    CreateDir: if ( WriteRestart .and. (this_image() == 1) ) then

#ifdef GNU
        call SYSTEM(command, status=istat)
#else
        istat = SYSTEM(command)
#endif
        dir_create_err: if ( istat /= 0 ) then
            write(*, *) 'Directory creating failed: iostat = ', istat
        end if dir_create_err

    end if CreateDir
    !* Restart *!

    !* 4th order Scattering Matrix elements write *!
    command = 'mkdir -p -v ./'//trim(adjustl(adjustr(SctrEl4thWDir)))

    CreateDir2: if ( WriteSctrEl4th .and. (this_image() == 1) ) then

#ifdef GNU
        call SYSTEM(command, status=istat)
#else
        istat = SYSTEM(command)
#endif
        dir_create_err2: if ( istat /= 0 ) then
            write(*, *) 'Directory creating failed: iostat = ', istat
        end if dir_create_err2

    end if CreateDir2
    !* 4th order Scattering Matrix elements write *!

    ! ================================= q1 with 0-centered mesh ================================= !

    q1mesh = qmesh 
    Nq1points = product( q1mesh )
    Nq1Perm = (Nq1points / 2) + 1
    Nq1q2Perm = ((Nq1points * Nq1points) / 6) + (Nq1points / 2) + 1

    call tetraQ%make_tetrahedron( sys, q1mesh, shift, .true. )

    Ntetrahedrons = Nq1points * 6

    ! ================================= q1 with 0-centered mesh ================================= !

    ! ================================= q0 from Irreducible BZ ================================== !

    call tetraQ%FindIrrQ( sys, q1mesh, shift, time_rev, symprec )

    Nq0points = size( tetraQ%irr_q_pnt_int, 2)

    ! ================================= q0 from Irreducible BZ ================================== !

    ! ========== q-points for tetrahedron integration in Green's Function calculatiuon ========== !

    PointDefCheck1: if ( PointDef ) then

        call FC2_pd%set_FC2(IFC2_file_pd)

        allocate( tau_q0pd_inv(Ndof) )
        tau_q0pd_inv = 0.0_dp

        allocate( ScatterAllq0_pd(Ndof, Nq0points) )

        outfile_pdScatter = 'ScatterRate_pd_IBZ_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.txt'

        call qPnts_pd%make_tetrahedron( sys, qMesh_pd, shift, .true. )
        NQpnts_pd = product( qMesh_pd )
        NTtrhdrn_pd = 6 * NQpnts_pd

        call DistributeQpoints( NQpnts_pd, my_Qsize_pd, my_offsetq_pd, my_edgeq_pd ) !Tetrahedron distribution for 
                                                                                     !Greens function integration
        allocate( V_pertbFC(Ndof, Ndof), V_pertbM(Ndof, Ndof) )
        V_pertbFC = 0.0_dp
        V_pertbM = 0.0_dp

        call PerturbationMatrix( Ndof, Nbasis, Mass_def, sys, FC2, FC2_pd, V_pertbFC, V_pertbM )

    end if PointDefCheck1

    ! ========== q-points for tetrahedron integration in Green's Function calculatiuon ========== !

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

    call ph_q0%set_Phon_mp( sys, FC2, EwaldConst, tetraQ%irr_q_pnt, T, freq_cut, Ndof, LongEW, flag_pd=.false. )

    call ph_q1%set_Phon_mp( sys, FC2, EwaldConst, tetraQ%q_pnt, T, freq_cut, Ndof, LongEW, flag_pd=.false. )

    if ( PointDef ) call ph_pd%set_Phon_mp( sys, FC2, EwaldConst, qPnts_pd%q_pnt, T, freq_cut, Ndof, LongEW, flag_pd=.true. )

    SYNC ALL

    ! ======================= find phonon freq., velocity and eigenvectors ====================== !

    ! ========================== distribute the q0's through the images ========================= !
    call DistributeQpoints( Nq0points, my_Qsizeq0, my_offsetq0, my_edgeq0 )
    ! ========================== distribute the q0's through the images ========================= !

    ! ========================== distribute the q1's through the images ========================= !
    ! For each images there are my_Qsize number of qs. So the shape of the allocatable array is   !
    ! my_Qsize with index from 1 to my_Qsize. To recover the q-point add my_offset with the index !

    !*  New Distribute  *!
    call DistributeQpoints(Nq1points, my_QsizeAllq, my_offsetAllq, my_edgeAllq) !For Isotopic scattering 
                                                                                !and FourPhonon Integration
    !*  New Distribute  *!

    !*  New Distribute  *!
    call DistributeQpoints(Nq1Perm, my_QsizePermq, my_offsetPermq, my_edgePermq) !Three-Phonon scattering
    !*  New Distribute  *!

    if ( .not. ReadRestart ) &
                             & call ShowProcessorDistribution( "q-points distribution (3-ph)", &
                             &                                 my_edgePermq, my_QsizePermq, my_offsetPermq )

    ! ========================== distribute the q1's through the images ========================= !

    ! ==================== distribute the q's through the images for 4th-order ================== !

    FourPhonon1: if ( FourPhonon ) then

        allocate( QinWhichProcs_4th(2, Nq1q2Perm), STAT=istat, ERRMSG=msg )
        if ( istat /= 0 ) then
            write(*, 20) this_image(), msg
            20 FORMAT( " Memory allocation ERROR in image : ", I5, ". Error message: ", A128 )
            ERROR STOP
        end if

        call DistributeQpoints(Nq1q2Perm, my_Qsize4th, my_offset4th, my_edge4th, &
                              & QinWhichProcs_4th, QPerProcsMax_4th)

        if ( .not. ReadRestart ) & 
                                 & call ShowProcessorDistribution( "q-points distribution (4-ph)", &
                                 &                                  my_edge4th, my_Qsize4th, my_offset4th )
        SYNC ALL

        ! ** Debug ** !
        !-! if ( this_image() == 1 ) then
        !-!     write(*, *) "QPerProcsMax_4th = ", QPerProcsMax_4th
        !-!     do i0 = 1, Nq1q2Perm
        !-!         write(*, *) i0, " => ", QinWhichProcs_4th(:, i0)
        !-!     end do
        !-! end if
        ! ** Debug ** !
    end if FourPhonon1

    ! ==================== distribute the q's through the images for 4th-order ================== !

    ! ====================== distribute the tetrahedrons through the images ===================== !
    ! For each images there are my_Qsize number of qs. So the shape of the allocatable array is   !
    ! my_Qsize with index from 1 to my_Qsize. To recover the q-point add my_offset with the index !

    !*  New Distribute  *!
    call DistributeQpoints(Ntetrahedrons, my_QsizeT, my_offsetT, my_edgeT)
    !*  New Distribute  *!

    if ( .not. ReadRestart ) & 
                             & call ShowProcessorDistribution( "q-points distribution (for tetrahedron integration)", &
                             &                                  my_edgeT, my_QsizeT, my_offsetT )

    if ( this_image() == 1 ) then
        time_2nd = time_passed%elapsed_time()
        call time_passed%start_timer()
    end if

    ! ====================== distribute the tetrahedrons through the images ===================== !

    SYNC ALL
    allocate( LW_q0(Ndof)[*] )
    SYNC ALL
    LW_q0 = 0.0_dp

    allocate( P_q0(Ndof)[*] )
    SYNC ALL
    P_q0 = 0.0_dp

    allocate( LW_Allq0(Ndof, Nq0points)[*] )
    SYNC ALL
    LW_Allq0(:, :) = 0.0_dp

    allocate( P_Allq0(Ndof, Nq0points)[*] )
    SYNC ALL
    P_Allq0(:, :) = 0.0_dp

    allocate( Sctrq1(Ndof, Ndof, Ndof, Nq1points)[*], STAT=istat, ERRMSG=msg )
    if ( istat /= 0 ) then
        write(*, 25) this_image(), msg
        25 FORMAT( " Memory allocation ERROR in image : ", I5, ". Error message: ", A128 )
        ERROR STOP
    end if
    SYNC ALL

    PhaseSpace3 = 0.0_dp

    Isotope1: if ( isotope ) then
        allocate( W_iso(Ndof, Ndof, Nq1points)[*], STAT=istat, ERRMSG=msg )
        if ( istat /= 0 ) then
            write(*, 26) this_image(), msg
            26 FORMAT( " Memory allocation ERROR in image : ", I5, ". Error message: ", A128 )
            ERROR STOP
        end if
        SYNC ALL
    end if Isotope1

    FourPhonon2: if ( FourPhonon ) then

        allocate( Sctrq1q2(Ndof, Ndof, Ndof, Ndof, QPerProcsMax_4th)[*], STAT=istat, ERRMSG=msg )
        if ( istat /= 0 ) then
            write(*, 32) this_image(), msg
            32 FORMAT( " Memory allocation ERROR in image : ", I5, ". Error message: ", A128 )
            ERROR STOP
        end if
        SYNC ALL

        allocate( Allq1q2q3(3, Nq1points, Nq1points), STAT=istat, ERRMSG=msg )
        if ( istat /= 0 ) then
            write(*, 28) this_image(), msg
            28 FORMAT( " Memory allocation ERROR in image : ", I5, ". Error message: ", A128 )
            ERROR STOP
        end if

        allocate( q1q2q3_permuted(3, Nq1q2Perm), STAT=istat, ERRMSG=msg )
        if ( istat /= 0 ) then
            write(*, 30) this_image(), msg
            30 FORMAT( " Memory allocation ERROR in image : ", I5, ". Error message: ", A128 )
            ERROR STOP
        end if

    end if FourPhonon2

    ! ===================================== Phonon Linewidth ==================================== !

    RestartRead: if ( ReadRestart ) then

        if ( this_image() == 1 ) then
            call r_RestartLW(LW_Restartfile, Ndof, Nq0points, progress, LW_Allq0)
        end if

        SYNC ALL
        call co_broadcast( LW_Allq0, source_image=1, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_broadcast : ', A128 )") msg
        call co_broadcast( progress, source_image=1, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_broadcast : ', A128 )") msg
        SYNC ALL

        if ( progress(2) == 1 ) then

            start = progress(1) + 1
            StartFromBegin4th = .true.

            HeaderWrite = .false.
            HeaderWrite_pd = .false.
            
        else

            start = progress(1)
            StartFromBegin4th = .false.

            if ( start == 1 ) then
                HeaderWrite = .true.
                HeaderWrite_pd = .true.
            else
                HeaderWrite = .false.
                HeaderWrite_pd = .false.
            end if

        end if

    else RestartRead

        start = 1
        StartFromBegin4th = .true.
        HeaderWrite = .true.
        HeaderWrite_pd = .true.

    end if RestartRead

    allocate( Allq1q2(4, Nq1points) )
    allocate( q1q2_permuted(2, Nq1Perm) )
    allocate( lw(Ndof) )
    allocate( P3_q0(Ndof) )

    q0Loop: do i0 = start, Nq0points

        if ( this_image() == 1 ) then

            write(*, 55) 
            55 FORMAT( 5X, "Starting Calculation of  Third-Order Scattering-Matrix elements ...")

        end if

        Sctrq1 = 0.0_dp

        call MomentumConsrv3rd( tetraQ, Nq1points, Nq1Perm, i0, q1mesh, Allq1q2, q1q2_permuted )

        call ScatterProbAllq1( sys, FC3, ph_q0, ph_q1, tetraQ, q1mesh, &
                             & i0, Nq1points, Nq1Perm, Nbasis, Ndof, my_QsizePermq, &
                             & my_offsetPermq, q1q2_permuted, Sctrq1 )

        SYNC ALL
        call co_sum( Sctrq1, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
        SYNC ALL

        if ( this_image() == 1 ) then

            write(*, 65) i0, Nq0points
            65 FORMAT( 6X, "Calculation of  Third-Order Scattering-Matrix elements completed. Progress = ", &
                     & I6, " / ", I6 )

        end if

        call CalculateLwTetrahedron( Ndof, Nq1points, my_QsizeT, my_offsetT, &
                                   & tetraQ, ph_q0, ph_q1, Sctrq1, &
                                   & Polarization3rd, Allq1q2, i0, P3, lw, P3_q0 )
        LW_q0 = lw !** This is third-order linewidth, so no 0.5 multiplication **!
        P_q0 = P3_q0

        PhaseSpace3(1) = PhaseSpace3(1) + dble( tetraQ%multiplicity(i0) ) * P3

        IsotopeScatter: if ( isotope ) then

            W_iso = 0.0_dp

            call CalculateIsoW( ph_q0, ph_q1, i0, Nq1points, &
                              & my_QsizeAllq, my_offsetAllq, Nbasis, Ndof, g2iso, W_iso)

            SYNC ALL
            call co_sum( W_iso, stat=istat, ERRMSG=msg )
            if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
            SYNC ALL

            call CalculateIsoLwTetra( Ndof, Nq1points, my_QsizeT, my_offsetT, &
                                    & tetraQ, ph_q0, ph_q1, &
                                    & W_iso, i0, lw )

            LW_q0 = LW_q0 + (0.5_dp * lw) !** this is scattering-rate, 0.5 multiplication to convert it into linewidth **!

        end if IsotopeScatter

        FourPhonon3: if ( FourPhonon ) then

            if ( this_image() == 1 ) then

                write(*, 45) 
                45 FORMAT( 5X, "Starting Calculation of Fourth-Order Scattering-Matrix elements ...")

            end if

            Sctrq1q2 = 0.0_dp

            call MomentumConsrv4th( tetraQ, Nq1points, Nq1q2Perm, i0, q1mesh, Allq1q2q3, q1q2q3_permuted )

            call ScatterProbq1q2_per( sys, FC4, ph_q0, ph_q1, tetraQ, restrt_timer, RestartDir, SctrEl4thWDir, &
                                    & SctrEl4thRDir, timelimit, i0, Nq1points, Nq1q2Perm, Nbasis, Ndof, my_Qsize4th, &
                                    & my_offset4th, QPerProcsMax_4th, q1q2q3_permuted, WriteRestart, ReadRestart, &
                                    & StartFromBegin4th, WriteSctrEl4th, CalSctrEl4th, TimeUp, Sctrq1q2 )

            SYNC ALL

            if ( WriteRestart .and. ((restrt_timer%elapsed_time() / 60.0_dp) > timelimit) ) TimeUp = .true.

            TimeUpWrite: if ( TimeUp ) then

                if ( this_image() == 1 ) then
                    
                    call w_RestartLW(LW_Restartfile, Ndof, Nq0points, i0, 0, LW_Allq0)
                    write(*, 258)

                end if

                SYNC ALL
                ERROR STOP

            end if TimeUpWrite

            SYNC ALL
            if ( this_image() == 1 ) then

                write(*, 75) i0, Nq0points
                75 FORMAT( 6X, "Calculation of Fourth-order Scattering-Matrix elements completed. Progress = ", &
                         & I6, " / ", I6 )

            end if

            call TetrahedronIntegration4th( Nq1points, Ndof, QPerProcsMax_4th, Nq1q2Perm, &
                                          & my_QsizeAllq, my_offsetAllq, i0, tetraQ, &
                                          & ph_q0, ph_q1, Sctrq1q2, Polarization4th, Allq1q2q3, &
                                          & QinWhichProcs_4th, lw )

            LW_q0 = LW_q0 + (0.5_dp * lw) !** this is scattering-rate, 0.5 multiplication to convert it into linewidth **!

        end if FourPhonon3

        PointDefCheck2: if ( PointDef ) then

            call ScatterRatePointDefect( i0, Ndof, NTtrhdrn_pd, my_Qsize_pd, my_offsetq_pd, vol, rho_def, &
                                       & V_pertbFC, V_pertbM, tetraQ, qPnts_pd, ph_q0, ph_pd, &
                                       & MatInvFull, OptcTherm, tau_q0pd_inv )

            if ( this_image() == 1 ) call WriteAsciiOutFile( outfile_pdScatter, outfile_frmt, &
                                                           & ph_q0%omega(:, i0), tau_q0pd_inv, Ndof, HeaderWrite_pd )
            ScatterAllq0_pd(:, i0) = tau_q0pd_inv

        end if PointDefCheck2

        SYNC ALL
        call co_sum( LW_q0, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

        SYNC ALL
        call co_sum( P_q0, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
        SYNC ALL

        if ( PointDef ) LW_q0 = LW_q0 + (0.5_dp * tau_q0pd_inv)
                                        !** this is scattering-rate, 0.5 multiplication to convert it into linewidth **!

        LW_Allq0(:, i0) = LW_q0(:)
        P_Allq0(:, i0) = P_q0(:) * 1000.0_dp / ( 9.0_dp * dble( Ndof**2 ) * dble( Nq1points ) ) ! in fs

        img1_chk3: if ( this_image() == 1 ) then

            call WriteAsciiOutFile(outfile_ascii, outfile_frmt, &
                                 & ph_q0%omega(:, i0), LW_q0, Ndof, HeaderWrite)

        end if img1_chk3

        StartFromBegin4th = .true.

        SYNC ALL
        RestartWrite2: if ( WriteRestart .and. &
                         & ( (restrt_timer%elapsed_time() / 60.0_dp) > timelimit ) ) then

            if ( this_image() == 1 ) then
                call w_RestartLW(LW_Restartfile, Ndof, Nq0points, i0, 1, LW_Allq0)
                write(*, 258)
                TimeUp = .true.
            end if

            SYNC ALL
            ERROR STOP

        end if RestartWrite2

    end do q0Loop

    ! ===================================== Phonon Linewidth ==================================== !

    ! ================================== Scattering Phase Space ================================= !

    SYNC ALL
    call co_sum( PhaseSpace3, stat=istat, ERRMSG=msg )
    if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
    SYNC ALL
    
    PhaseSpace3 = PhaseSpace3 * 1000.0_dp / ( 9.0_dp * dble( Ndof**3 ) * dble( Nq1points**2 ) ) ! in fs

    ! ================================== Scattering Phase Space ================================= !

    ! ========================== Thermal conductivity (non iterative)  ========================== !

    SymmVelCheck: if ( Symm ) then

        call GrVelKronProdSymm( Nq0points, Ndof, my_Qsizeq0, my_offsetq0, &
                              & SpGr, tetraQ, ph_q0, KrnGrpVelSym )

        SYNC ALL
        call co_sum( KrnGrpVelSym, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
        SYNC ALL

        call kappa_noIter_symm( Nq1points, Nq0points, Ndof, my_Qsizeq0, my_offsetq0, &
                              & tetraQ, ph_q0, vol, T, LW_Allq0, KrnGrpVelSym, Polarization3rd, &
                              & Polarization4th, SpectralResolve, ModeResolve, FourPhonon, &
                              & Spectral_k, Mode_k, kappa )

    else SymmVelCheck

        call kappa_noIter( Nq1points, Nq0points, Ndof, my_QsizeAllq, my_offsetAllq, &
                         & tetraQ, ph_q1, vol, T, LW_Allq0, Polarization3rd, Polarization4th, &
                         & SpectralResolve, ModeResolve, FourPhonon, Spectral_k, Mode_k, kappa )

    end if SymmVelCheck

    SYNC ALL
    call co_sum( kappa, stat=istat, ERRMSG=msg )
    if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

    SYNC ALL
    if ( SpectralResolve ) then
        call co_sum( Spectral_k, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
    end if

    SYNC ALL
    if ( ModeResolve ) then
        call co_sum( Mode_k, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
    end if

    ! ========================== Thermal conductivity (non iterative)  ========================== !

    img1_chk10: if ( this_image() == 1 ) then

        call w1dh5_no_grp(kappa, kappa_outFile)
        call w2dh5_no_grp_type2( Ndof, Nq0points, ph_q0%omega, LW_Allq0, kappa_outFile, "LineWidth" )
        call w2dh5_no_grp_type2( Ndof, Nq0points, ph_q0%omega, P_Allq0, kappa_outFile, "PhaseSpace3rd" )
        if ( PointDef ) call w2dh5_no_grp_type2( Ndof, Nq0points, ph_q0%omega, ScatterAllq0_pd, kappa_outFile, "Scatter_pd" )

        Spectral: if ( SpectralResolve ) then

            bound = (/Nq0points, Ndof/)
            q0loop2: do qi0 = 1, Nq0points
                s_loop: do s = 1, Ndof

                    current_address = (/qi0, s/)
                    current_indx = indx2num( bound, current_address )
                    Spectral_k(1, current_indx) = dble( s )
                    Spectral_k(2, current_indx) = ph_q0%omega(s, qi0)

                end do s_loop
            end do q0loop2

            call w2dh5_no_grp(Spectral_k, kappa_outFile, "SpectralResolvedKappa")

        end if Spectral

        if ( ModeResolve ) call w2dh5_no_grp(Mode_k, kappa_outFile, "ModeResolvedKappa")

        if ( nd_Select_3rd ) write(*, 264) SelectionRule3
        264 FORMAT( ' ********  You have chosen non-defult selection rule: (', 2(A2, ' '), A2, ' ) ******** ' )

        if ( FourPhonon .and. nd_Select_4th ) write(*, 244) SelectionRule4
        244 FORMAT( ' ******** (FourPhonon) You have chosen non-defult selection rule: (', 3(A2, ' '), A2, ' ) ******** ' )

        write(*, *)
        write(*, "( 'Thermal conductivity information written in file: ', A64 )") kappa_outFile

        write(*, *)
        write(*, "( 37X, '|kxx = ', G20.12, ', kxy = ', G20.12, ', kxz = ', G20.12, '|' )") kappa(1:3)
        write(*, "( 37X, '|', 82X, '|')")
        write(*, 265) kappa(4:6)
        265 FORMAT( 'Thermal conductivity tensor (W/m.K): ', '|kyx = ', G20.12, ', kyy = ', G20.12, ', kyz = ', G20.12, '|' )
        write(*, "( 37X, '|', 82X, '|')")
        write(*, "( 37X, '|kzx = ', G20.12, ', kzy = ', G20.12, ', kzz = ', G20.12, '|' )") kappa(7:9)
        write(*, *)

        write(*, 285) PhaseSpace3(1)
        285 FORMAT( '3-phonon scattering phase-space: ', G20.12, ' fs' )

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

    258 FORMAT( "RESTART: Timelimit exceeded, Restart files written, STOPPING ..." )

end program main

