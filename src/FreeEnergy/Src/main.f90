
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
    use FC4_mod,                    only : FC4type
    use EwaldMod,                   only : EwaldParam
    use CreateMesh,                 only : mesh_points, q_points_highsymm
    use DistributeQs,               only : DistributeQpoints, &
                                         & ShowProcessorDistribution
    use phonon_m,                   only : Phon

    use Irr_q_point,                only : q_points_data
    use HarmonicFreeEnergy,         only : FreeEnergy2nd
    use FreeEnergy3rdTerm1,         only : AllQ
    use FreeEnergy3rdTerm2,         only : ScatterMatAllqq2, FindFreeEngergy3rdT2
    use FreeEnergy4th,              only : FindFreeEnergy4th

    implicit none

    type(cell)                                              :: sys
    type(timer)                                             :: time_passed
    type(FC2type)                                           :: FC2
    type(FC3type)                                           :: FC3
    type(FC4type)                                           :: FC4

    type(q_points_data)                                     :: Qpoints
    type(q_points_data)                                     :: Qpnt_Har

    type(Phon)                                              :: phonon_dat
    type(Phon)                                              :: phonon_Har
    type(EwaldParam)                                        :: EwaldConst

    type(AllQ)                                              :: my_data

    complex(dp)                                             :: FreeEng3rd_t2
    complex(dp)                                             :: FreeEng4th

    real(dp)                                                :: t_phonon, t_harmonic, &
                                                             & t_3rd, t_4th

    real(dp)                                                :: FreeEng_h
    real(dp)                                                :: FreeEng3rd_t1
    real(dp)                                                :: vol
    real(dp)                                                :: T, freq_cut

    ! ..............:::::::::::::: Coindexed Object ::::::::::::::.............. !
    real(dp), dimension(1), codimension[*]                          :: FreeEng2_Large
    real(dp), dimension(1), codimension[*]                          :: FreeEng2
    complex(dp), dimension(1), codimension[*]                       :: FreeEng3
    complex(dp), dimension(1), codimension[*]                       :: FreeEng4
    complex(dp), allocatable, dimension(:,:,:,:), codimension[:]    :: Phiqq2
    ! ..............:::::::::::::: Coindexed Object ::::::::::::::.............. !
    
    character(len=256)                                      :: filename
    character(len=24)                                       :: Temp_char
    character(len=8)                                        :: frmt
    character(len=64)                                       :: fc2file, fc3file, fc4file

    integer, dimension(3)                                   :: mesh, mesh_Har, shift
    integer                                                 :: time_rev

    integer                                                 :: Nqpoints, Nqpoints_Har, &
                                                             & num_irr_q, Nbasis, Ndof
    integer                                                 :: my_Qsize, my_offset
    integer, dimension(2)                                   :: my_edge
    !integer                                                 :: err
    !integer                                                 :: istat

    integer                                                 :: Num_qq1, qMeshL
    integer                                                 :: i0, i1, count_qq1
    integer, dimension(:, :), allocatable                   :: All_qq1
    integer, dimension(:, :), allocatable                   :: my_qq1

    logical                                                 :: LongEW, AnharFreeEng, HarmonicLarge

    if ( this_image() == 1 ) then
        call ShowWelcomeBanner()
        call time_passed%start_timer()
    end if

    call get_command_line(filename, T, freq_cut, qMeshL, AnharFreeEng, HarmonicLarge)
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

    fc4file = 'FC4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'
    call FC4%set_FC4(sys, fc4file)

    !! ****************** pre-calculate the Ewald Constants for polar materials ***************** !!

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

    !! ****************** pre-calculate the Ewald Constants for polar materials ***************** !!

    !! *********************** Harmonic Free-Energy from larger Mesh-grid *********************** !!

    mesh_Har(:) = qMeshL
    Nqpoints_Har = product( mesh_Har )

    shift(:) = 0
    time_rev = 0

    call Qpnt_Har%make_tetrahedron(sys, mesh_Har, shift, time_rev)

    call phonon_Har%set_Phon_mp( sys, FC2, EwaldConst, Qpnt_Har%q_pnt, T, freq_cut, Ndof, LongEW )

    if ( this_image() == 1 ) then
        t_phonon = time_passed%elapsed_time()
        call time_passed%start_timer()
    end if

    SYNC ALL

    !*  New Distribute  *!
    call DistributeQpoints( Nqpoints_Har, my_Qsize, my_offset, my_edge )
    !*  New Distribute  *!

    call FreeEnergy2nd( Qpnt_Har, Nbasis, my_Qsize, my_offset, &
                      & phonon_Har, T, FreeEng_h )

    FreeEng2_Large = FreeEng_h

    SYNC ALL
    call co_sum( FreeEng2_Large )
    SYNC ALL

    img1_chk1: if ( this_image() == 1 ) then

        write(*, *)
        call execute_command_line(' ')

        write(*, 505) mesh_Har, FreeEng2_Large(1)
        call execute_command_line(' ')

        write(*, *)
        call execute_command_line(' ')

        505 FORMAT("Harmonic Free Energy ( ", 2(I3, ' X '), I3, " q-mesh ): ", E18.10, " eV")

    end if img1_chk1

    SYNC ALL

    call Qpnt_Har%clean_Qdata()
    call phonon_Har%clean_Phon()

    !! *********************** Harmonic Free-Energy from larger Mesh-grid *********************** !!

    mesh = sys%sup_cell
    Nqpoints = product( mesh )

    shift(:) = 0
    time_rev = 0

    call Qpoints%make_tetrahedron(sys, mesh, shift, time_rev)
    num_irr_q = size( Qpoints%irr_q_int, 2 )

    ! ======================= find phonon freq., velocity and eigenvectors ====================== !

    call phonon_dat%set_Phon_mp( sys, FC2, EwaldConst, Qpoints%q_pnt, T, freq_cut, Ndof, LongEW )

    SYNC ALL

    ! ======================= find phonon freq., velocity and eigenvectors ====================== !

    ! =============================== Harmonic part of Free Energy ============================== !

    !*  New Distribute  *!
    call DistributeQpoints( Nqpoints, my_Qsize, my_offset, my_edge )
    !*  New Distribute  *!

    call FreeEnergy2nd( Qpoints, Nbasis, my_Qsize, my_offset, &
                      & phonon_dat, T, FreeEng_h )

    FreeEng2 = FreeEng_h

    SYNC ALL
    call co_sum( FreeEng2 )
    SYNC ALL

    img1_chk2: if ( this_image() == 1 ) then

        t_harmonic = time_passed%elapsed_time()
        call time_passed%start_timer()

        write(*, *)
        call execute_command_line(' ')

        write(*, 500) mesh, FreeEng2(1)
        call execute_command_line(' ')

        write(*, *)
        call execute_command_line(' ')

        500 FORMAT("Harmonic Free Energy ( ", 2(I3, ' X '), I3, " q-mesh ): ", E18.10, " eV")

    end if img1_chk2

    SYNC ALL

    ! =============================== Harmonic part of Free Energy ============================== !

    FreeEng3 = dcmplx(0.0_dp, 0.0_dp)
    FreeEng4 = dcmplx(0.0_dp, 0.0_dp)
    t_3rd = 0.0_dp
    t_4th = 0.0_dp

    AnHarmonic: if ( AnharFreeEng ) then

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~ First Term of Free Energy Third Order ~~~~~~~~~~~~~~~~~~~~~~~~~~ !

        ! ========================== distribute the q's through the images ========================== !
        ! For each images there are my_Qsize number of qs. So the shape of the allocatable array is   !
        ! my_Qsize with index from 1 to my_Qsize. To recover the q-point add my_offset with the index !

        !*  New Distribute  *!
        call DistributeQpoints(num_irr_q, my_Qsize, my_offset, my_edge)
        !*  New Distribute  *!

        call ShowProcessorDistribution( "q-point distribution for Free Energy 3rd (Term-1)", &
                                      & my_edge, my_Qsize, my_offset )

        ! ========================== distribute the q's through the images ========================== !

        if ( this_image() == 1 ) write(*, "('Calculating Free Energy 3rd-order (term 1)')")
        call my_data%QconsrvMomentum(Qpoints, mesh, my_Qsize, my_offset, &
                                   & sys, FC3, phonon_dat, FreeEng3rd_t1) 

        FreeEng3 = FreeEng3 + dcmplx( FreeEng3rd_t1 )
        SYNC ALL

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~ First Term of Free Energy Third Order ~~~~~~~~~~~~~~~~~~~~~~~~~~ !


        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~ Second Term of Free Energy Third Order ~~~~~~~~~~~~~~~~~~~~~~~~~ !

        ! =============================== coindexed object allocation =============================== !
        allocate( Phiqq2(Ndof, Ndof, Ndof, Nqpoints)[*] )
        SYNC ALL
        Phiqq2 = dcmplx(0.0_dp, 0.0_dp)

        ! =============================== coindexed object allocation =============================== !

        !*  New Distribute  *!
        call DistributeQpoints(Nqpoints, my_Qsize, my_offset, my_edge)
        !*  New Distribute  *!

        call ShowProcessorDistribution( "q-point distribution for Scatter Mat. elements (Free Energy Term-2)", &
                                      & my_edge, my_Qsize, my_offset )

        if ( this_image() == 1 ) then

            write(*, 120)
            call execute_command_line(' ')
            
            write(*, *)
            call execute_command_line(' ')

            120 FORMAT("Calculating Scattering Matrix elements for Free Energy 3rd Order (term 2)")
        end if

        call ScatterMatAllqq2( Qpoints, mesh, Nbasis, my_Qsize, my_offset, &
                             & sys, FC3, phonon_dat, Phiqq2 )

        SYNC ALL
        call co_sum( Phiqq2 )
        SYNC ALL

        Num_qq1 = Nqpoints * Nqpoints

        allocate( All_qq1(2, Num_qq1) )

        count_qq1 = 1
        q_loop: do i0 = 1, Nqpoints
            q1_loop: do i1 = 1, Nqpoints

                All_qq1(:, count_qq1) = (/i0, i1/)
                count_qq1 = count_qq1 + 1

            end do q1_loop
        end do q_loop

        !*  New Distribute  *!
        call DistributeQpoints(Num_qq1, my_Qsize, my_offset, my_edge)
        !*  New Distribute  *!

        call ShowProcessorDistribution( "q-point distribution for Free Energy 3rd (Term-2) and 4th order", &
                                      & my_edge, my_Qsize, my_offset )

        allocate( my_qq1(2, my_Qsize) )
        my_qq1 = All_qq1( :, my_edge(1) : my_edge(2) )
        deallocate( All_qq1 )

        call FindFreeEngergy3rdT2( Qpoints, mesh, Nbasis, my_Qsize, my_qq1, &
                                 & phonon_dat, Phiqq2, FreeEng3rd_t2 )

        FreeEng3 = FreeEng3 + FreeEng3rd_t2

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~ Second Term of Free Energy Third Order ~~~~~~~~~~~~~~~~~~~~~~~~~ !
        
        SYNC ALL
        call co_sum( FreeEng3 )
        SYNC ALL

        img1_chk9: if ( this_image() == 1 ) then

            write(*, *)
            call execute_command_line(' ')

            write(*, 50) dble(FreeEng3(1))
            call execute_command_line(' ')

            write(*, 70) dimag(FreeEng3(1))
            call execute_command_line(' ')

            write(*, *)
            call execute_command_line(' ')

            50 FORMAT("Free enegry 3-rd order: ", E14.6, " eV")

            70 FORMAT("Imaginary part of Free Energy 3rd Order( Ideally 0) : ", E14.6)

            t_3rd = time_passed%elapsed_time()
            call time_passed%start_timer()

        end if img1_chk9

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Free Energy Fourth Order ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
        if ( this_image() == 1 ) then

            write(*, 130)
            call execute_command_line(' ')
            write(*, *)
            call execute_command_line(' ')

            130 FORMAT("Calculating Scattering Matrix elements for Free Energy 4th Order ")

        end if

        call FindFreeEnergy4th( Qpoints, mesh, Nbasis, my_Qsize, my_qq1, &
                              & sys, FC4, phonon_dat, FreeEng4th )

        FreeEng4 = FreeEng4th
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Free Energy Fourth Order ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

        SYNC ALL
        call co_sum( FreeEng4 )
        SYNC ALL

    end if AnHarmonic

    img1_chk10: if ( this_image() == 1 ) then

        t_4th = time_passed%elapsed_time()

        write(*, *)
        write(*, 11)
        write(*, 22)
        write(*, 32) vol
        write(*, 33) FreeEng2(1)
        write(*, 44) dble( FreeEng3(1) )
        write(*, 55) dble( FreeEng4(1) )
        write(*, 22)
        write(*, 11)
        write(*, *)

        write(*, 75) dimag(FreeEng3(1))
        write(*, 80) dimag(FreeEng4(1))

        75 FORMAT(5X, "Imaginary part of Free Energy 3rd Order( Ideally 0) : ", E14.6)
        80 FORMAT(5X, "Imaginary part of Free Energy 4th Order( Ideally 0) : ", E14.6)

        write(*, *)
        write(*, 111)
        write(*, 222)
        write(*, 666) t_phonon / 60.0_dp
        write(*, 333) t_harmonic 
        write(*, 444) t_3rd / 60.0_dp
        write(*, 555) t_4th / 60.0_dp
        write(*, 222)
        write(*, 111)
        write(*, *)

    end if img1_chk10

    if ( this_image() == 1) write(*, '(/)')

    11 FORMAT(10X,'|', 23('='), ' Free Energy Distribution ', 24('='), '|')
    22 FORMAT(10X,'|',73X,'|')

    32 FORMAT(10X, '|', 6X, '                   Volume of unit-cell: ', E18.10, " A^3", 5X, '|')
    33 FORMAT(10X, '|', 6X, '  Harmonic contribution to Free-Energy: ', E18.10, " eV", 6X, '|')
    44 FORMAT(10X, '|', 6X, ' 3rd-Order contribution to Free-Energy: ', E18.10, " eV", 6X, '|')
    55 FORMAT(10X, '|', 6X, ' 4th-Order contribution to Free-Energy: ', E18.10, " eV", 6X, '|')

    666 FORMAT(10X, '|', 23X, '    Phonon calculation: ', F10.4, ' min', 24X, '|')
    333 FORMAT(10X, '|', 23X, '  Harmonic Free-Energy: ', F10.4, ' sec', 24X, '|')
    444 FORMAT(10X, '|', 23X, ' 3rd-Order Free-Energy: ', F10.4, ' min', 24X, '|')
    555 FORMAT(10X, '|', 23X, ' 4th-Order Free-Energy: ', F10.4, ' min', 24X, '|')

    111 FORMAT(10X,'|', 23('='), ' Time taken for different subroutines ', 24('='), '|')
    222 FORMAT(10X,'|',85X,'|')

end program main

