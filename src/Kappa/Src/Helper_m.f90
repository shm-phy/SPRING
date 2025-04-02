
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

    use kinds,          only : dp
    use constants,      only : EPS, MAX_DEF
    use unit_cell,      only : cell

    implicit none
    private

    public          :: DistributeQpoints, ShowProcessorDistribution, &
                     & WriteAsciiOutFile, unique_sort, ReadNmlKappa, &
                     & ReadPointDefVar, FindPositionAndMultiplicity, &
                     & PolarizationStrtEnd, kron_prod

contains

    subroutine ReadNmlKappa( filename, prefix, SelectionRule3, SelectionRule4, SctrEl4thRDir, &
                           & Sigma, symprec, timelimit, g2iso, qmesh, shift, time_rev, &
                           & isotope, FourPhonon, WriteRestart, ReadRestart, WriteSctrEl4th, &
                           & ReadSctrEl4th )

        implicit none

        character(len=*), intent(in)                :: filename, prefix
        character(len=1), dimension(3), intent(out) :: SelectionRule3
        character(len=1), dimension(4), intent(out) :: SelectionRule4
        character(len=512), intent(out)             :: SctrEl4thRDir

        real(dp), intent(out)                       :: Sigma, symprec, timelimit
        real(dp), dimension(:), intent(out)         :: g2iso
        integer, dimension(3), intent(out)          :: qmesh, shift
        integer, intent(out)                        :: time_rev
        logical, intent(out)                        :: isotope, FourPhonon, WriteRestart, &
                                                     & ReadRestart, WriteSctrEl4th, ReadSctrEl4th

        namelist   /AnhInfo/   SelectionRule3, SelectionRule4, SctrEl4thRDir, Sigma, symprec, timelimit, g2iso, &
                             & qmesh, shift, time_rev, isotope, FourPhonon, WriteRestart, ReadRestart, &
                             & WriteSctrEl4th, ReadSctrEl4th

        ! ================================ Local Variables  ================================ !

        character(len=256)                                      :: err_msg
        integer                                                 :: err

        ! ================================ Local Variables  ================================ !

        ! ** default values ** !
        SelectionRule3(:) = (/'D', 'D', 'D'/)
        SelectionRule4(:) = (/'D', 'D', 'D', 'D'/)
        SctrEl4thRDir = './ScatterMatEl4th_'//trim(prefix)
        Sigma = 5.1_dp
        symprec = 1.0E-7_dp
        timelimit = 365.0_dp * 24.0_dp * 60.0_dp !Minute
        g2iso(:) = 0.0_dp
        qmesh = (/7, 7, 7/)
        shift = (/0, 0, 0/)
        time_rev = 0
        isotope = .false.
        FourPhonon = .false.
        WriteRestart = .false.
        ReadRestart = .false.
        WriteSctrEl4th = .false.
        ReadSctrEl4th = .false.
        ! ** default values ** !

        open(unit=5, file=filename, status='OLD', iostat=err, iomsg=err_msg, &
             action='READ', delim='APOSTROPHE')

            if ( (err /= 0) .and. (this_image() == 1) ) then
                write(*, 100) filename
                write(*, *) err_msg
                ERROR STOP
                100 FORMAT(" ERROR in opening namelist file: ", A128, "ERROR MESSAGE ==> ")
            end if

            !-! open_chk: if ( err /= 0 ) then
            !-!     write(*, 300) err
            !-!     300 FORMAT('Input file OPEN failed: iostat =', I3)
            !-!     write(*, 400) err_msg
            !-!     400 FORMAT('Error message = ', A)
            !-! end if open_chk

            read(unit=5, nml=AnhInfo, iostat=err, iomsg=err_msg)        ! ** !

            if ( (err /= 0) .and. (this_image() == 1) ) then
                write(*, 101) filename
                write(*, *) err_msg
                ERROR STOP
                101 FORMAT(" ERROR in reading AnhInfo namelist in file: ", A128, "ERROR MESSAGE ==> ")
            end if

        close(unit=5, status='KEEP', iostat=err, iomsg=err_msg)

        if ( (err /= 0) .and. (this_image() == 1) ) then
            write(*, 102) filename
            write(*, *) err_msg
            ERROR STOP
            102 FORMAT(" ERROR in closing namelist file: ", A128, "ERROR MESSAGE ==> ")
        end if

        if ( WriteSctrEl4th .and. ReadSctrEl4th .and. (this_image()==1) ) write(*, 400) 
        400 FORMAT( 'ERROR: both WriteSctrEl4th and ReadSctrEl4th can not be .true.' )
        if ( ReadSctrEl4th .and. WriteRestart .and. (this_image()==1) ) write(*, 401)
        401 FORMAT( 'WARNING: No Restart information will be written if ReadSctrEl4th = .true.' )
        if ( (.not. FourPhonon) .and. WriteRestart .and. (this_image()==1) ) write(*, 402)
        402 FORMAT( 'WARNING: Restart information writing is only supported for FourPhonon = .true.' )
        if ( (.not. FourPhonon) .and. WriteSctrEl4th .and. (this_image()==1) ) write(*, 403)
        403 FORMAT( 'WARNING: Writing 4th order scattering matrix elements is supported only for FourPhonon = .true.' )
        if ( (.not. FourPhonon) .and. ReadSctrEl4th .and. (this_image()==1) ) write(*, 404)
        404 FORMAT( 'WARNING: Reading 4th order scattering matrix elements is supported only for FourPhonon = .true.' )

    end subroutine ReadNmlKappa


    subroutine ReadPointDefVar( filename, Temp_char, Nbasis, sys, IFC2_file_pd, Mass_def, rho_def, &
                              & qMesh_pd, NumDef_sup, cell_bsatm, MatInvFull, OptcTherm, PointDef )

        implicit none

        real(dp), parameter         :: FILL=99999.0_dp

        character(len=*), intent(in)                        :: filename, Temp_char
        integer, intent(in)                                 :: Nbasis
        type(cell), intent(in)                              :: sys

        character(len=128), intent(out)                     :: IFC2_file_pd
        real(dp), dimension(Nbasis), intent(out)            :: Mass_def
        real(dp), intent(out)                               :: rho_def

        integer, dimension(3), intent(out)                  :: qMesh_pd
        integer, intent(out)                                :: NumDef_sup
        integer, dimension(:, :), allocatable, intent(out)  :: cell_bsatm

        logical, intent(out)                                :: MatInvFull, OptcTherm, PointDef

        ! ================================ Local Variables  ================================ !
        character(len=256)                      :: err_msg
        real(dp), dimension(Nbasis)             :: Mass_def_c
        integer                                 :: err, mu, cnt
        integer, dimension(Nbasis)              :: def_pos
        integer, dimension(4, MAX_DEF)          :: cell_atom
        ! ================================ Local Variables  ================================ !

        namelist   /PointDefInfo/  IFC2_file_pd, Mass_def, rho_def, qMesh_pd, def_pos, NumDef_sup, &
                                 & cell_atom, MatInvFull, OptcTherm

        !** Default Values **!
        IFC2_file_pd = 'FC2nd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F_pd.h5'
        Mass_def(:) = FILL
        rho_def = 1.0E-6_dp !in A^-3 (=10^18 cm^-3)
        qMesh_pd = (/7, 7, 7/)
        def_pos = 0
        NumDef_sup = 0
        cell_atom(:, :) = int( FILL )
        MatInvFull = .true.
        OptcTherm = .true.
        !** Default Values **!

        PointDef = .false.

        open(unit=5, file=filename, status='OLD', iostat=err, iomsg=err_msg, &
             action='READ', delim='APOSTROPHE')

            if ( (err /= 0) .and. (this_image() == 1) ) then
                write(*, 100) filename
                write(*, *) err_msg
                ERROR STOP
                100 FORMAT(" ERROR in opening namelist file: ", A128, "ERROR MESSAGE ==> ")
            end if

            read(unit=5, nml=PointDefInfo, iostat=err, iomsg=err_msg)        ! ** !

            if (err /= 0)  then

                PointDef = .false.

                if ( this_image() == 1 ) then
                    write(*, *)
                    write(*, 101) filename
                    write(*, 107) err_msg
                    write(*, 103) PointDef
                    write(*, *)
                    101 FORMAT( " No PointDefInfo detected in namelist file or variables not provided correctly: ", A56 )
                    107 FORMAT( " Error message: ", A128 )
                    103 FORMAT( " Point-defect calculation will be set to: ", L7 )
                end if

            else

                PointDef = .true.
                if ( this_image() == 1 ) write(*, 105) qMesh_pd
                105 FORMAT( "q-grid for the tetrahedron integration in Green's function: ( ", 2(I3, ' X'), I3, " )" )

            end if

        close(unit=5, status='KEEP', iostat=err, iomsg=err_msg)

        if ( (err /= 0) .and. (this_image() == 1) ) then
            write(*, 102) filename
            write(*, *) err_msg
            ERROR STOP
            102 FORMAT(" ERROR in closing namelist file: ", A128, "ERROR MESSAGE ==> ")
        end if

        !** Check for correctness of point defect variables **!
        ChecckPdInput: if ( PointDef ) then

            out_check: if ( all(def_pos == 0 ) ) then
                write(*, 200)
                write(*, 230)
                write(*, 200)
                ERROR STOP

            else out_check

                Mass_def_c = Mass_def
                cnt = 0
                basis_loop: do mu = 1, Nbasis
                    Mass_def(mu) = sys%mass(mu)
                    if ( (def_pos(mu) /= 0) ) then

                        cnt = cnt + 1
                        in_check1: if ( dabs(Mass_def_c(cnt) - FILL) < EPS ) then
                            write(*, 200) 
                            write(*, 210) cnt
                            write(*, 220) (cnt-1)
                            write(*, 200) 
                            ERROR STOP
                        else in_check1
                            Mass_def(mu) = Mass_def_c(cnt)
                            Mass_def_c(cnt) = 0.0_dp
                        end if in_check1

                    end if
                end do basis_loop

            end if out_check

            do mu = 1, Nbasis
                if ( .not. ((Mass_def_c(mu) < EPS) .or. (dabs(Mass_def_c(mu)-FILL) < EPS) ) ) then
                    write(*, 200)
                    write(*, 240)
                    write(*, 200)
                    ERROR STOP
                end if
            end do

            if ( NumDef_sup == 0 ) then
                write(*, 200)
                write(*, 250)
                write(*, 200)
                ERROR STOP
            else

                allocate( cell_bsatm(4, NumDef_sup) )
                do mu = 1, NumDef_sup
                    if ( .not. all( (cell_atom(:, mu) - int(FILL)) == 0 ) ) then
                        cell_bsatm(:, mu) = cell_atom(:, mu)
                    else
                        write(*, 200)
                        write(*, 260)
                        write(*, 200)
                        ERROR STOP
                    end if
                end do

            end if

        else ChecckPdInput

            !** Just to avoid error of unallocated **!
            NumDef_sup = 1
            allocate( cell_bsatm(4, NumDef_sup) )
            cell_bsatm(:, 1) = int(FILL)

        end if ChecckPdInput

        200 FORMAT( "=============== ERROR(in PointDefInfo) ===============")
        210 FORMAT( 5X, "Defect position (def_pos) is specified for ", I3, " defects.")
        220 FORMAT( 5X, "Defect mass (Mass_def) is given for ", I3, " positions." )
        230 FORMAT( 5X, "def_pos variable (mandatory) is not provided for point defect valculation " )
        240 FORMAT( 5X, "Mismatch between def_pos and Mass_def" )
        250 FORMAT( 5X, "Number of defect sites in supercell (NumDef_sup) is 0 " )
        260 FORMAT( 5X, "Number of defect sites in supercell (NumDef_sup) and 'cell_atom' of defects do not match")
        !** Check for correctness of point defect variables **!

    end subroutine ReadPointDefVar


    subroutine PolarizationStrtEnd( order, Ndof, character_arr, out_int_arr, nd_selection )

        implicit none

        integer, intent(in)                                 :: order, Ndof
        character(len=1), dimension(order), intent(in)      :: character_arr
        integer, dimension(2, order), intent(out)           :: out_int_arr
        logical, intent(out)                                :: nd_selection

        ! ================================ Local Variables ================================ !
        character(len=1)                    :: Branch
        integer                             :: oi
        ! ================================ Local Variables ================================ !

        ! ** Default ** !
        out_int_arr(1, :) = 1
        out_int_arr(2, :) = Ndof
        nd_selection = .false.
        ! ** Default ** !

        order_loop: do oi = 1, order

            Branch = character_arr(oi)

            SelectPolarization: if ( (Branch == 'A') .or. (Branch == 'a') ) then
                out_int_arr(1, oi) = 1
                out_int_arr(2, oi) = 3
                nd_selection = nd_selection .or. .true.

            else if ( (Branch == 'O') .or. (Branch == 'o') ) then SelectPolarization
                out_int_arr(1, oi) = 4
                out_int_arr(2, oi) = Ndof
                nd_selection = nd_selection .or. .true.

            else SelectPolarization
                out_int_arr(1, oi) = 1
                out_int_arr(2, oi) = Ndof
                nd_selection = nd_selection .or. .false.

            end if SelectPolarization

        end do order_loop

    end subroutine PolarizationStrtEnd


    subroutine DistributeQpoints(Nqs, my_Qsize, my_offset, my_edge, QinWhichProcs, QPerProcsMax)

        implicit none

        integer, intent(in)                                 :: Nqs
        integer, intent(out)                                :: my_Qsize, my_offset
        integer, dimension(2), intent(out)                  :: my_edge

        integer, dimension(2, Nqs), intent(out), optional   :: QinWhichProcs
        integer, intent(out), optional                      :: QPerProcsMax

        ! ================================= Local variable =================================== !

        character(len=128)                                  :: msg
        integer, dimension(:), allocatable                  :: QinEachProcs
        integer                                             :: Nprocs, div_int, rem, i, istat

        ! ================================= Local variable =================================== !

        Nprocs = num_images()
        if ( present(QinWhichProcs) ) QinWhichProcs(:, :) = 0

        allocate( QinEachProcs( Nprocs ) )
        QinEachProcs = 0

        OnlyImage1: if ( this_image() == 1 ) then

            div_int = Nqs / Nprocs
            rem = mod( Nqs, Nprocs )

            Warn: if ( div_int == 0 ) then

                write(*, 25) Nqs, Nprocs
                25 FORMAT('Number of q-points ( ', I5, ' ) is smaller than number of processors ( ', I5, ' ) ')

            end if Warn

            QinEachProcs( : ) = div_int
            QinEachProcs( 1:rem ) = QinEachProcs( 1:rem ) + 1

            Check: if ( sum(QinEachProcs) /= Nqs ) then

                write(*, 35)
                35 FORMAT( 'Something wrong in DistributeQpoints ')
                ERROR STOP

            end if Check

        end if OnlyImage1

        SYNC ALL
        call co_broadcast( QinEachProcs, source_image=1, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_broadcast : ', A128 )") msg
        SYNC ALL

        my_Qsize = QinEachProcs(this_image())

        my_edge(1) = sum( QinEachProcs( 1 : (this_image() - 1) ) ) + 1
        my_edge(2) = my_edge(1) + my_Qsize - 1

        my_offset = my_edge(1) - 1

        if ( present(QinWhichProcs) ) then

            QinWhichProcs(1, my_edge(1) : my_edge(2)) = this_image()
            QinWhichProcs(2, my_edge(1) : my_edge(2)) = [ ( i, i=(my_edge(1)-my_offset),(my_edge(2)-my_offset) ) ]

            SYNC ALL
            call co_sum( QinWhichProcs,  stat=istat, ERRMSG=msg )
            if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg
            SYNC ALL

        end if

        if ( present(QPerProcsMax) ) QPerProcsMax = maxval( QinEachProcs )

        deallocate( QinEachProcs )

    end subroutine DistributeQpoints


    subroutine ShowProcessorDistribution( character_to_print, &
                                        & my_edge, my_Qsize, my_offset )

        implicit none

        character(len=*)                                :: character_to_print
        integer, dimension(2), intent(in)               :: my_edge
        integer, intent(in)                             :: my_Qsize, my_offset

        ! =============================== Local variables =============================== !

        integer                             :: char_len, left, right

        character(len=8)                    :: frmt
        character(len=8)                    :: int1_to_char, int2_to_char
        character(len=256)                  :: write_frmt

        ! =============================== Local variables =============================== !

        char_len = LEN_TRIM( character_to_print )

        left = (111 - char_len - 2) / 2
        right = 111 - (left + char_len + 2)

        frmt = '(I5)'
        
        write(int1_to_char, frmt) left
        write(int2_to_char, frmt) right

        write_frmt = "(10X, '|', "//trim(adjustl(adjustr(int1_to_char)))//"('='), ' ', '" &
                   & //trim(adjustl(adjustr(character_to_print)))//"', ' ', "&
                   & //trim(adjustl(adjustr(int2_to_char)))//"('='), '|')" 

        ! ////////////////////  Show the qpoints distributions among processors \\\\\\\\\\\\\\\\\\\\\ !
        img1_chk1: if ( this_image() == 1 ) then
            write(*, *)
            call execute_command_line(' ')

            write(*, write_frmt)
            call execute_command_line(' ')

            write(*, 25)
            call execute_command_line(' ')

        end if img1_chk1

        SYNC ALL

        write(*, 30) this_image(), my_edge, my_Qsize, my_offset
        call execute_command_line(' ')

        SYNC ALL

        img1_chk2: if ( this_image() == 1 ) then

            write(*, 25)
            call execute_command_line(' ')

            write(*, write_frmt)
            call execute_command_line(' ')
            write(*, *)
            call execute_command_line(' ')

        end if img1_chk2

        25 FORMAT(10X,'|',111X,'|')
        30 FORMAT(10X,'| **', 1X,' Image No. = ', I4, ', Q-points covered: (', I9, '=> ', I9, &
                & '). my_Qsize = ', I9, ', my_offset = ', I7, 1X,'** |')
        ! \\\\\\\\\\\\\\\\\\\\  Show the qpoints distributions among processors ///////////////////// !

    end subroutine ShowProcessorDistribution


    subroutine unique_sort( input_arr, output_arr )

        implicit none

        integer, dimension(:), intent(in)                   :: input_arr
        integer, dimension(:), allocatable, intent(out)     :: output_arr

        ! ============================ Local Variables ============================ !

        integer                                             :: num_count, arr_size, &
                                                             & min_val, max_val
        integer, dimension(:), allocatable                  :: unique

        ! ============================ Local Variables ============================ !

        arr_size = size( input_arr )

        allocate( unique(arr_size) )

        num_count = 0
        unique(:) = 0
    
        min_val = minval(input_arr)-1
        max_val = maxval(input_arr)

        do while (min_val < max_val)

            num_count = num_count + 1
            min_val = minval(input_arr, mask=input_arr>min_val)
            unique(num_count) = min_val

        end do

        allocate( output_arr(num_count) )

        output_arr = unique(1:num_count)
        deallocate( unique )

        ! ** Yes!, this is taken (copied) from stackoverflow ** !
        ! https://stackoverflow.com/questions/44198212/a-fortran-equivalent-to-unique !

    end subroutine unique_sort


    subroutine FindPositionAndMultiplicity( num_grid_pnt, num_irr_q, map, unique_irrq, &
                                          & pos, multiplicity )

        implicit none

        integer, intent(in)                                 :: num_grid_pnt, num_irr_q
        integer, dimension(num_grid_pnt), intent(in)        :: map
        integer, dimension(num_irr_q), intent(in)           :: unique_irrq

        integer, dimension(num_grid_pnt), intent(out)       :: pos
        integer, dimension(num_irr_q), intent(out)          :: multiplicity

        ! ============================== Local Variables ============================== !

        integer                                             :: nq, q_grid_indx
        integer, dimension(1)                               :: position_indx

        ! ============================== Local Variables ============================== !

        pos(:) = 0
        multiplicity(:) = 0

        AllQloop: do nq = 1, num_grid_pnt

            q_grid_indx = map(nq)
            position_indx = findloc( unique_irrq, q_grid_indx )

            pos(nq) = position_indx(1)
            multiplicity( position_indx(1) ) = multiplicity( position_indx(1) ) + 1

        end do AllQloop

        Check: if ( sum(multiplicity) /= num_grid_pnt ) then

            write(*, 80) sum(multiplicity), num_grid_pnt
            80 FORMAT( "ERROR in FindPositionAndMultiplicity: sum(multiplicity) /= num_grid_pnt. ", I7, " /= ", I7 )
            ERROR STOP

        end if Check

    end subroutine FindPositionAndMultiplicity


    subroutine WriteAsciiOutFile(filename, write_frmt, freq, LW_q0, Ndof, HeaderWrite)

        implicit none

        character(len=*), intent(in)                :: filename, write_frmt
        real(dp), dimension(:), intent(in)          :: freq, LW_q0
        integer, intent(in)                         :: Ndof
        logical, intent(inout)                      :: HeaderWrite

        ! =================================== Local Variables =================================== !

        integer                                     :: ui, err
        character(len=512)                          :: err_msg

        ! =================================== Local Variables =================================== !

        ui = 5

        FirstWrite: if ( HeaderWrite ) Then

            open( unit=ui, file=filename, status='REPLACE', action='WRITE', &
                & iostat=err, iomsg=err_msg )

                open_chk1: if ( err /= 0 ) then
                    write(*, *) 'File OPEN failed: iostat = ', err
                    write(*, *) 'Error message = ', err_msg

                else open_chk1

                    write(unit=ui, fmt=20) Ndof, Ndof
                    write(unit=ui, fmt=write_frmt) freq, LW_q0

                    HeaderWrite = .false.

                    20 FORMAT( '#--------- phonon_frequency(rad.THz) ( ', I4, ' columns) ---------      ', &
                             & '---------- linewidth(rad.THz) ( ', I4, ' columns) ----------')

                end if open_chk1

            close(unit=ui, status='KEEP', iostat=err, iomsg=err_msg)

            close_chk1: if ( err /= 0 ) then
                 write(*, *) 'File close failed: iostat = ', err
                 write(*, *) 'Error message = ', err_msg
            end if close_chk1

        else FirstWrite

            open( unit=ui, file=filename, status='OLD', action='WRITE', POSITION='APPEND', &
                & iostat=err, iomsg=err_msg )

                open_chk2: if ( err /= 0 ) then
                    write(*, *) 'File OPEN failed: iostat = ', err
                    write(*, *) 'Error message = ', err_msg

                else open_chk2

                    write(unit=ui, fmt=write_frmt) freq, LW_q0

                end if open_chk2

            close(unit=ui, status='KEEP', iostat=err, iomsg=err_msg)

            close_chk2: if ( err /= 0 ) then
                 write(*, *) 'File close failed: iostat = ', err
                 write(*, *) 'Error message = ', err_msg
            end if close_chk2

        end if FirstWrite

    end subroutine WriteAsciiOutFile


    Function kron_prod(invec1, invec2) Result(outvec)

        implicit none

        integer, parameter                                  :: len_v1=3, len_v2=3

        real(dp), dimension(3), intent(in)                  :: invec1, invec2
        real(dp), dimension(9)                              :: outvec !Result

        ! ================================== Local Variables ================================== !

        integer                                             :: vi, strtp, endp

        ! ================================== Local Variables ================================== !

        strtp = 1
        endp = len_v2

        invec1_loop: do vi = 1, len_v1

            outvec(strtp : endp) = invec1(vi) * invec2(:)

            strtp = endp + 1
            endp = endp + len_v2

        end do invec1_loop

    end Function kron_prod

end module Helper

