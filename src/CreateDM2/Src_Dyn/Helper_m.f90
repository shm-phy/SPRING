
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
    PRIVATE
    public  :: ReadPointDefVar

contains

    subroutine ReadPointDefVar( filename, Nbasis, sys, qMesh_pd, Mass_def, PointDef, NumDef_sup, &
                              & cell_bsatm )

        implicit none

        real(dp), parameter         :: FILL=99999.0_dp

        character(len=*), intent(in)                        :: filename
        integer, intent(in)                                 :: Nbasis
        type(cell), intent(in)                              :: sys

        integer, dimension(3), intent(out)                  :: qMesh_pd
        real(dp), dimension(Nbasis), intent(out)            :: Mass_def
        logical, intent(out)                                :: PointDef
        integer, intent(out)                                :: NumDef_sup
        integer, dimension(:, :), allocatable, intent(out)  :: cell_bsatm

        ! ================================ Local Variables  ================================ !
        character(len=256)                      :: err_msg
        real(dp), dimension(Nbasis)             :: Mass_def_c
        integer                                 :: err, mu, cnt
        integer, dimension(Nbasis)              :: def_pos
        integer, dimension(4, MAX_DEF)          :: cell_atom
        ! ================================ Local Variables  ================================ !

        namelist   /PointDefInfo/   qMesh_pd, Mass_def, def_pos, NumDef_sup, cell_atom

        !** Default Values **!
        PointDef = .false.
        qMesh_pd = (/7, 7, 7/)
        Mass_def(:) = FILL
        def_pos = 0
        NumDef_sup = 0
        cell_atom(:, :) = int( FILL )
        !** Default Values **!
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
                if ( this_image() == 1 ) write(*, 105) PointDef
                105 FORMAT( "Point defect related calculation is set to: ", L7 )

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

end module Helper

