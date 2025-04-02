
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
    use CreateMat,                  only : MakeLsqData, MakeLsqDataRenorm, LSqDataPointDef
    use timer_class,                only : timer
    use mklWrap,                    only : RankDefLSqr, RankDefLSqr2
    use hdf5_wrap,                  only : WriteIndFC

    implicit none

    character(len=256)                          :: filename
    real(dp)                                    :: T, SingValCut
    logical                                     :: Renorm, &
                                                 & PointDef !** Variable for point-defect calculation **!

    type(cell)                                  :: sys

    real(dp), dimension(:,:), allocatable       :: Dmat
    real(dp), dimension(:), allocatable         :: Frow
    integer, dimension(2,3)                     :: MatShp
    logical, dimension(3)                       :: OrdArr

    type(timer)                                 :: TElsp
    real(dp), dimension(:), allocatable         :: Xout

    integer                                     :: ord, strtCol, endCol
    real(dp), dimension(:), allocatable         :: ind_fc
    character(len=128)                          :: FDfilename, grp_name
    character(len=24)                           :: Temp_char
    character(len=8)                            :: frmt


    ! ------------------------------------------------------------------------------------- !

    call ShowWelcomeBanner()

    call get_command_line(filename, T, SingValCut, Renorm, PointDef)

    frmt = '(F6.1)'
    write(Temp_char, frmt) T

    call sys%init_data(filename)

    RenormChk: if ( .not. Renorm ) then

        call MakeLsqData(sys, T, Dmat, Frow, MatShp, OrdArr)
        if ( PointDef ) call LSqDataPointDef( sys, T, Dmat, Frow )

        call TElsp%start_timer()

        call RankDefLSqr( Dmat, Frow, Xout, SingValCut )
        deallocate( Dmat )
        deallocate( Frow )

        write(*, 80) ( TElsp%elapsed_time() / 60.0_dp )

        grp_name = 'Ind_FCs'

        OrdLoop: do ord = 1, 3

            OutIf: if ( OrdArr(ord) ) then

                InIf: if ( ord == 1 ) then

                    FDfilename = 'FD_2nd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

                    strtCol = 1
                    endCol = MatShp(2, 1)
                    allocate( ind_fc(MatShp(2, 1)) )
                    ind_fc = Xout(strtCol:endCol)

                    call WriteIndFC(FDfilename, grp_name, ind_fc)

                    deallocate( ind_fc )

                else if ( ord == 2 ) then InIf

                    FDfilename = 'FD_3rd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

                    strtCol = MatShp(2, 1) + 1
                    endCol = strtCol + MatShp(2, 2) - 1
                    allocate( ind_fc(MatShp(2, 2)) )
                    ind_fc = Xout(strtCol:endCol)

                    call WriteIndFC(FDfilename, grp_name, ind_fc)

                    deallocate( ind_fc )

                else InIf

                    FDfilename = 'FD_4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

                    strtCol = MatShp(2, 1) + MatShp(2,2) + 1
                    endCol = strtCol + MatShp(2, 3) - 1
                    allocate( ind_fc(MatShp(2, 3)) )
                    ind_fc = Xout(strtCol:endCol)

                    call WriteIndFC(FDfilename, grp_name, ind_fc)

                    deallocate( ind_fc )

                end if InIf

            end if OutIf

        end do OrdLoop
        
    else RenormChk

        call MakeLsqDataRenorm(sys, T, Dmat, Frow, MatShp, OrdArr)
        if ( PointDef ) call LSqDataPointDef( sys, T, Dmat, Frow )

        call TElsp%start_timer()

        call RankDefLSqr( Dmat, Frow, Xout, SingValCut )
        deallocate( Dmat )
        deallocate( Frow )

        write(*, 80) ( TElsp%elapsed_time() / 60.0_dp )

        grp_name = 'Ind_FCs'

        OrdLoop2: do ord = 2, 3

            OutIf2: if ( OrdArr(ord) ) then

                InIf2: if ( ord == 2 ) then 

                    FDfilename = 'FD_3rd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

                    strtCol = 1
                    endCol = MatShp(2, 2) 
                    allocate( ind_fc(MatShp(2, 2)) )
                    ind_fc = Xout(strtCol:endCol)

                    call WriteIndFC(FDfilename, grp_name, ind_fc)

                    deallocate( ind_fc )

                else InIf2

                    FDfilename = 'FD_4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

                    strtCol = MatShp(2,2) + 1
                    endCol = strtCol + MatShp(2, 3) - 1
                    allocate( ind_fc(MatShp(2, 3)) )
                    ind_fc = Xout(strtCol:endCol)

                    call WriteIndFC(FDfilename, grp_name, ind_fc)

                    deallocate( ind_fc )

                end if InIf2

            end if OutIf2

        end do OrdLoop2

    end if RenormChk

    write(*, '(/)')

    80 FORMAT("Elapsed time for Linear least squares routine: ", F7.4, " min.")

end program main

