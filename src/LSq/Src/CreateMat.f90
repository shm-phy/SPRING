
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


module CreateMat

    use kinds,                  only : dp
    use constants,              only : EPS_nzero
    use unit_cell,              only : cell
    use hdf5_wrap,              only : ReadDispMat, ReadForcedd, ReadForceRenorm, &
                                     & ReadFlags

    implicit none
    private

    public      :: MakeLsqData, MakeLsqDataRenorm, LSqDataPointDef

contains

    subroutine MakeLsqData(sys, T, Dmat, Frow, MatShp, OrdArr)

        implicit none

        type(cell), intent(in)                              :: sys
        real(dp), intent(in)                                :: T
        real(dp), dimension(:,:), allocatable, intent(out)  :: Dmat
        real(dp), dimension(:), allocatable, intent(out)    :: Frow
        integer, dimension(2,3), intent(out)                :: MatShp
        logical, dimension(3), intent(out)                  :: OrdArr

        ! ================================== Local Variables ================================== !

        character(len=24)                                   :: Temp_char
        character(len=8)                                    :: frmt

        real(dp), dimension(:,:), allocatable               :: Mat2, Mat3, Mat4
        real(dp), dimension(:), allocatable                 :: Frow2, Frow3, Frow4, &
                                                             & Frow_dd, Frow_md

        integer, dimension(2)                               :: MatShape2, MatShape3, MatShape4
        character(len=128)                                  :: FD2file, FD3file, FD4file, F_dd_file

        integer                                             :: NRow, NCol, strtCol, endCol, ord

        ! ================================== Local Variables ================================== !

        write(*, *)

        MatShp(:, :) = 0
        OrdArr(:) = .false.

        frmt = '(F6.1)'
        write(Temp_char, frmt) T

        OuterIf: if ( sys%SecondFC ) then

            OrdArr(1) = .true.
            FD2file = 'FD_2nd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

            call ReadDispMat(FD2file, Mat2, Frow2, MatShape2)
            MatShp(:, 1) = MatShape2

            InerIf1: if( sys%ThirdFC ) then
                
                OrdArr(2) = .true.
                FD3file = 'FD_3rd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

                call ReadDispMat(FD3file, Mat3, Frow3, MatShape3)
                MatShp(:, 2) = MatShape3

            end if InerIf1

            InerIf2: if( sys%FourthFC ) then

                OrdArr(3) = .true.
                FD4file = 'FD_4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

                call ReadDispMat(FD4file, Mat4, Frow4, MatShape4)
                MatShp(:, 3) = MatShape4

            end if InerIf2

            InerIf3: if ( sys%Force_dd ) then

                F_dd_file = 'FC2_dd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

                call ReadForcedd(F_dd_file, Frow_dd, Frow_md)

            end if InerIf3

        else OuterIf

            write(*, 125)
            STOP

        end if OuterIf

        NRow = MatShp(1,1)
        NCol = sum( MatShp(2, :) )

        allocate( Dmat(NRow, NCol) )
        Dmat = 0.0_dp

        allocate( Frow(NRow) )
        Frow = Frow2

        deallocate( Frow2 )

        OrdLoop: do ord = 1, 3

            OutIf: if ( OrdArr(ord) ) then

                Inif: if ( ord == 1 ) then

                    strtCol = 1
                    endCol = MatShp(2, 1)
                    Dmat(:, strtCol:endCol) = Mat2

                    deallocate( Mat2 )

                else if ( ord == 2 ) then Inif

                    strtCol = MatShp(2, 1) + 1
                    endCol = strtCol + MatShp(2, 2) - 1
                    Dmat(:, strtCol:endCol) = Mat3

                    deallocate( Mat3 )

                    debug1: if ( any(dabs(Frow-Frow3) > EPS_nzero) ) then
                        write(*, 150) 
                        STOP
                    end if debug1

                    deallocate( Frow3 )

                else Inif

                    strtCol = MatShp(2, 1) + MatShp(2,2) + 1
                    endCol = strtCol + MatShp(2, 3) - 1
                    Dmat(:, strtCol:endCol) = Mat4

                    deallocate( Mat4 )

                    debug2: if ( any(dabs(Frow-Frow4) > EPS_nzero) ) then
                        write(*, 175) 
                        STOP
                    end if debug2

                    deallocate( Frow4 )

                end if Inif

            end if OutIf

        end do OrdLoop

        FrcddChk: if ( sys%Force_dd ) then

            debug3: if ( any(dabs(Frow-Frow_md) > EPS_nzero) ) then
                write(*, 200) 
                STOP
            end if debug3

            deallocate( Frow_md )

            Frow = Frow - Frow_dd

            deallocate( Frow_dd )

        end if FrcddChk

        125 FORMAT("Can not proceed without 2nd Order FC. SecondFC is set to .false.")
        150 FORMAT("ERROR(in MakeLsqData) : Frow2 and Frow3 are not same")
        175 FORMAT("ERROR(in MakeLsqData) : Frow2 and Frow4 are not same")
        200 FORMAT("ERROR(in MakeLsqData) : Frow2 and Frow_md are not same")

    end subroutine MakeLsqData


    subroutine MakeLsqDataRenorm(sys, T, Dmat, Frow, MatShp, OrdArr)

        implicit none

        type(cell), intent(in)                              :: sys
        real(dp), intent(in)                                :: T
        real(dp), dimension(:,:), allocatable, intent(out)  :: Dmat
        real(dp), dimension(:), allocatable, intent(out)    :: Frow
        integer, dimension(2,3), intent(out)                :: MatShp
        logical, dimension(3), intent(out)                  :: OrdArr

        ! ================================== Local Variables ================================== !

        character(len=24)                                   :: Temp_char
        character(len=8)                                    :: frmt

        real(dp), dimension(:,:), allocatable               :: Mat3, Mat4
        real(dp), dimension(:), allocatable                 :: Force2Renorm, Frow3, Frow4, &
                                                             & Frow_dd, Frow_md1, Frow_md2

        integer, dimension(2)                               :: MatShape3, MatShape4
        character(len=128)                                  :: Force2RenormFile, FD3file, FD4file, & 
                                                             & F_dd_file

        integer                                             :: NRow, NCol, strtCol, endCol, ord

        ! ================================== Local Variables ================================== !

        write(*, *)

        MatShp(:, :) = 0
        OrdArr(:) = .false.

        frmt = '(F6.1)'
        write(Temp_char, frmt) T

        ! ------------------------------ Read Force from the renormalized 2nd-IFCs ------------------------------ !
        Force2RenormFile = 'FC2ndRenorm_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K_F.h5'

        call ReadForceRenorm(Force2RenormFile, Force2Renorm, Frow_md1)
        ! ------------------------------ Read Force from the renormalized 2nd-IFCs ------------------------------ !

        InerIf1: if( sys%ThirdFC ) then
            
            OrdArr(2) = .true.
            FD3file = 'FD_3rd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

            call ReadDispMat(FD3file, Mat3, Frow3, MatShape3)
            MatShp(:, 2) = MatShape3

        else InerIf1

            write(*, 13)
            13 FORMAT('WARNING: ThirdFC is set false when calculating renormalized anharmonic IFCs')

        end if InerIf1

        InerIf2: if( sys%FourthFC ) then

            OrdArr(3) = .true.
            FD4file = 'FD_4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

            call ReadDispMat(FD4file, Mat4, Frow4, MatShape4)
            MatShp(:, 3) = MatShape4

        else InerIf2

            write(*, 23)
            23 FORMAT('WARNING: FourthFC is set false when calculating renormalized anharmonic IFCs')

        end if InerIf2

        InerIf3: if ( sys%Force_dd ) then

            F_dd_file = 'FC2_dd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

            call ReadForcedd(F_dd_file, Frow_dd, Frow_md2)

        end if InerIf3

        debug4: if ( MatShp(1,2) /= MatShp(1,3) ) then

            write(*, 33)
            33 FORMAT('WARNING: Number of rows in Mat3 and Mat4 are not same')

        end if debug4

        NRow = MatShp(1,2)
        NCol = sum( MatShp(2, 2:) )

        allocate( Dmat(NRow, NCol) )
        Dmat = 0.0_dp

        allocate( Frow(NRow) )
        Frow = Frow_md1

        deallocate( Frow_md1 )

        OrdLoop: do ord = 2, 3

            OutIf: if ( OrdArr(ord) ) then

                Inif: if ( ord == 2 ) then 

                    strtCol = 1
                    endCol = MatShp(2, 2) 
                    Dmat(:, strtCol:endCol) = Mat3

                    deallocate( Mat3 )

                    debug1: if ( any(dabs(Frow-Frow3) > EPS_nzero) ) then
                        write(*, 150) 
                        STOP
                    end if debug1

                    deallocate( Frow3 )

                else Inif

                    strtCol = MatShp(2,2) + 1
                    endCol = strtCol + MatShp(2, 3) - 1
                    Dmat(:, strtCol:endCol) = Mat4

                    deallocate( Mat4 )

                    debug2: if ( any(dabs(Frow-Frow4) > EPS_nzero) ) then
                        write(*, 175) 
                        STOP
                    end if debug2

                    deallocate( Frow4 )

                end if Inif

            end if OutIf

        end do OrdLoop

        FrcddChk: if ( sys%Force_dd ) then

            debug3: if ( any(dabs(Frow-Frow_md2) > EPS_nzero) ) then
                write(*, 200) 
                STOP
            end if debug3

            deallocate( Frow_md2 )

            Frow = Frow - Force2Renorm - Frow_dd

            deallocate( Frow_dd )

        else FrcddChk

            Frow = Frow - Force2Renorm

        end if FrcddChk

        150 FORMAT("ERROR(in MakeLsqDataRenorm) : Frow_md1 and Frow3 are not same")
        175 FORMAT("ERROR(in MakeLsqDataRenorm) : Frow_md1 and Frow4 are not same")
        200 FORMAT("ERROR(in MakeLsqDataRenorm) : Frow_md1 and Frow_md2 are not same")

    end subroutine MakeLsqDataRenorm


    subroutine LSqDataPointDef( sys, T, Dmat, Frow )

        implicit none

        type(cell), intent(in)                                  :: sys
        real(dp), intent(in)                                    :: T
        real(dp), dimension(:,:), allocatable, intent(inout)    :: Dmat
        real(dp), dimension(:), allocatable, intent(inout)      :: Frow

        ! ================================== Local Variables ================================== !

        character(len=24)                                   :: Temp_char
        character(len=8)                                    :: frmt
        character(len=128)                                  :: FD2file, FD3file, FD4file

        real(dp), dimension(:,:), allocatable               :: Dmat_tmp, Dmat_tmp2
        real(dp), dimension(:), allocatable                 :: Frow_tmp

        integer, dimension(:), allocatable                  :: flags, flags_3rd, flags_4th
        integer                                             :: NumIFC, NumForces, ii, forces_count

        ! ================================== Local Variables ================================== !

        write(*, *)

        frmt = '(F6.1)'
        write(Temp_char, frmt) T

        OuterIf: if ( sys%SecondFC ) then

            FD2file = 'FD_2nd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

            call ReadFlags( FD2file, flags )

            InerIf1: if( sys%ThirdFC ) then
                
                FD3file = 'FD_3rd_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

                call ReadFlags( FD3file, flags_3rd )
                flags = flags + flags_3rd

                deallocate( flags_3rd )

            end if InerIf1

            InerIf2: if( sys%FourthFC ) then

                FD4file = 'FD_4th_'//trim(sys%prefix)//trim(adjustl(adjustr(Temp_char)))//'K.h5'

                call ReadFlags( FD4file, flags_4th )
                flags = flags + flags_4th

                deallocate( flags_4th )

            end if InerIf2

        else OuterIf

            write(*, 125)
            STOP

        end if OuterIf

        Allocation: if ( (.not. ALLOCATED(Dmat)) .or. (.not. ALLOCATED(Frow)) ) then

            write(*, 130)
            STOP
            
        else Allocation

            forces_count = 0

            NumIFC = size(Dmat, 2)
            NumForces = size(Dmat, 1)

            allocate( Dmat_tmp( NumIFC, NumForces) )
            Dmat_tmp = TRANSPOSE( Dmat )
            allocate( Dmat_tmp2( NumIFC, NumForces) )
            Dmat_tmp2 = 0.0_dp

            if ( size(Frow) /= NumForces ) then
                write(*, 135)
                STOP
            else
                allocate( Frow_tmp(NumForces) )
                Frow_tmp = 0.0_dp
            end if

            ValidForcesLoop: do ii = 1, NumForces

                check: if ( flags(ii) == 0 ) then

                    forces_count = forces_count + 1

                    Dmat_tmp2(:, forces_count) = Dmat_tmp(:, ii)
                    Frow_tmp(forces_count) = Frow( ii )

                end if check

            end do ValidForcesLoop

            !call move_alloc( Dmat_tmp2(:, 1:forces_count), Dmat_tmp )
            deallocate( Dmat_tmp )
            allocate( Dmat_tmp( NumIFC, forces_count ) )
            Dmat_tmp(:, :) = Dmat_tmp2(:, 1:forces_count)
            deallocate( Dmat_tmp2 )
            
            !call move_alloc( TRANSPOSE(Dmat_tmp), Dmat )
            deallocate( Dmat )
            allocate( Dmat(forces_count, NumIFC) )
            Dmat = TRANSPOSE( Dmat_tmp )
            deallocate( Dmat_tmp )

            !call move_alloc( Frow_tmp(1:forces_count), Frow )
            deallocate( Frow )
            allocate( Frow( forces_count ) )
            Frow(:) = Frow_tmp(1:forces_count)
            deallocate( Frow_tmp )

            write(*, 140) forces_count, NumIFC

        end if Allocation

        deallocate( flags )

        125 FORMAT("Can not proceed without 2nd Order FC. SecondFC is set to .false.")
        130 FORMAT( "Either 'Dmat' or 'Frow' or both are not allocated on entry" )
        135 FORMAT( "Something wrong: size(Frow) /= NumForces" )
        140 FORMAT( "Point-Defect => Dimension of the matrix for least-square : (", I6, "X ", I5, " )" )

    end subroutine LSqDataPointDef

end module CreateMat

