
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


subroutine ReducePntPerASR(this, Nbasis, ZeroEPS, Cut_Off)

    implicit none

    class(FCtyp)                            :: this

    integer, intent(in)                     :: Nbasis
    real(dp), intent(in)                    :: ZeroEPS
    type(CutOff), intent(in)                :: Cut_Off

    ! =========================================== Local Variables =========================================== !

    real(dp), dimension(:,:), allocatable   :: InterFC

    integer                                 :: mu, Natm, NumCol
    integer                                 :: num, nf, nel, ind, nf2, istat

    integer                                 :: numfc
    integer, dimension(2*ordfc3)            :: FCindx

    real(dp)                                :: pivot, cf, pivot2, indx_chk
    real(dp), dimension(:), allocatable     :: after, row, coef_lin, row2, zero_chk

    integer, dimension(:,:), allocatable    :: PosIndxFlat_tmp

    character(len=256)                      :: msg

    ! =========================================== Local Variables =========================================== !

    write(*, 40)

    numfc = this%ind_f
    ! --------------------------------------- Initialize & Allocation --------------------------------------- !

    !-! allocate( this%RedInFC(numfc, numfc), STAT=istat, ERRMSG=msg )
    !-! if ( istat /= 0 ) call this%AllocationError( istat, msg )
    !-! this%RedInFC = 0.0_dp

    this%Nzrofin = 0 
    this%Depfin = 0 
    this%Indfin = 0 
    this%Ignorefin = 0

    allocate( this%PosIndx_fin(Nbasis) )
    mu_loop: do mu = 1, Nbasis

        Natm = Cut_Off%numAtom(mu)
        allocate( this%PosIndx_fin(mu)%Pos_mu(3,3,3,Natm,Natm), STAT=istat, ERRMSG=msg )
        if ( istat /= 0 ) call this%AllocationError( istat, msg )

        this%PosIndx_fin(mu)%Pos_mu = FILLVAL

    end do mu_loop

    allocate( this%PosIndx_flat_fin( (2*ordfc3), numfc ) )
    this%PosIndx_flat_fin = 0

    this%PosIndCntr_fin = 0

    ! --------------------------------------- Initialize & Allocation --------------------------------------- !

    call this%ReadChunks(InterFC)

    write(*, *)
    num = 0
    fc_loop: do nf = numfc, 1, -1

        num = num + 1

        !(/mu_ind, N2_ind, N3_ind, alpha, beta, gama/)
        FCindx = this%PosIndx_flat(:, nf)

        pivot = InterFC(nf, nf)

        nel = numfc - nf
        allocate( after(nel) )
        after = InterFC((nf+1):numfc, nf)

        PivotChk: if ( dabs(pivot) > ZeroEPS ) then

            AftrChk: if ( (size(after) == 0) .or. (all(dabs(after) < ZeroEPS)) ) then

                write(*, 50) 
                write(*, 80) FCindx

                this%F(FCindx(1))%FCmu(1, FCindx(6), FCindx(5), &
                                        & FCindx(4), FCindx(3), FCindx(2)) = 8.0_dp
                InterFC(nf, nf) = 8.0_dp
                this%Nzrofin = this%Nzrofin + 1

            else AftrChk !if ( (size(after) /= 0) .and. (any(dabs(after) > ZeroEPS)) ) then AftrChk

                allocate( row(nel) )
                allocate( coef_lin(nel) )

                row = 0.0_dp
                coef_lin = after / ( -1.0_dp * pivot )

                CoefLinloop: do ind = 1, nel

                    cf = coef_lin(ind)
                    NonZeroCf: if ( dabs(cf) > ZeroEPS ) then

                        nf2 = nf + ind
                        pivot2 = InterFC(nf2, nf2) !*!

                        Pivot2Chk: if ( dabs(pivot2) < EPS ) then

                            row(ind) = row(ind) + cf

                        else if ( dabs(pivot2 - 1.0_dp) < EPS ) then Pivot2Chk

                            allocate( row2(nel) )
                            row2 = 0.0_dp
                            !strt = ind + 1 !nf2 - nf + 1
                            row2( (ind+1):nel ) = InterFC( (nf2+1):numfc, nf2 ) !*!

                            row = row + (cf * row2)

                            !-! Daxpy: if ( nel < 1024 ) then
                            !-!     row = row + (cf * row2)
                            !-! else Daxpy
                            !-!     call daxpy(nel, cf, row2, 1, row, 1)
                            !-! end if Daxpy

                            deallocate( row2 )

                        else Pivot2Chk

                            not8: if ( dabs(pivot2-8.0_dp) > EPS ) then
                                write(*, 100) pivot2
                                STOP
                            end if not8

                        end if Pivot2Chk

                    end if NonZeroCf

                end do CoefLinloop

                InterFC(nf, nf) = 1.0_dp !*!
                InterFC( (nf+1):numfc, nf ) = row !*!

                this%Depfin = this%Depfin + 1
                this%F(FCindx(1))%FCmu(1, FCindx(6), FCindx(5), &
                                        & FCindx(4), FCindx(3), FCindx(2)) = -9.0_dp

                deallocate( row, coef_lin )

            end if AftrChk

        else if (dabs(pivot) < ZeroEPS) then PivotChk

            this%Indfin = this%Indfin + 1
            call this%MakeIndIndxFin(FCindx)

            DebugIgnoreChk: if ( (size( after ) /= 0) .and. (any(dabs(after) > ZeroEPS)) ) then

                this%Ignorefin = this%Ignorefin + 1

            end if DebugIgnoreChk

            !-! InterFC(nf, nf) = 0.0_dp
            InterFC( nf:numfc, nf ) = 0.0_dp
            !-!

        else PivotChk
            write(*, 120) pivot
            STOP

        end if PivotChk

        allocate( zero_chk(nel) )
        zero_chk = InterFC( (nf+1):numfc, nf )
        indx_chk = InterFC(nf, nf)

        ZeroReduce: if ( (all(dabs(zero_chk) < ZeroEPS)) .and. (dabs(indx_chk-1.0_dp)<EPS) ) then

            write(*, 150)
            write(*, 80) FCindx

            this%F(FCindx(1))%FCmu(1, FCindx(6), FCindx(5), &
                                    & FCindx(4), FCindx(3), FCindx(2)) = 8.0_dp
            InterFC(nf, nf) = 8.0_dp !*!
            this%Nzrofin = this%Nzrofin + 1

        end if ZeroReduce

        deallocate( zero_chk )

        deallocate( after )

        NonAdvancingOut: if ( this%NonAdvancWrite ) then
            write(6, 500, advance="no") char(13), num, numfc, this%Indfin, this%Depfin, this%Nzrofin
            FLUSH(6)

        else NonAdvancingOut
            write(*, 550) num, numfc, this%Indfin, this%Depfin, this%Nzrofin

        end if NonAdvancingOut

    end do fc_loop

    write(*, *)

    ResizePosIndxFlat: if ( this%PosIndCntr_fin == this%Indfin ) then

        NumCol = this%Indfin

        allocate( PosIndxFlat_tmp((2*ordfc3), NumCol) )
        PosIndxFlat_tmp(:, :) = this%PosIndx_flat_fin(:, 1:NumCol)
        deallocate( this%PosIndx_flat_fin )

        !call move_alloc(PosIndxFlat_tmp, this%PosIndx_flat_fin)
        allocate( this%PosIndx_flat_fin((2*ordfc3), NumCol) )
        this%PosIndx_flat_fin = PosIndxFlat_tmp
        deallocate( PosIndxFlat_tmp )

    else ResizePosIndxFlat

        write(*, 60) this%PosIndCntr_fin, this%Indfin
        STOP

    end if ResizePosIndxFlat
    !~!
    call move_alloc( InterFC, this%RedInFC )
    !~!

    40 FORMAT("Reducing from Inter-FC dependent matrix (Found from Point Group, Permutation and ASR)")
    50 FORMAT("Zero value of FC found ( Potential for Error )")
    60 FORMAT("ERROR( in ReducePntPerASR ): this%PosIndCntr_fin /= this%Indfin ", I5, I5)
    80 FORMAT("The FC corresponfing to the indx: [", 5(I3, ', '), I3, "] is zero.")
    100 FORMAT("ERROR( in ReducePntPerASR ): Illegal value of pivot2 ( can not be other than 0, 1, 8) ", F6.3)
    120 FORMAT("ERROR( in ReducePntPerASR ): Illegal value of pivot ", F6.3)
    150 FORMAT("Zero value of FC found from reduction ( Potential for Error )")

    500 FORMAT(1a1, "FC Covered: ", I5, " /", I5, " /*\ ind_fc= ", I5, ", dep_fc= ", I5, ", zero_fc= ", I5)
    550 FORMAT("FC Covered: ", I5, " /", I5, " /*\ ind_fc= ", I5, ", dep_fc= ", I5, ", zero_fc= ", I5)

end subroutine ReducePntPerASR


subroutine MakeIndIndxFin(this, FCindx)

    implicit none

    class(FCtyp)                                :: this

    integer, dimension(2*ordfc3), intent(in)    :: FCindx


    ! ====================================== Local Variables ====================================== !

    integer                                     :: val

    ! ====================================== Local Variables ====================================== !

    val = this%PosIndx_fin(FCindx(1))%Pos_mu(FCindx(6), FCindx(5), &
                                           & FCindx(4), FCindx(3), FCindx(2)) 

    ValChk: if ( val /= FILLVAL ) then
        write(*, 180) val
        STOP

    else ValChk

        this%PosIndCntr_fin = this%PosIndCntr_fin + 1

        ! ...........::::::::::::: For flat writing :::::::::::::........... !
        SizeChk: if ( this%PosIndCntr_fin <= this%ind_f ) then

            this%PosIndx_flat_fin(:, this%PosIndCntr_fin) = FCindx

        else SizeChk

            write(*, 200) this%PosIndCntr_fin, this%ind_f
            STOP

        end if SizeChk
        ! ...........::::::::::::: For flat writing :::::::::::::........... !

        this%PosIndx_fin(FCindx(1))%Pos_mu(FCindx(6), FCindx(5), &
                                         & FCindx(4), FCindx(3), FCindx(2)) = this%PosIndCntr_fin 

    end if ValChk

    180 FORMAT("ERROR( in MakeIndIndxFin ): Initially this%PosIndx_fin must be FILLVAL ", I6)
    200 FORMAT("ERROR( in MakeIndIndxFin ): PosIndCntr_fin can not exceed ind_f ", I6, I6)

end subroutine MakeIndIndxFin

