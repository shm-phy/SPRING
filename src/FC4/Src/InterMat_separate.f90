
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


subroutine Extend_Mat(this)

    implicit none

    class(FCtyp)                                :: this

    ! ================================= Local Variables ================================= !

    real(dp), dimension(:,:), allocatable       :: InterMat_tmp
    integer                                     :: NumCol, NumRow, InterMatCol, InterMatRow

    ! ================================= Local Variables ================================= !

    allocated_chk: if ( .not. (allocated(this%InterMat)) ) then

        NumRow = this%ind_f
        NumCol = (this%NumInterPer + this%NumInterPnt) * (3**ordfc4)

        allocate( this%InterMat(NumRow, NumCol) )
        this%InterMat = 0.0_dp

        this%InterMatPos = 1

    else allocated_chk

        InterMatRow = size( this%InterMat, 1 )
        InterMatCol = size( this%InterMat, 2 )

        !call move_alloc(this%InterMat, InterMat_tmp)
        allocate( InterMat_tmp(InterMatRow, InterMatCol) )
        InterMat_tmp = this%InterMat
        deallocate( this%InterMat )

        NumRow = this%ind_f ! InterMatRow + (this%ind_f - InterMatRow)
        NumCol = InterMatCol + ( (this%NumInterPer + this%NumInterPnt) * (3**ordfc4) )

        allocate( this%InterMat(NumRow, NumCol) )
        this%InterMat = 0.0_dp

        this%InterMat(1:InterMatRow, 1:InterMatCol) = InterMat_tmp
        deallocate( InterMat_tmp )

        this%InterMatPos = InterMatCol + 1

    end if allocated_chk

end subroutine Extend_Mat


subroutine CreateInterGmat(this, mu_dep, dep_N2, dep_N3, dep_N4, dep_var, &
                               & mu_ind, N2_ind, N3_ind, N4_ind, mat_num, PntPer, comb)

    implicit none

    class(FCtyp)                                :: this

    integer, intent(in)                         :: mu_dep, dep_N2, dep_N3, dep_N4, dep_var, &
                                                 & mu_ind, N2_ind, N3_ind, N4_ind, mat_num
    type(PntPer_typ), intent(in)                :: PntPer

    integer, dimension(ordfc4), intent(in)      :: comb

    ! ===================================== Local Variables ===================================== !

    real(dp), dimension(:), allocatable         :: row

    integer                                     :: a, b, c, d
    integer, dimension(ordfc4)                  :: dep_cart

    ! ===================================== Local Variables ===================================== !

    allocate( row(this%ind_f) )
    row = 0.0_dp

    a_loop: do a = 1, 3
        b_loop: do b = 1, 3
            c_loop: do c = 1, 3
                d_loop: do d = 1, 3

                    !row = 0.0_dp
                    dep_var_chk: if ( dep_var == -7 ) then
                        call this%LeftSidePer(mu_ind, N2_ind, N3_ind, N4_ind, &
                                            & a, b, c, d, comb, dep_cart, row)
                    else dep_var_chk
                        call this%LeftSidePntGr(mu_ind, N2_ind, N3_ind, N4_ind, &
                                              & a, b, c, d, mat_num, PntPer, dep_cart, row)
                    end if dep_var_chk

                    call this%RightSide(mu_dep, dep_N2, dep_N3, dep_N4, &
                                      & dep_cart(1), dep_cart(2), dep_cart(3), dep_cart(4), row)

                    !*! debug: if ( any(dabs(this%InterMat(:, this%InterMatPos)) > EPS) ) then
                    !*!     write(*, 20)
                    !*!     STOP
                    !*! end if debug

                    !^! Not necessary !^!
                    WHERE( dabs(row) < EPS )  row = 0.0_dp
                    !^! Not necessary !^!

                    this%InterMat(:, this%InterMatPos) = row
                    this%InterMatPos = this%InterMatPos + 1

                    row = 0.0_dp

                end do d_loop
            end do c_loop
        end do b_loop
    end do a_loop

    deallocate( row )

    !*! 20 FORMAT("ERROR( in CreateInterGmat ): Initially the marix column is not zero")

end subroutine CreateInterGmat


subroutine LeftSidePer(this, mu_ind, N2_ind, N3_ind, N4_ind, &
                     & a, b, c, d, comb, dep_cart, row)

    implicit none
    
    class(FCtyp)                                :: this

    integer, intent(in)                         :: mu_ind, N2_ind, N3_ind, N4_ind, &
                                                 & a, b, c, d
    integer, dimension(ordfc4), intent(in)      :: comb

    integer, dimension(ordfc4), intent(out)     :: dep_cart
    real(dp), dimension(:), intent(inout)       :: row

    ! ===================================== Local Variables ===================================== !

    real(dp), dimension(3**ordfc4)              :: coef
    integer, dimension(ordfc4)                  :: ind_cart
    integer                                     :: af, i_ind

    ! ===================================== Local Variables ===================================== !

    coef = 0.0_dp

    ind_cart = (/a, b, c, d/)

    af_loop: do af = 1, ordfc4

        dep_cart(af) = ind_cart( comb(af) )

    end do af_loop

    i_ind = 27*(ind_cart(1) - 1) + 9*(ind_cart(2) -1) + 3*(ind_cart(3) - 1) + ind_cart(4) !-1 + 1

    coef(i_ind) = coef(i_ind) + 1.0_dp
    coef = -1.0_dp * coef !important -1

    call this%Find_Lin_Dep(mu_ind, N2_ind, N3_ind, N4_ind, coef, row)

end subroutine LeftSidePer


subroutine LeftSidePntGr(this, mu_ind, N2_ind, N3_ind, N4_ind, &
                       & a, b, c, d, mat_num, PntPer, dep_cart, row)

    implicit none

    class(FCtyp)                                :: this

    integer, intent(in)                         :: mu_ind, N2_ind, N3_ind, N4_ind, &
                                                 & a, b, c, d, mat_num
    type(PntPer_typ), intent(in)                :: PntPer

    integer, dimension(ordfc4), intent(out)     :: dep_cart
    real(dp), dimension(:), intent(inout)       :: row

    ! ===================================== Local Variables ===================================== !

    real(dp), dimension(3,3)                    :: matxyz
    real(dp), dimension(3)                      :: vec_a, vec_b, vec_c, vec_d
    real(dp), dimension(:), allocatable         :: outkron9_1, outkron9_2
    real(dp), dimension(:), allocatable         :: outkron81_fnl
    real(dp), dimension(3**ordfc4)              :: coef

    ! ===================================== Local Variables ===================================== !

    matxyz = PntPer%Rxyz(:, :, mat_num)
    vec_a = matxyz(a, :)
    vec_b = matxyz(b, :)
    call kron_prod(vec_a, vec_b, outkron9_1)

    vec_c = matxyz(c, :)
    vec_d = matxyz(d, :)
    call kron_prod(vec_c, vec_d, outkron9_2)

    call kron_prod(outkron9_1, outkron9_2, outkron81_fnl)
    coef = -1.0_dp * outkron81_fnl !important -1

    deallocate( outkron9_1, outkron9_2, outkron81_fnl )

    dep_cart = (/a, b, c, d/)

    call this%Find_Lin_Dep(mu_ind, N2_ind, N3_ind, N4_ind, coef, row)

end subroutine LeftSidePntGr


subroutine RightSide(this, mu_dep, dep_N2, dep_N3, dep_N4, &
                   & ar, br, cr, dr, row)

    implicit none

    class(FCtyp)                                :: this

    integer, intent(in)                         :: mu_dep, dep_N2, dep_N3, dep_N4, &
                                                 & ar, br, cr, dr
    
    real(dp), dimension(:), intent(inout)       :: row

    ! =================================== Local Variables =================================== !

    real(dp)                                    :: dvar_tst
    real(dp), dimension(3**ordfc4)              :: coef
    integer, dimension(ordfc4)                  :: ind_Indx

    ! =================================== Local Variables =================================== !

    dvar_tst = this%F(mu_dep)%FCmu(1, dr, cr, br, ar, dep_N4, dep_N3, dep_N2)

    debug: if ( (dabs(dvar_tst-1.0_dp) < EPS_nzero) .or. (dabs(dvar_tst+5.0_dp) < EPS_nzero) .or. &
              & (dabs(dvar_tst-8.0_dp) < EPS_nzero) .or. (dabs(dvar_tst) < EPS_nzero) ) then

        write(*, 80) dvar_tst
        STOP

    end if debug

    ind_Indx = int( this%F(mu_dep)%FCmu(3:6, dr, cr, br, ar, dep_N4, dep_N3, dep_N2) )

    coef = this%F(mu_dep)%FCmu(7:, dr, cr, br, ar, dep_N4, dep_N3, dep_N2)

    call this%Find_Lin_Dep(ind_Indx(1), ind_Indx(2), ind_Indx(3), ind_Indx(4), coef, row)

    80 FORMAT("ERROR( in RightSide ): dvar_tst can not be other than -1 or -7 ", F6.3)

end subroutine RightSide


subroutine Find_lin_Dep(this, mu_ind, ind_N2, ind_N3, ind_N4, coef, row)

    implicit none

    class(FCtyp)                                :: this

    integer, intent(in)                         :: mu_ind, ind_N2, ind_N3, ind_N4
    real(dp), dimension(3**ordfc4), intent(in)  :: coef

    real(dp), dimension(:), intent(inout)       :: row

    ! ==================================== Local Variales ==================================== !

    real(dp)                                    :: dep_var1, cf
    integer                                     :: a_ind, res1, b_ind, res2, c_ind, d_ind, colno
    integer                                     :: i, ip

    ! ==================================== Local Variales ==================================== !

    coef_loop: do i = 1, (3**ordfc4)

        cf = coef(i)
        non_zero: if ( dabs(cf) > EPS ) then

            ip = (i - 1)

            a_ind = (ip / 27) + 1       !**!
            res1 = mod( ip, 27 )
            b_ind = (res1 / 9) + 1      !**!
            res2 = mod( res1, 9 )
            c_ind = (res2 / 3) + 1      !**!
            d_ind = mod( res2, 3 ) + 1  !**!

            dep_var1 = this%F(mu_ind)%FCmu(1, d_ind, c_ind, b_ind, a_ind, ind_N4, ind_N3, ind_N2)

            dep_var_chk: if ( dabs(dep_var1-1.0_dp) < EPS_nzero ) then
                call this%find_col_ind(mu_ind, ind_N2, ind_N3, ind_N4, a_ind, b_ind, c_ind, d_ind, colno)
                row(colno) = row(colno) + cf

            else if ( dabs(dep_var1+5.0_dp) < EPS_nzero ) then dep_var_chk
                call this%find_ind5(mu_ind, ind_N2, ind_N3, ind_N4, a_ind, b_ind, c_ind, d_ind, cf, row)

            else dep_var_chk

                not8: if ( dabs(dep_var1-8.0_dp) > EPS_nzero ) then
                    write(*, 90) dep_var1
                    STOP

                end if not8

            end if dep_var_chk

        end if non_zero

    end do coef_loop

    90 FORMAT("ERROR( in Find_lin_Dep ): dep_var1 can not be other than 1, -5 or 8 ", F6.3)

end subroutine Find_lin_Dep


subroutine find_ind5(this, mu, N2, N3, N4, a, b, c, d, cof, row)

    implicit none
    class(FCtyp)                                :: this

    integer, intent(in)                         :: mu, N2, N3, N4, a, b, c, d
    real(dp), intent(in)                        :: cof
    real(dp), dimension(:), intent(inout)       :: row

    ! ==================================== Local Variales ==================================== !

    real(dp), dimension(3**ordfc4)              :: coef_dep
    real(dp)                                    :: cf
    integer                                     :: a_ind, res1, b_ind, res2, c_ind, d_ind, colno
    integer                                     :: i, ip

    ! ==================================== Local Variales ==================================== !

    coef_dep = this%F(mu)%FCmu(7:, d, c, b, a, N4, N3, N2)

    debug: if ( all( dabs(coef_dep) < EPS_nzero ) ) then
        write(*, 120) 
        STOP

    end if debug

    coef_depLoop: do i = 1, 3**ordfc4

        cf = coef_dep(i)

        nonZeroChk: if ( dabs(cf) > EPS ) then
            ip = (i - 1)

            a_ind = (ip / 27) + 1       !**!
            res1 = mod( ip, 27 )
            b_ind = (res1 / 9) + 1      !**!
            res2 = mod( res1, 9 )
            c_ind = (res2 / 3) + 1      !**!
            d_ind = mod( res2, 3 ) + 1  !**!

            call this%find_col_ind(mu, N2, N3, N4, a_ind, b_ind, c_ind, d_ind, colno)
            row(colno) = row(colno) + (cof*cf)

        end if nonZeroChk

    end do coef_depLoop

    120 FORMAT("ERROR( in find_ind5 ): All coefficients can not be zero for dep_var=-5")

end subroutine find_ind5


subroutine find_col_ind(this, mu, N2, N3, N4, a_ind, b_ind, c_ind, d_ind, colno)

    implicit none
    class(FCtyp)                                :: this

    integer, intent(in)                         :: mu, N2, N3, N4, a_ind, b_ind, c_ind, d_ind
    integer, intent(out)                        :: colno

    ! ==================================== Local Variales ==================================== !

    real(dp)                                    :: dep_var1

    ! ==================================== Local Variales ==================================== !

    dep_var1 = this%F(mu)%FCmu(1, d_ind, c_ind, b_ind, a_ind, N4, N3, N2)

    depVarChk: if ( dabs(dep_var1-1.0_dp) < EPS_nzero ) then
        colno = this%IndFCPos(mu)%Pos_mu(d_ind, c_ind, b_ind, a_ind, N4, N3, N2)

        FillvalChk: if ( colno == FILLVAL ) then
            write(*, 145) 
            ERROR STOP

        end if FillvalChk

    else depVarChk
        write(*, 165) dep_var1
        ERROR STOP

    end if depVarChk

    145 FORMAT("ERROR( in find_col_ind ): Position of independent FC can not be FILLVAL")
    165 FORMAT("ERROR( in find_col_ind ): dep_var or independent FC can not be other than 1", F6.3)

end subroutine find_col_ind


subroutine ChunkSizeCheck(this)

    implicit none

    class(FCtyp)                                :: this

    ! ==================================== Local Variales ==================================== !

    real(dp)                                    :: InterMat_Sz
    integer, dimension(2)                       :: MatShape

    ! ==================================== Local Variales ==================================== !

    MatShape = shape( this%InterMat )
    InterMat_Sz = dble(MatShape(1)) * dble(MatShape(2)) * 8.0_dp &
                & / (1024.0_dp ** 2) !in MB

    SizeChk: if ( InterMat_Sz > this%InterMatChnkSize ) then

        this%InterMatChunk = this%InterMatChunk + 1
        !Write the chunk
        call this%WriteInterMatChunk(this%InterMat, this%InterMatChunk, this%InterMatFile)

        deallocate( this%InterMat )

    end if SizeChk

end subroutine ChunkSizeCheck

