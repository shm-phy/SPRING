
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


module FCmod


    use hdf5,               only : hid_t, hsize_t, H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER,&
                                 & H5F_ACC_TRUNC_F, H5F_ACC_RDWR_F, h5open_f, h5fcreate_f, &
                                 & h5screate_simple_f, h5dcreate_f,  h5dwrite_f, &
                                 & h5dclose_f, h5sclose_f, h5fclose_f, h5close_f, h5tclose_f, &
                                 & HID_T, HSIZE_T, H5F_ACC_RDONLY_F, h5dget_type_f, &
                                 & h5dget_space_f, h5sget_simple_extent_dims_f, h5dread_f, &
                                 & h5dopen_f, h5fopen_f, h5gcreate_f, h5gopen_f, h5gclose_f
    use kinds,              only : dp
    use constants,          only : EPS_nzero, ordfc2, &
                                 & FILLVAL, EPS, initMatsz2, &
                                 & Iden3, dummy_Comb
    use timer_class,        only : timer
    use unit_cell,          only : cell
    use Symmetry,           only : PntPer_typ
    use CutInfo,            only : CutOff, CellBasis
    use Helper,             only : TransformAtoms, get_atom_num, &
                                 & kron_prod
    use mklWrap,            only : LU_decompose

    implicit none
    private

    type, public        :: FCarr_mu

        real(dp), dimension(:,:,:,:), allocatable               :: FCmu

    end type FCarr_mu

    type, public        :: PosIndFC_mu

        integer, dimension(:,:,:), allocatable                  :: Pos_mu

    end type PosIndFC_mu

    type, public        :: FCtyp

        PRIVATE
        type(FCarr_mu), dimension(:), allocatable               :: F
        type(PosIndFC_mu), dimension(:), allocatable            :: IndFCPos

        integer                                                 :: fcount, ind_f, dep_f, &
                                                                 & zero_f, ignoreterms, &
                                                                 & IdenCount, InterTerm, &
                                                                 & InterTerm2, TotalFCNum
        logical                                                 :: NonAdvancWrite

        real(dp), dimension(:,:), allocatable                   :: GaussElMat 

        integer, dimension(:), allocatable                      :: PntGrInter
        integer, dimension(1)                                   :: PerInter
        integer                                                 :: NumInterPer, NumInterPnt

        integer, dimension(:), allocatable                      :: PntGrIter
        integer, dimension(1)                                   :: PerIter

        integer                                                 :: Pos_flat_size, PosIndCntr
        integer, dimension(:,:), allocatable                    :: PosIndx_flat

        integer                                                 :: InterMatChunk, InterMatPos
        real(dp)                                                :: InterMatChnkSize
        real(dp), dimension(:,:), allocatable                   :: InterMat
        character(len=128)                                      :: InterMatFile 

        real(dp), dimension(:,:), allocatable                   :: asr_mat
        integer                                                 :: N_ASR, asr_zero

        !-! type(FCarr_mu), dimension(:), allocatable               :: F_fin

        real(dp), dimension(:,:), allocatable                   :: RedInFC

        integer                                                 :: Nzrofin, Depfin, Indfin, Ignorefin

        type(PosIndFC_mu), dimension(:), allocatable            :: PosIndx_fin
        integer, dimension(:,:), allocatable                    :: PosIndx_flat_fin
        integer                                                 :: PosIndCntr_fin

        contains

            procedure, public, pass         :: init_FCtyp, Find_Independent, ApplyASR, ReducePntPerASR, &
                                             & Saveh5_F

            procedure, private, pass        :: PrimitiveExists, ApplySymm, Fillup_cartesian, Permute, &
                                             & MakeGaussMatTyp3, Exists, Fillup_cart_dep_perm, MakeGaussMatTyp4, &
                                             & MakeGaussMatTyp2, PointGroup, MakeGaussMatTyp1, Fillup_cart_dep, &
                                             & cart_dep, Create_Index_ind, Extend_Mat, CreateInterGmat, LeftSidePer, &
                                             & LeftSidePntGr, RightSide, Find_lin_Dep, find_ind5, find_col_ind, &
                                             & ChunkSizeCheck, IterFC, find_ind2, find_ind5Red, FindLC_ASR, &
                                             & PntPerPreReduced, DepIndx5_ASR, ReadChunks, MakeIndIndxFin
            
            procedure, private, nopass      :: WriteInterMatChunk, AllocationError

    end type FCtyp

contains

    subroutine init_FCtyp(this, ChnkSz, AdvncOut, sys, Cut_Off, PntPer)

        implicit none

        class(FCtyp)                                :: this

        real(dp), intent(in)                        :: ChnkSz
        logical, intent(in)                         :: AdvncOut
        type(cell), intent(in)                      :: sys
        type(CutOff), intent(in)                    :: Cut_Off
        type(PntPer_typ), intent(in)                :: PntPer

        ! =========================== Local variables ============================ !

        integer                                     :: Nbasis, Natm, cnst, NSymApply
        integer                                     :: mu, i

        ! =========================== Local variables ============================ !

        Nbasis = sys%natm

        allocate( this%F(Nbasis) )

        cnst = 2+ordfc2+(3**ordfc2)

        this%TotalFCNum = 0
        mu_loop1: do mu = 1, Nbasis

            Natm = Cut_Off%numAtom(mu)
            allocate( this%F(mu)%FCmu(cnst,3,3,Natm) )

            this%F(mu)%FCmu = 0.0_dp

            this%TotalFCNum = this%TotalFCNum + (Natm * (3**ordfc2))

        end do mu_loop1

        allocate( this%IndFCPos(Nbasis) )

        mu_loop2: do mu = 1, Nbasis

            Natm = Cut_Off%numAtom(mu)
            allocate( this%IndFCPos(mu)%Pos_mu(3,3,Natm) )

            this%IndFCPos(mu)%Pos_mu = FILLVAL

        end do mu_loop2

        this%fcount = 0
        this%ind_f = 0
        this%dep_f = 0
        this%zero_f = 0
        this%ignoreterms = 0
        this%IdenCount = 0
        this%InterTerm = 0
        this%InterTerm2 = 0

        this%NonAdvancWrite = ( .not. (AdvncOut) )

        NSymApply = (PntPer%Nsym - 1)

        allocate( this%PntGrInter( NSymApply ) )
        this%PntGrInter = [ (i, i=2,PntPer%Nsym) ]
        this%PerInter = [ (i, i=2,2) ]

        allocate( this%PntGrIter( NSymApply ) )
        this%PntGrIter = [ (i, i=2,PntPer%Nsym) ]
        this%PerIter = [ (i, i=2,2) ]

        this%Pos_flat_size = initMatsz2
        this%PosIndCntr = 0
        allocate( this%PosIndx_flat( (2*ordfc2), this%Pos_flat_size ) )
        this%PosIndx_flat = 0

        this%InterMatChunk = 0
        this%InterMatPos = 0
        this%InterMatChnkSize = ChnkSz ! in MB
        this%InterMatFile = 'InterMatChunk_2nd_'//trim(sys%prefix)//'.h5'

    end subroutine init_FCtyp


    subroutine Find_Independent(this, sys, Cut_Off, PntPer)

        implicit none

        class(FCtyp)                                :: this

        type(cell), intent(in)                      :: sys
        type(CutOff), intent(in)                    :: Cut_Off
        type(PntPer_typ), intent(in)                :: PntPer

        ! =========================== Local variables ============================ !

        integer                                     :: Nbasis, Natm
        integer                                     :: mu, N2, NumCol
        logical                                     :: PrmtvExst
        integer, dimension(:,:), allocatable        :: PosIndxFlat_tmp

        ! =========================== Local variables ============================ !

        Nbasis = sys%natm

        write(*, *)
        mu_loop: do mu = 1, Nbasis

            Natm = Cut_Off%numAtom(mu)
            N2_loop: do N2 = 1, Natm

                call this%PrimitiveExists(mu, N2, PrmtvExst)

                apply_pntper: if (.not. PrmtvExst) then

                    call this%ApplySymm(mu, N2, PntPer, Cut_Off, sys)

                end if apply_pntper

            end do N2_loop

        end do mu_loop

        dealoc_chk: if ( allocated(this%InterMat) ) then
            
            this%InterMatChunk = this%InterMatChunk + 1
            call this%WriteInterMatChunk(this%InterMat, this%InterMatChunk, this%InterMatFile)
            deallocate( this%InterMat )

        end if dealoc_chk

        ResizePosIndxFlat: if ( this%PosIndCntr == this%ind_f ) then

            NumCol = this%ind_f

            allocate( PosIndxFlat_tmp((2*ordfc2), NumCol) )
            PosIndxFlat_tmp(:, :) = this%PosIndx_flat(:, 1:NumCol)
            deallocate( this%PosIndx_flat )

            !call move_alloc(PosIndxFlat_tmp, this%PosIndx_flat)
            allocate( this%PosIndx_flat((2*ordfc2), NumCol) )
            this%PosIndx_flat = PosIndxFlat_tmp
            deallocate( PosIndxFlat_tmp )

        else ResizePosIndxFlat

            write(*, 50) this%PosIndCntr, this%ind_f
            STOP

        end if ResizePosIndxFlat

        write(*, *)
        write(*, 40) this%InterMatChunk

        call this%IterFC(Cut_Off, Nbasis)

        40 FORMAT("Number of InterMat chunks written: ", I5)
        50 FORMAT("ERROR( in Find_Independent ): this%PosIndCntr /= this%ind_f ", I5, I5)

    end subroutine Find_Independent


    subroutine PrimitiveExists(this, mu, N2, PrmtvExst)

        implicit none
        class(FCtyp)                    :: this
        
        integer, intent(in)             :: mu, N2
        logical, intent(out)            :: PrmtvExst

        prmtv_chk: if (this%fcount == 0) then
            PrmtvExst = .false.

        else if ( all(dabs(this%F(mu)%FCmu(2,:,:,N2)) < EPS) ) then prmtv_chk
            PrmtvExst = .false.

        else if ( all(dabs(this%F(mu)%FCmu(2,:,:,N2) - 1.0_dp) < EPS) &
           & .or. all(dabs(this%F(mu)%FCmu(2,:,:,N2) + 1.0_dp) < EPS) &
           & .or. all(dabs(this%F(mu)%FCmu(2,:,:,N2) + 7.0_dp) < EPS) ) then prmtv_chk

           PrmtvExst = .true.

        else  prmtv_chk
            write(*, 12)
            STOP

        end if prmtv_chk

        12 FORMAT("ERROR(in PrimitiveExists) : Illegal value of dep_var encountered")

    end subroutine PrimitiveExists


    subroutine ApplySymm(this, mu_ind, N2_ind, PntPer, Cut_Off, sys)

        implicit none

        class(FCtyp)                                :: this

        integer, intent(in)                         :: mu_ind, N2_ind
        type(PntPer_typ), intent(in)                :: PntPer
        type(CutOff), intent(in)                    :: Cut_Off
        type(cell), intent(in)                      :: sys

        ! =================================== Local variables =================================== !

        logical                                 :: InterGMat, writetype2GMat

        ! =================================== Local variables =================================== !


        call this%Fillup_cartesian(mu_ind, N2_ind, 1)

        this%IdenCount = 0
        InterGMat = .false.
        writetype2GMat = .true.

        this%NumInterPer = 0
        !this%PerInter(:) = 0
        call this%Permute(mu_ind, N2_ind, writetype2GMat, InterGMat, PntPer, Cut_Off)

        this%NumInterPnt = 0
        !this%PntGrInter(:) = 0
        call this%PointGroup(mu_ind, N2_ind, writetype2GMat, InterGMat, PntPer, Cut_Off, sys)

        gaus_el_chk: if ( this%IdenCount /= 0 ) then

            call this%cart_dep(mu_ind, N2_ind)

            deallocate( this%GaussElMat )
            this%IdenCount = 0

        end if gaus_el_chk

        this%InterTerm2 = this%InterTerm2 + (this%NumInterPer + this%NumInterPnt)

        call this%Create_Index_ind(mu_ind, N2_ind)

        ! ----------------------------------- Second Iteration ----------------------------------- !
        InterGMat = .true.
        writetype2GMat = .false.

        chk1: if ( (this%NumInterPer + this%NumInterPnt) /= 0 ) then

            call this%Extend_Mat()

            chk2: if ( this%NumInterPer /= 0 ) then

                call this%Permute(mu_ind, N2_ind, writetype2GMat, InterGMat, PntPer, Cut_Off)

            end if chk2

            chk3: if ( this%NumInterPnt /= 0 ) then

                call this%PointGroup(mu_ind, N2_ind,  writetype2GMat, InterGMat, PntPer, Cut_Off, sys)

            end if chk3

            call this%ChunkSizeCheck()

        end if chk1
        ! ----------------------------------- Second Iteration ----------------------------------- !

        NonAdvancingOut: if ( this%NonAdvancWrite ) then
            write(6, 60, advance="no") char(13), this%fcount, this%TotalFCNum, this%ind_f
            FLUSH(6)

        else NonAdvancingOut
            write(*, 50) this%fcount, this%TotalFCNum, this%ind_f

        end if NonAdvancingOut

        50 FORMAT("No. of FC covered: ", I8, "/", I8, " /|\ No. of independent FC: ", I6)
        60 FORMAT(1a1, "No. of FC covered: ", I8, "/", I8, " /|\ No. of independent FC: ", I6)

    end subroutine ApplySymm


    subroutine Fillup_cartesian(this, mu_ind, N2_ind, dep_var)

        implicit none

        class(FCtyp)                                :: this

        integer, intent(in)                         :: mu_ind, N2_ind, dep_var

        ! =================================== Local variables =================================== !

        integer                             :: alpha, beta
        integer                             :: i_ind
        integer, dimension(ordfc2)          :: ind_Indx
        real(dp), dimension(3**ordfc2)      :: coef

        ! =================================== Local variables =================================== !

        ind_Indx = (/mu_ind, N2_ind/)

        alpha_loop: do alpha = 1, 3
            beta_loop: do beta = 1, 3

                this%F(mu_ind)%FCmu(1:2, beta, alpha, N2_ind) = dble( dep_var )
                this%F(mu_ind)%FCmu(3:4, beta, alpha, N2_ind) = dble( ind_Indx )

                coef = 0.0_dp

                i_ind = 3*(alpha-1) + beta !- 1) + 1

                coef(i_ind) = coef(i_ind) + 1.0_dp

                this%F(mu_ind)%FCmu(5:, beta, alpha, N2_ind) = coef

                this%fcount = this%fcount + 1
                this%ind_f = this%ind_f + 1

            end do beta_loop
        end do alpha_loop

    end subroutine Fillup_cartesian


    subroutine Permute(this, mu_ind, N2_ind, &
                     & writetype2GMat, InterGMat, PntPer, Cut_Off)

        implicit none

        class(FCtyp)                                :: this

        integer, intent(in)                         :: mu_ind, N2_ind
                                                    
        logical, intent(in)                         :: writetype2GMat
        logical, intent(in)                         :: InterGMat

        type(PntPer_typ), intent(in)                :: PntPer
        type(CutOff), intent(in)                    :: Cut_Off

        ! =================================== Local variables =================================== !

        integer                                     :: NumPerIndx, err_var
        integer, dimension(1)                       :: PerIndx
        integer, dimension(ordfc2)                  :: comb
        integer, dimension(4, ordfc2)               :: ind_array

        integer, dimension(4)                       :: dep_N1crt_mu, dep_N2crt_nu

        integer                                     :: mu_dep, dep_N2

        integer, dimension(3)                       :: N1cell
        logical                                     :: outside, exist_chk

        integer                                     :: i, indx

        ! =================================== Local variables =================================== !

        ind_array(1:3, 1) = (/0, 0, 0/)
        ind_array(4, 1) = mu_ind

        ind_array(:, 2) = Cut_Off%atmIndx(mu_ind)%cbIndx(:, N2_ind)

        chk_interGMat: if ( InterGMat ) then
            NumPerIndx = this%NumInterPer
            PerIndx = this%PerInter

        else chk_interGMat
            NumPerIndx = size( PntPer%Per, 2) - 1 !1
            PerIndx = this%PerIter

        end if chk_interGMat

        perindx_loop: do i = 1, NumPerIndx

            indx = PerIndx(i)
            comb = PntPer%Per(:, indx)

            dep_N1crt_mu = ind_array(:, comb(1))
            dep_N2crt_nu = ind_array(:, comb(2))

            N1cell = dep_N1crt_mu(1:3)
            mu_dep = dep_N1crt_mu(4)

            dep_N1crt_mu(1:3) = dep_N1crt_mu(1:3) - N1cell
            dep_N2crt_nu(1:3) = dep_N2crt_nu(1:3) - N1cell

            dep_N2 = get_atom_num(mu_dep, dep_N2crt_nu, Cut_Off)

            outside = (dep_N2 == 0) 

            out_chk: if ( .not. outside ) then

                same_map: if ( (mu_dep == mu_ind) .and. (dep_N2 == N2_ind) ) then
                    
                    chk_interGmat2: if ( .not. InterGMat ) then

                        call this%MakeGaussMatTyp3(comb)

                    else chk_interGmat2
                        write(*, 35)
                        STOP

                    end if chk_interGmat2

                else same_map

                    call this%Exists(mu_dep, dep_N2, -7, mu_ind, N2_ind, comb, indx, &
                                   & InterGMat, writetype2GMat, PntPer, exist_chk, err_var)

                    err_chk: if ( err_var /= 0 ) then
                        write(*, 76) err_var
                        STOP
                    end if err_chk

                    fillupdep_Per: if ( (.not. exist_chk) .and. (.not. InterGMat) ) then

                        call this%Fillup_cart_dep_perm(mu_dep, dep_N2, -7, mu_ind, N2_ind, comb)

                    end if fillupdep_Per

                end if same_map

            end if out_chk

        end do perindx_loop

        35 FORMAT("ERROR(in Permute): Second time enters for Gauss Mat writing")
        76 FORMAT("ERROR(in Permute): Unchecked condition in subroutine Exists ", I2)

    end subroutine Permute


    subroutine MakeGaussMatTyp3(this, comb)

        implicit none

        class(FCtyp)                                :: this

        integer, dimension(ordfc2), intent(in)      :: comb

        ! ===================================== Local Variables ===================================== !

        integer                                     :: alpha, beta, &
                                                     & i_ind, i_dep, af
                                            
        integer, dimension(ordfc2)                  :: ind_cart, dep_cart
        real(dp), dimension(3**ordfc2, 3**ordfc2)   :: GaussMatloc
        real(dp), dimension(3**ordfc2)              :: coef

        real(dp), dimension(:, :), allocatable      :: GaussMatcopy
        integer                                     :: Matsize2

        ! ===================================== Local Variables ===================================== !

        GaussMatloc(:,:) = 0.0_dp

        alpha_loop: do alpha = 1, 3
            beta_loop: do beta = 1, 3

                ind_cart = (/alpha, beta/)

                atm4_loop: do af = 1, ordfc2
                    dep_cart(af) = ind_cart( comb(af) )
                end do atm4_loop

                coef = 0.0_dp

                i_ind = 3*(ind_cart(1) - 1) + ind_cart(2) !-1 + 1
                coef(i_ind) = coef(i_ind) + 1.0_dp

                i_dep = 3*(dep_cart(1) - 1) + dep_cart(2) !-1 + 1
                coef(i_dep) = coef(i_dep) - 1.0_dp

                GaussMatloc(:, i_ind) = coef

            end do beta_loop
        end do alpha_loop

        gaussMatExist: if ( this%IdenCount == 0) then
            allocate( this%GaussElMat(3**ordfc2, 3**ordfc2) )
            this%GaussElMat = GaussMatloc

        else gaussMatExist

            Matsize2 = size(this%GaussElMat, 2)

            !call move_alloc(this%GaussElMat, GaussMatcopy)
            allocate( GaussMatcopy( 3**ordfc2, Matsize2 ) )
            GaussMatcopy = this%GaussElMat
            deallocate( this%GaussElMat )

            allocate( this%GaussElMat( 3**ordfc2, (Matsize2+3**ordfc2) ) )
            this%GaussElMat(:, 1:Matsize2) = GaussMatcopy
            this%GaussElMat(:, (Matsize2+1): ) = GaussMatloc

            deallocate( GaussMatcopy )

        end if gaussMatExist

        this%IdenCount = this%IdenCount + 1

    end subroutine MakeGaussMatTyp3


    subroutine Exists(this, mu_dep, dep_N2, dep_var, mu_ind, N2_ind, comb, mat_num, &
                          & InterGMat, writetype2GMat, PntPer, exist_chk, info)

        implicit none

        class(FCtyp)                                :: this

        integer, intent(in)                         :: mu_dep, dep_N2, dep_var, &
                                                     & mu_ind, N2_ind, mat_num

        !-! real(dp), dimension(3,3), intent(in)        :: mat
        integer, dimension(ordfc2), intent(in)      :: comb

        logical, intent(in)                         :: InterGMat, writetype2GMat

        type(PntPer_typ), intent(in)                :: PntPer
        logical, intent(out)                        :: exist_chk
        integer, intent(out)                        :: info

        ! ====================================== Local Variables ====================================== !

        real(dp)                                    :: var
        integer, dimension(ordfc2)                  :: ind_Indx, indx_chk

        ! ====================================== Local Variables ====================================== !

        var = this%F(mu_dep)%FCmu(2, 1, 1, dep_N2)

        indx_chk = int( this%F(mu_dep)%FCmu(3:4, 1, 1, dep_N2) )

        ind_Indx = (/mu_ind, N2_ind/)

        info = 1

        out_if: if ( dabs(var) < EPS ) then

            exist_chk = .false.
            info = 0

        else out_if

            out2_if: if ( all((indx_chk-ind_Indx) == 0) ) then

                inner_if1: if ( dep_var == -7 ) then

                    inner_if2: if ( dabs(var-dble(dep_var)) < EPS ) then

                        inner_if3: if ( .not. InterGMat ) then

                            inner_if4: if (writetype2GMat) then

                                call this%MakeGaussMatTyp4(mu_dep, dep_N2, comb)
                                
                            end if inner_if4

                            exist_chk = .true.
                            info = 0

                        else inner_if3
                            write(*, 32)
                            STOP

                        end if inner_if3

                    else inner_if2 !( .not. ( dabs(var-dble(dep_var)) < EPS ) )
                        write(*, 42) var
                        STOP

                    end if inner_if2

                else inner_if1 !( dep_var == -1 )

                    iner2_if2: if ( ( dabs(var-dble(dep_var)) < EPS ) .or. &
                                  & ( dabs(var+7.0_dp) < EPS ) ) then

                        iner2_if3: if ( .not. InterGMat ) then

                            iner2_if4: if (writetype2GMat) then

                                call this%MakeGaussMatTyp2(mu_dep, dep_N2, mat_num, PntPer)

                            end if iner2_if4

                            exist_chk = .true.
                            info = 0

                        else iner2_if3
                            write(*, 32)
                            STOP

                        end if iner2_if3

                    !*! debug !*!
                    !*! else iner2_if2
                    !*!     write(*, 52) var
                    !*!     STOP
                    !*! debug !*!

                    end if iner2_if2

                end if inner_if1

            !else if ( .not. ( all((indx_chk-ind_Indx) == 0) ) ) then out2_if
            else out2_if !( .not. ( all((indx_chk-ind_Indx) == 0) ) )

                iner3_if1: if ( ( dabs(var+1.0_dp) < EPS ) .or. &
                              & ( dabs(var+7.0_dp) < EPS ) ) then

                    iner3_if2: if ( .not. InterGMat ) then

                        iner3_if3: if ( dep_var == -7 ) then

                            this%NumInterPer = this%NumInterPer + 1
                            this%PerInter(this%NumInterPer) = mat_num

                            info = 0

                        else iner3_if3 !( dep_var == -1 )

                            this%NumInterPnt = this%NumInterPnt + 1
                            this%PntGrInter(this%NumInterPnt) = mat_num

                            info = 0

                        end if iner3_if3

                        exist_chk = .true. 

                    else iner3_if2 !( InterGMat )

                        call this%CreateInterGmat(mu_dep, dep_N2, dep_var, &
                                                & mu_ind, N2_ind, mat_num, PntPer, comb)

                        info = 0

                    end if iner3_if2
                
                !*! debug !*!
                !*! else iner3_if1
                !*!     write(*, 52) var
                !*!     STOP
                !*! debug !*!

                end if iner3_if1

                exist_chk = .true.

            end if out2_if

        end if out_if

        32 FORMAT("ERROR( in Exists ): Second time of iteration should not enter here")
        42 FORMAT("ERROR( in Exists ): Permutation is applied first, So var can not be other than -7", F5.2)
        !*! 52 FORMAT("ERROR( in Exists ): Illegal value of var", I3)

    end subroutine Exists


    subroutine Fillup_cart_dep_perm(this, mu_dep, dep_N2, dep_var, &
                                        & mu_ind, N2_ind, comb)

        implicit none

        class(FCtyp)                                :: this

        integer, intent(in)                         :: mu_dep, dep_N2, dep_var, &
                                                     & mu_ind, N2_ind 

        integer, dimension(ordfc2), intent(in)      :: comb

        ! ====================================== Local Variables ====================================== !

        integer                                     :: alpha, beta, &
                                                     & i_ind, af
                                            
        integer, dimension(ordfc2)                  :: ind_Indx
        integer, dimension(ordfc2)                  :: ind_cart, dep_cart
        real(dp), dimension(3**ordfc2)              :: coef

        ! ====================================== Local Variables ====================================== !

        ind_Indx = (/mu_ind, N2_ind/)

        alpha_loop: do alpha = 1, 3
            beta_loop: do beta = 1, 3

                ind_cart = (/alpha, beta/)
                coef = 0.0_dp

                i_ind = 3*(ind_cart(1) - 1) + ind_cart(2) !-1 + 1

                coef(i_ind) = coef(i_ind) + 1.0_dp

                atm4_loop: do af = 1, ordfc2
                    dep_cart(af) = ind_cart( comb(af) )
                end do atm4_loop

                this%F(mu_dep)%FCmu(1:2, dep_cart(2), dep_cart(1), dep_N2) = dble( dep_var )

                this%F(mu_dep)%FCmu(3:4, dep_cart(2), dep_cart(1), dep_N2) = dble( ind_Indx )

                this%F(mu_dep)%FCmu(5:, dep_cart(2), dep_cart(1), dep_N2) = coef

                this%fcount = this%fcount + 1
                this%dep_f = this%dep_f + 1

            end do beta_loop
        end do alpha_loop

    end subroutine Fillup_cart_dep_perm


    subroutine MakeGaussMatTyp4(this, mu_dep, dep_N2, comb)

        implicit none

        class(FCtyp)                                :: this

        integer, intent(in)                         :: mu_dep, dep_N2

        integer, dimension(ordfc2), intent(in)      :: comb

        ! ====================================== Local Variables ====================================== !

        integer                                     :: alpha, beta, &
                                                     & i_ind, af, i
                                            
        integer, dimension(ordfc2)                  :: ind_cart, dep_cart
        real(dp), dimension(3**ordfc2)              :: coef, coef1

        real(dp), dimension(3**ordfc2, 3**ordfc2)   :: GaussMatloc

        real(dp), dimension(:, :), allocatable      :: GaussMatcopy
        integer                                     :: Matsize2

        ! ====================================== Local Variables ====================================== !

        GaussMatloc(:,:) = 0.0_dp

        alpha_loop: do alpha = 1, 3
            beta_loop: do beta = 1, 3

                ind_cart = (/alpha, beta/)
                coef = 0.0_dp

                i_ind = 3*(ind_cart(1) - 1) + ind_cart(2) !-1 + 1

                coef(i_ind) = coef(i_ind) + 1.0_dp

                atm4_loop: do af = 1, ordfc2
                    dep_cart(af) = ind_cart( comb(af) )
                end do atm4_loop

                coef1 = this%F(mu_dep)%FCmu(5:, dep_cart(2), dep_cart(1), dep_N2)

                coef = (coef - coef1)

                i = 3*(alpha-1) + beta !- 1) + 1

                GaussMatloc(:, i) = coef

            end do beta_loop
        end do alpha_loop

        gaussMatExist: if ( this%IdenCount == 0) then
            allocate( this%GaussElMat(3**ordfc2, 3**ordfc2) )
            this%GaussElMat = GaussMatloc

        else gaussMatExist

            Matsize2 = size(this%GaussElMat, 2)

            !call move_alloc(this%GaussElMat, GaussMatcopy)
            allocate( GaussMatcopy( 3**ordfc2, Matsize2 ) )
            GaussMatcopy = this%GaussElMat
            deallocate( this%GaussElMat )

            allocate( this%GaussElMat( 3**ordfc2, (Matsize2+3**ordfc2) ) )
            this%GaussElMat(:, 1:Matsize2) = GaussMatcopy
            this%GaussElMat(:, (Matsize2+1): ) = GaussMatloc

            deallocate( GaussMatcopy )

        end if gaussMatExist

        this%IdenCount = this%IdenCount + 1

    end subroutine MakeGaussMatTyp4


    subroutine MakeGaussMatTyp2(this, mu_dep, dep_N2, mat_num, PntPer)

        implicit none

        class(FCtyp)                                :: this

        integer, intent(in)                         :: mu_dep, dep_N2, &
                                                     & mat_num

        type(PntPer_typ), intent(in)                :: PntPer

        ! ====================================== Local Variables ====================================== !

        integer                                     :: alpha, beta, i
                                            
        real(dp), dimension(3**ordfc2)              :: coef1, coef

        real(dp), dimension(3**ordfc2, 3**ordfc2)   :: GaussMatloc

        real(dp), dimension(3,3)                    :: matxyz
        real(dp), dimension(3)                      :: vec_a, vec_b
        real(dp), dimension(:), allocatable         :: outkron9_1

        real(dp), dimension(:, :), allocatable      :: GaussMatcopy
        integer                                     :: Matsize2

        ! ====================================== Local Variables ====================================== !

        GaussMatloc(:,:) = 0.0_dp

        matxyz = PntPer%Rxyz(:, :, mat_num)

        alpha_loop: do alpha = 1, 3
            vec_a = matxyz(alpha, :)

            beta_loop: do beta = 1, 3

                vec_b = matxyz(beta, :)
                call kron_prod(vec_a, vec_b, outkron9_1) 

                coef1 = this%F(mu_dep)%FCmu(5:, beta, alpha, dep_N2)

                coef = outkron9_1 - coef1

                i = 3*(alpha-1) + beta !- 1) + 1

                GaussMatloc(:, i) = coef

                deallocate( outkron9_1 )

            end do beta_loop

        end do alpha_loop

        gaussMatExist: if ( this%IdenCount == 0) then
            allocate( this%GaussElMat(3**ordfc2, 3**ordfc2) )
            this%GaussElMat = GaussMatloc

        else gaussMatExist

            Matsize2 = size(this%GaussElMat, 2)

            !call move_alloc(this%GaussElMat, GaussMatcopy)
            allocate( GaussMatcopy( 3**ordfc2, Matsize2 ) )
            GaussMatcopy = this%GaussElMat
            deallocate( this%GaussElMat )

            allocate( this%GaussElMat( 3**ordfc2, (Matsize2+3**ordfc2) ) )
            this%GaussElMat(:, 1:Matsize2) = GaussMatcopy
            this%GaussElMat(:, (Matsize2+1): ) = GaussMatloc

            deallocate( GaussMatcopy )

        end if gaussMatExist

        this%IdenCount = this%IdenCount + 1

    end subroutine MakeGaussMatTyp2


    subroutine PointGroup(this, mu_ind, N2_ind, &
                        & writetype2GMat, InterGMat, PntPer, Cut_Off, sys)

        implicit none

        class(FCtyp)                                :: this

        integer, intent(in)                         :: mu_ind, N2_ind
                                                    
        logical, intent(in)                         :: writetype2GMat
        logical, intent(in)                         :: InterGMat

        type(PntPer_typ), intent(in)                :: PntPer
        type(CutOff), intent(in)                    :: Cut_Off
        type(cell), intent(in)                      :: sys

        ! =================================== Local variables =================================== !

        integer                                     :: NumPntIndx, err_var
        integer, dimension(:), allocatable          :: PntIndx

        real(dp), dimension(3,3)                    :: Rot
        real(dp), dimension(3)                      :: Trans

        integer, dimension(4)                       :: ind_N1crt_mu, ind_N2crt_nu

        integer, dimension(4)                       :: dep_N1crt_mu, dep_N2crt_nu

        integer                                     :: mu_dep, dep_N2

        logical                                     :: outside, iden_trans, exist_chk

        integer                                     :: i, indx

        ! =================================== Local variables =================================== !

        ind_N1crt_mu(1:3) = (/0, 0, 0/)
        ind_N1crt_mu(4) = mu_ind

        ind_N2crt_nu(:) = Cut_Off%atmIndx(mu_ind)%cbIndx(:, N2_ind)

        allocate( PntIndx( size(this%PntGrIter) ) )

        chk_interGMat: if ( InterGMat ) then
            NumPntIndx = this%NumInterPnt
            PntIndx = this%PntGrInter

        else chk_interGMat
            NumPntIndx = size( this%PntGrIter ) 
            PntIndx = this%PntGrIter

        end if chk_interGMat

        pntSym_loop: do i = 1, NumPntIndx

            indx = PntIndx(i)
            Rot = PntPer%R(:,:,indx)
            Trans = PntPer%T(:, indx)

            dep_N1crt_mu = ind_N1crt_mu
            dep_N2crt_nu = ind_N2crt_nu
            
            call TransformAtoms(sys, dep_N1crt_mu, dep_N2crt_nu, Rot, Trans)
            
            mu_dep = dep_N1crt_mu(4)

            dep_N2 = get_atom_num(mu_dep, dep_N2crt_nu, Cut_Off)

            outside = (dep_N2 == 0) 

            out_chk: if ( .not. outside ) then 

                iden_trans = ( (mu_dep == mu_ind) .and. (dep_N2 == N2_ind) .and. &
                             & (.not. all(dabs(Iden3-Rot) < EPS_nzero)) ) 

                same_map: if ( iden_trans ) then

                    chk_interGmat2: if ( .not. InterGMat ) then

                        call this%MakeGaussMatTyp1(indx, PntPer)

                    else chk_interGmat2
                        write(*, 35)
                        STOP

                    end if chk_interGmat2

                else same_map

                    call this%Exists(mu_dep, dep_N2, -1, mu_ind, N2_ind, dummy_Comb, indx, & 
                                   & InterGMat, writetype2GMat, PntPer, exist_chk, err_var)

                    err_chk: if ( err_var /= 0 ) then
                        write(*, 76) err_var
                        STOP
                    end if err_chk

                    fillupdep_Per: if ( (.not. exist_chk) .and. (.not. InterGMat) ) then

                        call this%Fillup_cart_dep(mu_dep, dep_N2, -1, &
                                                & mu_ind, N2_ind, indx, PntPer)

                    end if fillupdep_Per

                end if same_map

            end if out_chk

        end do pntSym_loop

        deallocate( PntIndx )

        35 FORMAT("ERROR(in PointGroup): Second time enters for Gauss Mat writing")
        76 FORMAT("ERROR(in PointGroup): Unchecked condition in subroutine Exists ", I2)

    end subroutine PointGroup


    subroutine MakeGaussMatTyp1(this, mat_num, PntPer)

        implicit none

        class(FCtyp)                                :: this

        integer, intent(in)                         :: mat_num
        type(PntPer_typ), intent(in)                :: PntPer

        ! ====================================== Local Variables ====================================== !

        integer                                     :: alpha, beta, i
                                            
        real(dp), dimension(3**ordfc2, 3**ordfc2)   :: GaussMatloc

        real(dp), dimension(3,3)                    :: matxyz
        real(dp), dimension(3)                      :: vec_a, vec_b
        real(dp), dimension(:), allocatable         :: outkron9_1

        real(dp), dimension(:, :), allocatable      :: GaussMatcopy
        integer                                     :: Matsize2

        ! ====================================== Local Variables ====================================== !

        GaussMatloc(:,:) = 0.0_dp

        matxyz = PntPer%Rxyz(:, :, mat_num)

        alpha_loop: do alpha = 1, 3
            vec_a = matxyz(alpha, :)

            beta_loop: do beta = 1, 3

                vec_b = matxyz(beta, :)
                call kron_prod(vec_a, vec_b, outkron9_1) 

                i = 3*(alpha-1) + beta !- 1) + 1

                outkron9_1(i) = outkron9_1(i) - 1.0_dp

                GaussMatloc(:, i) = outkron9_1

                deallocate( outkron9_1 )

            end do beta_loop

        end do alpha_loop

        gaussMatExist: if ( this%IdenCount == 0) then

            allocate( this%GaussElMat(3**ordfc2, 3**ordfc2) )
            this%GaussElMat = GaussMatloc

        else gaussMatExist

            Matsize2 = size(this%GaussElMat, 2)

            !call move_alloc(this%GaussElMat, GaussMatcopy)
            allocate( GaussMatcopy( 3**ordfc2, Matsize2 ) )
            GaussMatcopy = this%GaussElMat
            deallocate( this%GaussElMat )

            allocate( this%GaussElMat( 3**ordfc2, (Matsize2+3**ordfc2) ) )
            this%GaussElMat(:, 1:Matsize2) = GaussMatcopy
            this%GaussElMat(:, (Matsize2+1): ) = GaussMatloc

            deallocate( GaussMatcopy )

        end if gaussMatExist

        this%IdenCount = this%IdenCount + 1

    end subroutine MakeGaussMatTyp1


    subroutine Fillup_cart_dep(this, mu_dep, dep_N2, dep_var, &
                                   & mu_ind, N2_ind, mat_num, PntPer)


        implicit none

        class(FCtyp)                                :: this

        integer, intent(in)                         :: mu_dep, dep_N2, dep_var, &
                                                     & mu_ind, N2_ind, mat_num

        type(PntPer_typ), intent(in)                :: PntPer

        ! ====================================== Local Variables ====================================== !

        integer                                     :: alpha, beta
                                            
        integer, dimension(ordfc2)                  :: ind_Indx

        real(dp), dimension(3,3)                    :: matxyz
        real(dp), dimension(3)                      :: vec_a, vec_b
        real(dp), dimension(:), allocatable         :: outkron9_1

        ! ====================================== Local Variables ====================================== !

        ind_Indx = (/mu_ind, N2_ind/)

        matxyz = PntPer%Rxyz(:, :, mat_num)

        alpha_loop: do alpha = 1, 3
            vec_a = matxyz(alpha, :)

            beta_loop: do beta = 1, 3

                vec_b = matxyz(beta, :)
                call kron_prod(vec_a, vec_b, outkron9_1) 

                this%F(mu_dep)%FCmu(1:2, beta, alpha, dep_N2) = dble( dep_var )
                this%F(mu_dep)%FCmu(3:4, beta, alpha, dep_N2) = dble( ind_Indx )
                this%F(mu_dep)%FCmu(5:,  beta, alpha, dep_N2) = outkron9_1

                this%fcount = this%fcount + 1
                this%dep_f = this%dep_f + 1

                deallocate( outkron9_1 )

            end do beta_loop

        end do alpha_loop

    end subroutine Fillup_cart_dep


    subroutine cart_dep(this, mu_ind, ind_N2)

        implicit none

        class(FCtyp)                                :: this

        integer, intent(in)                         :: mu_ind, ind_N2

        ! ===================================== Local Variables ===================================== !

        real(dp), dimension(:,:), allocatable       :: GaussElMatin_T
        real(dp), dimension(:,:), allocatable       :: GausMatT

        integer                                     :: num_fc, alpha, beta, &
                                                     & alpha_lin, beta_lin

        integer, dimension(ordfc2)                  :: ind_Indx
        real(dp)                                    :: pivot, chk2, coeff_lin, chk_lin, indx_final, &
                                                     & chk

        real(dp), dimension(:), allocatable         :: after
        real(dp), dimension(3**ordfc2)              :: coef, coef_dep, coef_final

        integer                                     :: fc, fc_p, num_el, ind, fc_lin, fc_lin_p, p

        ! ===================================== Local Variables ===================================== !

        num_fc = 3**ordfc2
        ind_Indx = (/mu_ind, ind_N2/)

        allocate( GaussElMatin_T( size(this%GaussElMat, 2), size(this%GaussElMat, 1) ) )
        GaussElMatin_T = transpose( this%GaussElMat )
        !~! deallocate( this%GaussElMat )

        call LU_decompose(GaussElMatin_T, GausMatT)
        deallocate( GaussElMatin_T )

        fc_cart_loop: do fc = num_fc, 1, -1

            fc_p = fc - 1

            alpha = (fc_p / 3) + 1      !**!
            beta = mod( fc_p, 3 ) + 1   !**!

            pivot = GausMatT(fc, fc)

            last_fc: if ( fc == num_fc ) then
                
                pivot_chk1: if ( dabs(pivot) > EPS ) then

                    chk2 = this%F(mu_ind)%FCmu(1, beta, alpha, ind_N2)

                    write8_1: if ( dabs(chk2-1.0_dp) < EPS ) then

                        this%F(mu_ind)%FCmu(1, beta, alpha, ind_N2) = 8.0_dp
                        this%F(mu_ind)%FCmu(5:, beta, alpha, ind_N2) = 0.0_dp

                        this%zero_f = this%zero_f + 1
                        this%ind_f = this%ind_f - 1

                    else if ( dabs(chk2-8.0_dp) < EPS ) then write8_1

                        write(*, 44)

                    else write8_1

                        write(*, 54) chk2
                        STOP

                    end if write8_1

                else if ( dabs(pivot) < EPS ) then pivot_chk1

                    chk2 = this%F(mu_ind)%FCmu(1, beta, alpha, ind_N2)

                    not_ind_chk1: if ( dabs(chk2-1.0_dp) > EPS ) then

                        write(*, 54) chk2
                        STOP

                    else not_ind_chk1

                        ! ~ Overwrite ~ !
                        this%F(mu_ind)%FCmu(3:4, beta, alpha, ind_N2) = dble( ind_Indx )

                    end if not_ind_chk1

                end if pivot_chk1

            else last_fc

                num_el = ( num_fc - fc )

                allocate( after(num_el) )

                after = GausMatT( (fc+1):num_fc, fc )

                pivot_chk2: if ( dabs(pivot) > EPS ) then

                    any_aftr: if ( any(dabs(after) > EPS) ) then

                        coef = 0.0_dp

                        all_el_loop: do ind = 1, num_el

                            nonzero_el: if ( dabs(after(ind)) > EPS ) then

                                fc_lin = fc + ind
                                fc_lin_p = (fc_lin - 1)

                                alpha_lin = (fc_lin_p / 3) + 1      !**!
                                beta_lin = mod( fc_lin_p, 3 ) + 1   !**!

                                coeff_lin = -1.0_dp * after(ind) / pivot

                                chk_lin = this%F(mu_ind)%FCmu(1, beta_lin, alpha_lin, ind_N2)

                                chk_lin_cond: if ( dabs(chk_lin - 1.0_dp) < EPS ) then

                                    coef(fc_lin) = coef(fc_lin) + coeff_lin

                                else if ( dabs(chk_lin - 8.0_dp) < EPS ) then chk_lin_cond

                                    coef(fc_lin) = coef(fc_lin) + (coeff_lin * 0.0_dp)

                                else if ( dabs(chk_lin + 5.0_dp) < EPS ) then chk_lin_cond

                                    coef_dep = this%F(mu_ind)%FCmu(5:, beta_lin, alpha_lin, ind_N2)

                                    coef_dep_el: do p = 1, num_fc

                                        nonzro_coef_dep: if ( dabs(coef_dep(p)) > EPS ) then

                                            debug1: if ( .not. (p > fc) ) then

                                                write(*, 76) p, fc
                                                STOP

                                            else debug1

                                                coef(p) = coef(p) + ( coeff_lin * coef_dep(p) )

                                            end if debug1

                                        end if nonzro_coef_dep

                                    end do coef_dep_el

                                else chk_lin_cond

                                    write(*, 54) chk_lin
                                    STOP

                                end if chk_lin_cond

                            end if nonzero_el

                        end do all_el_loop

                        this%F(mu_ind)%FCmu(1, beta, alpha, ind_N2) = -5.0_dp
                        ! ~ Overwrite ~ !
                        this%F(mu_ind)%FCmu(3:4, beta, alpha, ind_N2) = dble( ind_Indx )
                        this%F(mu_ind)%FCmu(5:,  beta, alpha, ind_N2) = coef

                        this%ind_f = this%ind_f - 1
                        this%dep_f = this%dep_f + 1

                    !else any_aftr !( all( dabs(after) < EPS ) )
                    else if ( all( dabs(after) < EPS ) ) then any_aftr

                        chk2 = this%F(mu_ind)%FCmu(1, beta, alpha, ind_N2)
                        
                        write8_2: if ( dabs(chk2-1.0_dp) < EPS ) then
                        
                            this%F(mu_ind)%FCmu(1, beta, alpha, ind_N2) = 8.0_dp
                            this%F(mu_ind)%FCmu(5:,beta, alpha, ind_N2) = 0.0_dp
                        
                            this%zero_f = this%zero_f + 1
                            this%ind_f = this%ind_f - 1
                        
                        else if ( dabs(chk2-8.0_dp) < EPS ) then write8_2
                        
                            write(*, 44)
                        
                        else write8_2
                        
                            write(*, 54) chk2
                            STOP
                        
                        end if write8_2

                    end if any_aftr

                !else pivot_chk !( dabs(pivot) < EPS )
                else if ( dabs(pivot) < EPS ) then pivot_chk2

                    chk = this%F(mu_ind)%FCmu(1, beta, alpha, ind_N2)

                    chk_not_ind: if ( dabs(chk-1.0_dp) > EPS ) then
                        write(*, 54) chk
                        STOP

                    else chk_not_ind
                        ! ~ Overwrite ~ !
                        this%F(mu_ind)%FCmu(3:4, beta, alpha, ind_N2) = dble( ind_Indx )

                    end if chk_not_ind

                    ignore_cond: if ( any( dabs(after)  > EPS ) ) then

                        this%ignoreterms = this%ignoreterms + 1

                    end if ignore_cond

                end if pivot_chk2

                deallocate( after )

                indx_final = this%F(mu_ind)%FCmu(1, beta, alpha, ind_N2)
                coef_final = this%F(mu_ind)%FCmu(5:,beta, alpha, ind_N2)

                final_zro_chk: if ( (dabs(indx_final + 5.0_dp) < EPS) .and. &
                                  & (all( dabs(coef_final) < EPS )) ) then

                    this%F(mu_ind)%FCmu(1, beta, alpha, ind_N2) = 8.0_dp

                    this%zero_f = this%zero_f + 1
                    this%dep_f = this%dep_f - 1

                end if final_zro_chk

            end if last_fc

        end do fc_cart_loop

        44 FORMAT("Force constant already set to zero")
        54 FORMAT("ERROR( in cart_dep ): dep_var can not be -1 or -7 for independent", F5.2)
        76 FORMAT("ERROR( in cart_dep ): p > fc", I4, I4)

        deallocate( GausMatT )

    end subroutine cart_dep


    subroutine Create_Index_ind(this, mu_ind, N2_ind)

        implicit none

        class(FCtyp)                                :: this

        integer, intent(in)                         :: mu_ind, N2_ind

        ! =================================== Local Variables =================================== !

        real(dp)                                    :: chk_indx1
        integer                                     :: val, PosIndx_flatSz2
        integer, dimension(2*ordfc2)                :: indindx
        integer, dimension(:,:), allocatable        :: tmpPosIndx_flat
        integer                                     :: alpha, beta

        ! =================================== Local Variables =================================== !

        alpha_loop: do alpha = 1, 3
            beta_loop: do beta = 1, 3

                chk_indx1 = this%F(mu_ind)%FCmu(1, beta, alpha, N2_ind)

                ind_chk: if ( dabs(chk_indx1-1.0_dp) < EPS ) then

                    val = this%IndFCPos(mu_ind)%Pos_mu(beta, alpha, N2_ind)

                    FILVAL_chk: if ( val /= FILLVAL ) then
                        write(*, 32) val
                        STOP

                    else FILVAL_chk

                        this%PosIndCntr = this%PosIndCntr + 1

                        ! ...........::::::::::::: For flat writing :::::::::::::........... !
                        indindx = (/mu_ind, N2_ind, alpha, beta/)

                        chk_Mat_size: if ( this%PosIndCntr <= this%Pos_flat_size ) then

                            this%PosIndx_flat(:, this%PosIndCntr) = indindx

                        else chk_Mat_size

                            PosIndx_flatSz2 = size( this%PosIndx_flat, 2 )

                            !call move_alloc(this%PosIndx_flat, tmpPosIndx_flat)
                            allocate( tmpPosIndx_flat( (2*ordfc2), PosIndx_flatSz2 ) )
                            tmpPosIndx_flat = this%PosIndx_flat
                            deallocate( this%PosIndx_flat )

                            ! ** !
                            this%Pos_flat_size = (this%Pos_flat_size + initMatsz2)
                            ! ** !

                            allocate( this%PosIndx_flat( (2*ordfc2), this%Pos_flat_size ) )
                            this%PosIndx_flat = 0
                            this%PosIndx_flat(:, 1:PosIndx_flatSz2) = tmpPosIndx_flat

                            this%PosIndx_flat(:, this%PosIndCntr) = indindx

                            !-! deallocate !-!
                            deallocate( tmpPosIndx_flat )
                            !-! deallocate !-!

                        end if chk_Mat_size
                        ! ...........::::::::::::: For flat writing :::::::::::::........... !

                        this%IndFCPos(mu_ind)%Pos_mu(beta, alpha, N2_ind) = this%PosIndCntr

                    end if FILVAL_chk

                end if ind_chk

            end do beta_loop
        end do alpha_loop

        32 FORMAT("ERROR( in Create_Index_ind ): Initially this%IndFCPos must be FILLVAL", I6)

    end subroutine Create_Index_ind

    include "InterMat_separate.f90"

    include "hdf5_wrap.f90"

    include "Reduce_separate.f90"

    include "ASR_separate.f90"

    include "ReducePntPerASR_separate.f90"

end module FCmod

