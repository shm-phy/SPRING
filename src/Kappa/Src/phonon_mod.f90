
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

module phonon_m

    use kinds,          only : dp
    use Helper,         only : DistributeQpoints
    use DynaMat,        only : get_phonon, get_phonon_pd
    use unit_cell,      only : cell
    use FC2_mod,        only : FC2type
    use EwaldMod,       only : EwaldParam

    implicit none
    private

    type, public        :: Phon

        complex(dp), dimension(:, :, :), allocatable    :: Evec
        complex(dp), dimension(:, :, :, :), allocatable :: SpcProjOp !* Spectral projection operator *!
        real(dp), dimension(:, :), allocatable          :: omega, OnebyOmega
        real(dp), dimension(:, :), allocatable          :: nBE
        real(dp), dimension(:, :, :), allocatable       :: grp_vel
        real(dp)                                        :: omega_min, omega_max

    contains
        procedure, public, pass                         :: set_Phon, set_Phon_mp, &
                                                         & clean_Phon
        !FINAL                                           :: clean_Phon
    end type Phon
    
contains

    subroutine set_Phon(this, sys, FC2, EwaldConst, qpoints, T, freq_cut, Ndof, LongEW, flag_pd)

        implicit none

        class(Phon)                                                 :: this
        type(cell), intent(in)                                      :: sys
        type(FC2type), intent(in)                                   :: FC2
        type(EwaldParam), intent(in)                                :: EwaldConst

        real(dp), dimension(:, :), intent(in)                       :: qpoints
        real(dp), intent(in)                                        :: T, freq_cut

        integer, intent(in)                                         :: Ndof

        logical, intent(in)                                         :: LongEW, flag_pd

        ! ================================ Local variables ================================ !

        complex(dp), dimension(:, :, :), allocatable    :: Evec
        complex(dp), dimension(:, :, :, :), allocatable :: ProjOp
        real(dp), dimension(:, :), allocatable          :: freq, OnebyFreq, nBE
        real(dp), dimension(:, :, :), allocatable       :: grp_vel
        real(dp)                                        :: freq_min, freq_max

        integer                                         :: Nq, myQsize, myOffset

        ! ================================ Local variables ================================ !

        Nq = size(qpoints, 2)
        myQsize = Nq
        myOffset = 0

        allocate( freq(Ndof, Nq) )
        freq = 0.0_dp

        allocate( Evec(Ndof, Ndof, Nq) )
        Evec = dcmplx(0.0_dp, 0.0_dp)

        if ( .not. flag_pd ) then

            allocate( OnebyFreq(Ndof, Nq) )
            OnebyFreq = 0.0_dp
            allocate( nBE(Ndof, Nq) )
            nBE = 0.0_dp
            allocate( grp_vel(3, Ndof, Nq) )
            grp_vel = 0.0_dp

            call get_phonon( sys, FC2, EwaldConst, Nq, Ndof, myQsize, myOffset, &
                           & qpoints, T, LongEw, Evec, freq, OnebyFreq, nBE, grp_vel, &
                           & freq_min, freq_max )

        else

            allocate( ProjOp(Ndof, Ndof, Ndof, Nq) )
            ProjOp = dcmplx( 0.0_dp, 0.0_dp )

            call get_phonon_pd( sys, FC2, EwaldConst, Nq, Ndof, myQsize, myOffset, &
                              & qpoints, LongEw, freq, ProjOp, Evec, freq_min, freq_max )

        end if

        !call move_alloc( freq, this%omega )
        allocate( this%omega(Ndof, Nq) )
        this%omega = freq

        !call move_alloc( Evec, this%Evec )
        allocate( this%Evec(Ndof, Ndof, Nq) )
        this%Evec = Evec
        
        if ( .not. flag_pd ) then
            
            !call move_alloc( OnebyFreq, this%OnebyOmega )
            allocate( this%OnebyOmega(Ndof, Nq) )
            this%OnebyOmega = OnebyFreq

            !call move_alloc( nBE, this%nBE )
            allocate( this%nBE(Ndof, Nq) )
            this%nBE = nBE

            !call move_alloc( grp_vel, this%grp_vel )
            allocate( this%grp_vel(3, Ndof, Nq) )
            this%grp_vel = grp_vel

        else

            !call move_alloc( ProjOp, this%SpcProjOp )
            allocate( this%SpcProjOp(Ndof, Ndof, Ndof, Nq) )
            this%SpcProjOp = ProjOp

        end if

        this%omega_min = freq_min
        this%omega_max = freq_max
        this%omega_min = freq_cut !**!

        deallocate( freq )
        deallocate( Evec )
        if ( .not. flag_pd ) then
            deallocate( OnebyFreq, nBE, grp_vel )
        else
            deallocate( ProjOp )
        end if

    end subroutine set_Phon


    subroutine set_Phon_mp(this, sys, FC2, EwaldConst, qpoints, T, freq_cut, Ndof, LongEW, flag_pd)

        implicit none

        class(Phon)                                                 :: this
        type(cell), intent(in)                                      :: sys
        type(FC2type), intent(in)                                   :: FC2
        type(EwaldParam), intent(in)                                :: EwaldConst

        real(dp), dimension(:, :), intent(in)                       :: qpoints
        real(dp), intent(in)                                        :: T, freq_cut

        integer, intent(in)                                         :: Ndof

        logical, intent(in)                                         :: LongEW, flag_pd

        ! ================================ Local variables ================================ !

        character(len=128)                                              :: msg

        complex(dp), allocatable, dimension(:, :, :), codimension[:]    :: Evec
        complex(dp), allocatable, dimension(:, :, :, :), codimension[:] :: ProjOp

        real(dp), allocatable, dimension(:, :), codimension[:]          :: freq, OnebyFreq, nBE
        real(dp), allocatable, dimension(:, :, :), codimension[:]       :: grp_vel
        real(dp)                                                        :: freq_min, freq_max

        integer                                                         :: Nq, myQsize, myOffset, &
                                                                         & istat
        integer, dimension(2)                                           :: my_edge

        ! ================================ Local variables ================================ !

        Nq = size(qpoints, 2)
        call DistributeQpoints(Nq, myQsize, myOffset, my_edge)
        SYNC ALL

        allocate( freq(Ndof, Nq)[*] )
        SYNC ALL
        freq = 0.0_dp

        allocate( Evec(Ndof, Ndof, Nq)[*] )
        SYNC ALL
        Evec = dcmplx(0.0_dp, 0.0_dp)

        PointDef1: if ( .not. flag_pd ) then

            allocate( OnebyFreq(Ndof, Nq)[*] )
            SYNC ALL
            OnebyFreq = 0.0_dp

            allocate( nBE(Ndof, Nq)[*] )
            SYNC ALL
            nBE = 0.0_dp

            allocate( grp_vel(3, Ndof, Nq)[*] )
            SYNC ALL
            grp_vel = 0.0_dp

            call get_phonon( sys, FC2, EwaldConst, Nq, Ndof, myQsize, myOffset, &
                           & qpoints, T, LongEw, Evec, freq, OnebyFreq, nBE, grp_vel, &
                           & freq_min, freq_max )
        
        else PointDef1

            allocate( ProjOp(Ndof, Ndof, Ndof, Nq)[*] )
            SYNC ALL
            ProjOp = dcmplx(0.0_dp, 0.0_dp)

            call get_phonon_pd( sys, FC2, EwaldConst, Nq, Ndof, myQsize, myOffset, &
                              & qpoints, LongEw, freq, ProjOp, Evec, freq_min, freq_max )

        end if PointDef1

        SYNC ALL
        call co_sum( freq, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

        SYNC ALL
        call co_sum( Evec, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

        PointDef2: if ( .not. flag_pd ) then

            SYNC ALL
            call co_sum( OnebyFreq, stat=istat, ERRMSG=msg )
            if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

            SYNC ALL
            call co_sum( nBE, stat=istat, ERRMSG=msg )
            if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

            SYNC ALL
            call co_sum( grp_vel, stat=istat, ERRMSG=msg )
            if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

        else PointDef2

            SYNC ALL
            call co_sum( ProjOp, stat=istat, ERRMSG=msg )
            if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

        end if PointDef2

        SYNC ALL
        call co_broadcast( freq_min, source_image=1, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_broadcast : ', A128 )") msg

        SYNC ALL
        call co_broadcast( freq_max, source_image=1, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_broadcast : ', A128 )") msg
        SYNC ALL

        !call move_alloc( freq, this%omega )
        allocate( this%omega(Ndof, Nq) )
        this%omega = freq

        !call move_alloc( Evec, this%Evec )
        allocate( this%Evec(Ndof, Ndof, Nq) )
        this%Evec = Evec

        PointDef3: if ( .not. flag_pd ) then

            !call move_alloc( OnebyFreq, this%OnebyOmega )
            allocate( this%OnebyOmega(Ndof, Nq) )
            this%OnebyOmega = OnebyFreq

            !call move_alloc( nBE, this%nBE )
            allocate( this%nBE(Ndof, Nq) )
            this%nBE = nBE

            !call move_alloc( grp_vel, this%grp_vel )
            allocate( this%grp_vel(3, Ndof, Nq) )
            this%grp_vel = grp_vel

        else PointDef3

            !call move_alloc( ProjOp, this%SpcProjOp )
            allocate( this%SpcProjOp(Ndof, Ndof, Ndof, Nq) )
            this%SpcProjOp = ProjOp

        end if PointDef3

        this%omega_min = freq_min
        this%omega_max = freq_max

        this%omega_min = freq_cut !**!

        SYNC ALL

        deallocate( freq )
        deallocate( Evec )

        PointDef4: if ( .not. flag_pd ) then
            deallocate( OnebyFreq, nBE, grp_vel )
        else PointDef4
            deallocate( ProjOp )
        end if PointDef4

    end subroutine set_Phon_mp


    subroutine clean_Phon( this, flag_pd )

        ! ** Finalizer ** ! (Not the standard Finalizer)
        implicit none

        CLASS(Phon)                     :: this !** TYPE(Phon)
        logical, intent(in)             :: flag_pd

        ! ================================ Local variables ================================ !

        integer                         :: istat

        ! ================================ Local variables ================================ !

        deallocate( this%omega, STAT=istat )
        deallocate( this%Evec, STAT=istat )

        if ( .not. flag_pd ) then
            deallocate( this%grp_vel, STAT=istat )
            deallocate( this%OnebyOmega, STAT=istat )
            deallocate( this%nBE, STAT=istat )

        else
            deallocate( this%SpcProjOp, STAT=istat )

        end if

    end subroutine clean_Phon

end module phonon_m

