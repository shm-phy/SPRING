
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
    use DistributeQs,   only : DistributeQpoints
    use DynaMat,        only : get_phonon
    use unit_cell,      only : cell
    use FC2_mod,        only : FC2type
    use EwaldMod,       only : EwaldParam

    implicit none
    private

    type, public        :: Phon

        complex(dp), dimension(:, :, :), allocatable    :: Evec
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

    subroutine set_Phon(this, sys, FC2, EwaldConst, qpoints, T, Ndof, LongEW)

        implicit none

        class(Phon)                                                 :: this
        type(cell), intent(in)                                      :: sys
        type(FC2type), intent(in)                                   :: FC2
        type(EwaldParam), intent(in)                                :: EwaldConst

        real(dp), dimension(:, :), intent(in)                       :: qpoints
        real(dp), intent(in)                                        :: T

        integer, intent(in)                                         :: Ndof

        logical, intent(in)                                         :: LongEW

        ! ================================ Local variables ================================ !

        complex(dp), dimension(:, :, :), allocatable   :: Evec
        real(dp), dimension(:, :), allocatable         :: freq, OnebyFreq, nBE
        real(dp), dimension(:, :, :), allocatable      :: grp_vel
        real(dp)                                       :: freq_min, freq_max

        integer                                         :: Nq, myQsize, myOffset

        ! ================================ Local variables ================================ !

        Nq = size(qpoints, 2)
        myQsize = Nq
        myOffset = 0

        allocate( Evec(Ndof, Ndof, Nq) )
        Evec = dcmplx(0.0_dp, 0.0_dp)
        allocate( freq(Ndof, Nq) )
        freq = 0.0_dp
        allocate( OnebyFreq(Ndof, Nq) )
        OnebyFreq = 0.0_dp
        allocate( nBE(Ndof, Nq) )
        nBE = 0.0_dp
        allocate( grp_vel(3, Ndof, Nq) )
        grp_vel = 0.0_dp

        call get_phonon( sys, FC2, EwaldConst, Nq, Ndof, myQsize, myOffset, &
                       & qpoints, T, LongEw, Evec, freq, OnebyFreq, nBE, grp_vel, &
                       & freq_min, freq_max )

        !call move_alloc( Evec, this%Evec )
        allocate( this%Evec(Ndof, Ndof, Nq) )
        this%Evec = Evec
        
        !call move_alloc( freq, this%omega )
        allocate( this%omega(Ndof, Nq) )
        this%omega = freq

        !call move_alloc( OnebyFreq, this%OnebyOmega )
        allocate( this%OnebyOmega(Ndof, Nq) )
        this%OnebyOmega = OnebyFreq

        !call move_alloc( nBE, this%nBE )
        allocate( this%nBE(Ndof, Nq) )
        this%nBE = nBE

        !call move_alloc( grp_vel, this%grp_vel )
        allocate( this%grp_vel(3, Ndof, Nq) )
        this%grp_vel = grp_vel

        this%omega_min = freq_min
        this%omega_max = freq_max

        deallocate( Evec )
        deallocate( freq, OnebyFreq, nBE, grp_vel )

    end subroutine set_Phon


    subroutine set_Phon_mp(this, sys, FC2, EwaldConst, qpoints, T, Ndof, LongEW)

        implicit none

        class(Phon)                                                 :: this
        type(cell), intent(in)                                      :: sys
        type(FC2type), intent(in)                                   :: FC2
        type(EwaldParam), intent(in)                                :: EwaldConst

        real(dp), dimension(:, :), intent(in)                       :: qpoints
        real(dp), intent(in)                                        :: T

        integer, intent(in)                                         :: Ndof

        logical, intent(in)                                         :: LongEW

        ! ================================ Local variables ================================ !

        character(len=128)                                              :: msg

        complex(dp), allocatable, dimension(:, :, :), codimension[:]    :: Evec

        real(dp), allocatable, dimension(:, :), codimension[:]          :: freq, OnebyFreq, nBE
        real(dp), allocatable, dimension(:, :, :), codimension[:]       :: grp_vel
        real(dp)                                                        :: freq_min, freq_max

        integer                                                         :: Nq, myQsize, myOffset, &
                                                                         & istat
        integer, dimension(2)                                           :: my_edge

        ! ================================ Local variables ================================ !

        Nq = size(qpoints, 2)

        allocate( Evec(Ndof, Ndof, Nq)[*] )
        SYNC ALL
        Evec = dcmplx(0.0_dp, 0.0_dp)

        allocate( freq(Ndof, Nq)[*] )
        SYNC ALL
        freq = 0.0_dp

        allocate( OnebyFreq(Ndof, Nq)[*] )
        SYNC ALL
        OnebyFreq = 0.0_dp

        allocate( nBE(Ndof, Nq)[*] )
        SYNC ALL
        nBE = 0.0_dp

        allocate( grp_vel(3, Ndof, Nq)[*] )
        SYNC ALL
        grp_vel = 0.0_dp

        call DistributeQpoints(Nq, myQsize, myOffset, my_edge)
        SYNC ALL

        call get_phonon( sys, FC2, EwaldConst, Nq, Ndof, myQsize, myOffset, &
                       & qpoints, T, LongEw, Evec, freq, OnebyFreq, nBE, grp_vel, &
                       & freq_min, freq_max )

        SYNC ALL
        call co_sum( Evec, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

        SYNC ALL
        call co_sum( freq, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

        SYNC ALL
        call co_sum( OnebyFreq, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

        SYNC ALL
        call co_sum( nBE, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

        SYNC ALL
        call co_sum( grp_vel, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_sum : ', A128 )") msg

        SYNC ALL
        call co_broadcast( freq_min, source_image=1, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_broadcast : ', A128 )") msg

        SYNC ALL
        call co_broadcast( freq_max, source_image=1, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_broadcast : ', A128 )") msg
        SYNC ALL

        !call move_alloc( Evec, this%Evec )
        allocate( this%Evec(Ndof, Ndof, Nq) )
        this%Evec = Evec
        
        !call move_alloc( freq, this%omega )
        allocate( this%omega(Ndof, Nq) )
        this%omega = freq

        !call move_alloc( OnebyFreq, this%OnebyOmega )
        allocate( this%OnebyOmega(Ndof, Nq) )
        this%OnebyOmega = OnebyFreq

        !call move_alloc( nBE, this%nBE )
        allocate( this%nBE(Ndof, Nq) )
        this%nBE = nBE

        !call move_alloc( grp_vel, this%grp_vel )
        allocate( this%grp_vel(3, Ndof, Nq) )
        this%grp_vel = grp_vel

        this%omega_min = freq_min
        this%omega_max = freq_max

        SYNC ALL
        deallocate( Evec )
        deallocate( freq, OnebyFreq, nBE, grp_vel )

    end subroutine set_Phon_mp


    subroutine clean_Phon( this )

        ! ** Finalizer ** ! (Not the standard Finalizer)
        implicit none

        CLASS(Phon)                     :: this !** TYPE(Phon)

        ! ================================ Local variables ================================ !

        integer                         :: istat

        ! ================================ Local variables ================================ !

        deallocate( this%omega, this%nBE, STAT=istat )
        deallocate( this%Evec, STAT=istat )
        deallocate( this%grp_vel, STAT=istat )

    end subroutine clean_Phon

end module phonon_m

