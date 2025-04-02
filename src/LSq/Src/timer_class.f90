
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


! A timer Class from "Fortran for Scientists and Engineers" by Stephen J. Chapman !

module timer_class

    use kinds,      only : dp

    implicit none
    private

    type, public    :: timer

        PRIVATE
        real(dp)    :: saved_time

    contains

        procedure, public, pass   :: start_timer => start_timer_sub
        procedure, public, pass   :: elapsed_time => elapsed_time_fn
        
    end type timer

    private     :: start_timer_sub, elapsed_time_fn

contains

    subroutine start_timer_sub(this)

        implicit none

        class(timer)            :: this

        ! ========================= Local Variables ========================= !

        integer, dimension(8)   :: val

        ! ========================= Local Variables ========================= !

        call date_and_time( VALUES=val )

        this%saved_time = 86400.0_dp * dble(val(3)) + 3600.0_dp * dble(val(5)) + &
                        & 60.0_dp * dble(val(6)) + dble(val(7)) + 0.001_dp * dble(val(8))

    end subroutine start_timer_sub


    Function elapsed_time_fn(this) Result(elps_time)

        implicit none

        class(timer)            :: this

        real(dp)                :: elps_time !Result

        ! ========================= Local Variables ========================= !

        integer, dimension(8)   :: val
        real(dp)                :: current_time

        ! ========================= Local Variables ========================= !

        call date_and_time( VALUES=val )

        current_time = 86400.0_dp * dble(val(3)) + 3600.0_dp * dble(val(5)) + &
                     & 60.0_dp * dble(val(6)) + dble(val(7)) + 0.001_dp * dble(val(8))

        elps_time = current_time - this%saved_time

    end Function elapsed_time_fn

end module timer_class

