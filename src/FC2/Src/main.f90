
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
    use Symmetry,                   only : PntPer_typ
    use CutInfo,                    only : CutOff
    use FCmod,                      only : FCtyp

    implicit none

    type(cell)                                  :: sys
    character(len=256)                          :: filename
    real(dp)                                    :: ChnkSz, ZeroEPS
    logical                                     :: AdvncOut

    type(PntPer_typ)                            :: PntPer
    real(dp)                                    :: symprec

    type(CutOff)                                :: Cut_Off
    integer, dimension(:), allocatable          :: near_neigh

    integer                                     :: Nbasis
    type(FCtyp)                                 :: FCObj

    ! ------------------------------------------------------------------------------------- !

    call ShowWelcomeBanner()

    call get_command_line(filename, ChnkSz, ZeroEPS, AdvncOut)

    call sys%init_data(filename)
    Nbasis = sys%natm

    FC4chk: if ( sys%SecondFC ) then

        symprec = 1.0E-7_dp

        call PntPer%set_PntPerSPGlib(sys, symprec)

        allocate( near_neigh(sys%natm) )
        near_neigh(:) = sys%nn_2nd

        call Cut_Off%set_CutOff(sys, near_neigh)

        call FCObj%init_FCtyp(ChnkSz, AdvncOut, sys, Cut_Off, PntPer)
        call FCObj%Find_Independent(sys, Cut_Off, PntPer)

        call FCObj%ApplyASR(Cut_Off, Nbasis)

        call FCObj%ReducePntPerASR(Nbasis, ZeroEPS, Cut_Off)

        call FCObj%Saveh5_F(ZeroEPS, sys, Cut_Off)

    else FC4chk

        write(*, 100)
        STOP

    end if FC4chk

    write(*, '(/)')
    
    ! ------------------------------------------------------------------------------------- !

    100 FORMAT("No input detected for Second-Order Force Constant calculation ( SecondFC = .false. )")

end program main

