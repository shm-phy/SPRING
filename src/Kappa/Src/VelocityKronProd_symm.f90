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

module SymmVelKronProd

    use kinds,              only : dp
    use Helper,             only : kron_prod
    use Symmetry,           only : SpaceGr_typ
    use TetrahedronQ,       only : q_tetrahedron
    use phonon_m,           only : Phon

    implicit none
    private
    public                      :: GrVelKronProdSymm

contains 

    subroutine GrVelKronProdSymm( Nq0points, Ndof, my_Qsizeq0, my_offsetq0, &
                                & SpGr, tetraQ, ph_q0, KrnGrpVelSym )

        implicit none

        integer, intent(in)                                                     :: Nq0points, Ndof, &
                                                                                 & my_Qsizeq0, my_offsetq0

        type(SpaceGr_typ), intent(in)                                           :: SpGr
        type(q_tetrahedron), intent(in)                                         :: tetraQ
        type(Phon), intent(in)                                                  :: ph_q0

        real(dp), allocatable, dimension(:,:,:), codimension[:], intent(inout)  :: KrnGrpVelSym

        ! ================================== Local Variables ================================== !

        character(len=128)                  :: msg
        real(dp)                            :: p_q0
        real(dp), dimension(3,3)            :: R_pnt
        real(dp), dimension(9)              :: RvgKronRvg
        real(dp), dimension(3)              :: vg, Rvg
        integer                             :: NumOfSym, ns, istat, ii, q0i, s

        ! ================================== Local Variables ================================== !

        allocate( KrnGrpVelSym(9, Ndof, Nq0points)[*], STAT=istat, ERRMSG=msg )
        if ( istat /= 0 ) then
            write(*, 24) this_image(), msg
            24 FORMAT( " Memory allocation ERROR in image : ", I5, ". Error message: ", A128 )
            ERROR STOP
        end if
        KrnGrpVelSym = 0.0_dp
        SYNC ALL

        NumOfSym = SpGr%Nsym

        q0loop: do ii = 1, my_Qsizeq0
            q0i = my_offsetq0 + ii

            p_q0 = 1.0_dp / dble( NumOfSym / tetraQ%multiplicity(q0i) )

            dof_loop: do s = 1, Ndof
                vg = ph_q0%grp_vel(:, s, q0i)

                PointSymmLoop: do ns = 1, NumOfSym
                    R_pnt = SpGr%Rxyz(:, :, ns)
                    Rvg = matmul( R_pnt, vg )
                    RvgKronRvg = kron_prod( Rvg, Rvg )

                    KrnGrpVelSym(:, s, q0i) = KrnGrpVelSym(:, s, q0i) + RvgKronRvg

                end do PointSymmLoop

                KrnGrpVelSym(:, s, q0i) = KrnGrpVelSym(:, s, q0i) * p_q0

            end do dof_loop

        end do q0loop

    end subroutine GrVelKronProdSymm

end module SymmVelKronProd

