
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


module CreateMesh

    use kinds,          only: dp
    use IndxChng,       only: fold_indx
    implicit none
    private

    public  :: mesh_points

contains

    subroutine mesh_points(cell_vec, mesh, div_cell, mesh_cord)
    
        implicit none
    
        real(dp), dimension(3,3), intent(in)    :: cell_vec
        integer, intent(in)                     :: mesh(3)
        logical, intent(in)                     :: div_cell
        real(dp), allocatable, intent(out)      :: mesh_cord(:, :)
    
        real(dp)                                :: cell_unit(3, 3)
        integer                                 :: i1, i2, i3, num, num_row
        integer                                 :: indx(3), fold(3)
    
        chk_div: if ( div_cell ) then
    
            do i1 = 1, 3
                cell_unit(:, i1) = cell_vec(:, i1) / dble(mesh(i1))
            end do
    
        else chk_div
    
            cell_unit = cell_vec
    
        end if chk_div
    
        num_row = product(mesh)
        allocate(mesh_cord(3, num_row))
        
        num = 1
        do i1 = 0, (mesh(1)-1)
            do i2 = 0, (mesh(2)-1)
                do i3 = 0, (mesh(3)-1)
    
                    indx = (/i1, i2, i3/)
                    call fold_indx(mesh, indx, fold)
    
                    !mesh_cord(:, num) = matmul(cell_unit, dble(fold)) 
                    mesh_cord(:, num) = matmul(cell_unit, fold)
                    num = num + 1
    
                end do
            end do
        end do
    
    end subroutine mesh_points

end module CreateMesh

