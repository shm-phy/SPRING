
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

    public  :: mesh_points, q_points_highsymm

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


subroutine q_points_highsymm(q_high_sym, G, qpoints) 

    implicit none

    real(dp), dimension(:, :), intent(in)                   :: q_high_sym
    real(dp), dimension(3, 3), intent(in)                   :: G
    real(dp), dimension(:, :), allocatable, intent(out)     :: qpoints

    integer, dimension(:), allocatable      :: num_points
    integer                                 :: num_row, num_nodes, n, ii

    real(dp), dimension(3)                  :: start_node, end_node, v_hat
    real(dp)                                :: v_len, step_len, sd, dist
    integer                                 :: sp, ep

    num_nodes = size(q_high_sym, 2)

    allocate(num_points(num_nodes))
    num_points(:) = int( q_high_sym(4, :) )

    num_row = sum(num_points)
    allocate(qpoints(4, num_row))

    sp = 1
    sd = 0.0
    dist = 0.0

    do n=1, (num_nodes-1)

        start_node = matmul(G, q_high_sym(1:3, n))
        end_node = matmul(G, q_high_sym(1:3, n+1))

        ep = sp + num_points(n+1)

        v_len = norm2(end_node-start_node)
        v_hat = (end_node-start_node) / v_len
        step_len = v_len / num_points(n+1)

        do ii=1, num_points(n+1)

            dist = sd + ii*step_len
            qpoints(1, sp+(ii-1)) = dist
            qpoints(2:4, sp+(ii-1)) = start_node + (ii*step_len) * v_hat

        end do

        sp = ep
        sd = dist

    end do

end subroutine q_points_highsymm

end module CreateMesh

