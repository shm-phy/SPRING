
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

    subroutine mesh_points(cell_vec, mesh, div_cell, ZeroCentered, mesh_cord)
    
        implicit none
    
        real(dp), dimension(3,3), intent(in)    :: cell_vec
        integer, intent(in)                     :: mesh(3)
        logical, intent(in)                     :: div_cell, ZeroCentered
        real(dp), allocatable, intent(out)      :: mesh_cord(:, :)
    
        ! ================================== Local variables ================================== !
        real(dp)                                :: cell_unit(3, 3)
        integer                                 :: i1, i2, i3, num, num_row
        integer                                 :: indx(3), fold(3)
        ! ================================== Local variables ================================== !
    
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
    
                    if ( ZeroCentered ) then
                        call fold_indx(mesh, indx, fold)
                    else
                        fold = indx
                    end if
    
                    mesh_cord(:, num) = matmul(cell_unit, dble(fold)) 
                    !mesh_cord(:, num) = matmul(cell_unit, fold)
                    num = num + 1
    
                end do
            end do
        end do
    
    end subroutine mesh_points

    
    subroutine q_points_highsymm(q_high_sym, G, qpoints, dist_arr, ZeroCenteredMesh) 
    
        implicit none
    
        real(dp), dimension(:, :), intent(in)                   :: q_high_sym
        real(dp), dimension(3, 3), intent(in)                   :: G
        real(dp), dimension(:, :), allocatable, intent(out)     :: qpoints
        real(dp), dimension(:), allocatable, intent(out)        :: dist_arr
        logical, dimension(:), allocatable, intent(out)         :: ZeroCenteredMesh
    
        ! ================================== Local variables ================================== !

        real(dp), dimension(3)                  :: start_node, end_node, v_hat
        real(dp)                                :: v_len, step_len, sd, dist
        real(dp), dimension(3)                  :: strtPoint, endPoint

        integer, dimension(:), allocatable      :: num_points
        integer                                 :: num_row, num_nodes, n, ii
    
        integer                                 :: sp, ep
        logical                                 :: IsZeroCentered

        ! ================================== Local variables ================================== !
    
        num_nodes = size(q_high_sym, 2)
    
        allocate(num_points(num_nodes))
        num_points(:) = int( q_high_sym(4, :) )
    
        num_row = sum(num_points)

        allocate(qpoints(3, num_row))
        allocate(dist_arr(num_row))

        allocate( ZeroCenteredMesh(num_row) )
        !ZeroCenteredMesh(:) = .true.
        ZeroCenteredMesh(:) = .false.
    
        sp = 1
        sd = 0.0
        dist = 0.0
    
        do n=1, (num_nodes-1)
    
            strtPoint = q_high_sym(1:3, n)
            endPoint = q_high_sym(1:3, n+1)

            if ( any(strtPoint < 0.0_dp) .or. any(endPoint < 0.0_dp) .or. &
               & all(dabs(strtPoint) <= 0.5_dp) .or. all(dabs(endPoint) <= 0.5_dp) ) then

                IsZeroCentered = .true.

            else

                IsZeroCentered = .false.

            end if

            start_node = matmul(G, strtPoint)
            end_node = matmul(G, endPoint)
    
            ep = sp + num_points(n+1)
    
            v_len = norm2(end_node - start_node)
            v_hat = (end_node - start_node) / v_len

            AvoidInf: if ( num_points(n+1) /= 0 ) then
                step_len = v_len / num_points(n+1)
            else AvoidInf
                step_len = 0.0_dp
            end if AvoidInf

            do ii=1, num_points(n+1)
    
                dist = sd + ii*step_len
                qpoints(1:3, sp+(ii-1)) = start_node + (ii*step_len) * v_hat
                dist_arr(sp+(ii-1)) = dist

                ZeroCenteredMesh(sp+(ii-1)) = IsZeroCentered
    
            end do
    
            sp = ep
            sd = dist
    
        end do
    
    end subroutine q_points_highsymm

end module CreateMesh

