
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


module FindForce

    use kinds,              only : dp
    use hdf5_wrap,          only : r_force_disp
    use IndxChng,           only : indx2num, PeriodicMap
    use constants,          only : cmp_prec
    
    implicit none
    private

    public  :: CalculateForce 

contains

    subroutine CalculateForce(fdfile, FC2_dd, RCb, force_dd, force_md)
    
        implicit none
    
        character(len=*)                                                :: fdfile
        real(dp), dimension(:,:,:,:), contiguous, intent(in)            :: FC2_dd
        integer, dimension(:,:,:), contiguous, intent(in)               :: RCb
        real(dp), dimension(:), allocatable, intent(out)                :: force_dd, force_md
    
        real(dp), dimension(:,:,:,:,:,:), allocatable   :: disp, force
        integer, dimension(3)                           :: SupDim
        integer                                         :: Nbasis, ts, NCb
    
        integer, dimension(6)                           :: bound, indx_strt, indx_end
        integer                                         :: strtN, endN
    
        integer, dimension(4)                           :: Cellnu
        integer                                         :: mu, alpha, N2, beta, tsSupNum
    
        real(dp), dimension(:), allocatable             :: disp_N2_beta, force_mu_alpha, md_frc_mu_alpha
        real(dp)                                        :: phi2_dd
    
        call r_force_disp(fdfile, disp, force)
    
        SupDim = (/size(disp, 6), size(disp, 5), size(disp, 4)/) 
        Nbasis = size(disp, 3)
        ts = size(disp, 1)
        NCb = size(RCb, 2)
    
        bound(1) = Nbasis
        bound(2) = 3
        bound(3) = ts
        bound(4:6) = SupDim
    
        indx_strt = 0
        indx_strt(3:6) = 1
    
        indx_end = 0
        indx_end(3:6) = bound(3:6)
    
        tsSupNum = ts * product(SupDim)
    
        allocate(force_dd(product(bound)))
        force_dd = 0.0_dp
    
        allocate(force_md(product(bound)))
        force_md = 0.0_dp
    
        mu_loop: do mu = 1, Nbasis
            alpha_loop: do alpha = 1, 3
            
                allocate(force_mu_alpha(tsSupNum))
                force_mu_alpha = 0.0_dp
        
                indx_strt(1) = mu
                indx_strt(2) = alpha
                strtN = indx2num(bound, indx_strt)
        
                indx_end(1) = mu
                indx_end(2) = alpha
                endN = indx2num(bound, indx_end)
        
                N2_loop: do N2 = 1, NCb
        
                    Cellnu = RCb(:, N2, mu)
                    beta_loop: do beta = 1, 3
                         
                        phi2_dd = FC2_dd(beta, alpha, N2, mu)
                        call GetDisplacements(disp, Cellnu, beta, SupDim, ts, disp_N2_beta)
        
                        force_mu_alpha = force_mu_alpha + (disp_N2_beta * phi2_dd)
        
                        deallocate(disp_N2_beta)
        
                    end do beta_loop
                end do N2_loop
        
                call GetMDForce(force, mu, alpha, SupDim, ts, md_frc_mu_alpha)
        
                zero_chk: if ( any(dabs(force_dd(strtN:endn)) > cmp_prec) .or. &
                             & any(dabs(force_md(strtN:endn)) > cmp_prec)) then
        
                    write(*, *) 'Error: Initially all the elements are not zero'
        
                end if zero_chk
        
                force_dd(strtN:endN) = force_mu_alpha
                force_md(strtN:endN) = md_frc_mu_alpha
        
                deallocate(force_mu_alpha)
                deallocate(md_frc_mu_alpha)
        
            end do alpha_loop
        end do mu_loop

        deallocate( disp, force )
    
    end subroutine CalculateForce
    
    
    subroutine GetDisplacements(disp, Cellnu, beta, SupDim, ts, disp_N2_beta)
    
        implicit none
    
        real(dp), dimension(:,:,:,:,:,:), contiguous, intent(in)        :: disp 
        integer, dimension(4), intent(in)                               :: Cellnu
        integer, intent(in)                                             :: beta
        integer, dimension(3), intent(in)                               :: SupDim
        integer, intent(in)                                             :: ts
        real(dp), dimension(:), allocatable, intent(out)                :: disp_N2_beta
    
        integer, dimension(3)       :: Cell, Center, CellFold
        integer                     :: nu
        integer                     :: tsSupN, num
        integer                     :: tstep, cx, cy, cz
    
        Cell = Cellnu(1:3)
        nu = Cellnu(4)
    
        tsSupN = ts * product(SupDim)
        allocate(disp_N2_beta(tsSupN))
        disp_N2_beta = 0.0_dp
    
        num = 1
    
        cx_loop: do cx = 1, SupDim(1)
            cy_loop: do cy = 1, SupDim(2)
                cz_loop: do cz = 1, SupDim(3)
                    Center = (/cx, cy, cz/)
                    CellFold = PeriodicMap(SupDim, Center, Cell)
        
                    ts_loop: do tstep = 1, ts
                        disp_N2_beta(num) = &
                        & -1.0_dp * disp(tstep, beta, nu, CellFold(3), CellFold(2), CellFold(1))
                        num = num + 1
                    end do ts_loop
        
                end do cz_loop
            end do cy_loop
        end do cx_loop
    
        chk_num: if ( (num-1) /= tsSupN ) then
            write(*, *) "Error: All the timesteps and cells are not covered"
        end if chk_num
    
    end subroutine GetDisplacements
    

    subroutine GetMDForce(force, mu, alpha, SupDim, ts, md_frc_mu_alpha)
    
        implicit none
    
        real(dp), dimension(:,:,:,:,:,:), contiguous, intent(in)        :: force
        integer, intent(in)                                             :: mu
        integer, intent(in)                                             :: alpha
        integer, dimension(3), intent(in)                               :: SupDim
        integer, intent(in)                                             :: ts
        real(dp), dimension(:), allocatable, intent(out)                :: md_frc_mu_alpha
    
        integer                 :: tsSupN, num
        integer                 :: tstep, cx, cy, cz
    
        tsSupN = ts * product(SupDim)
        allocate(md_frc_mu_alpha(tsSupN))
        md_frc_mu_alpha = 0.0_dp
    
        num = 1
    
        cx_loop: do cx = 1, SupDim(1)
            cy_loop: do cy = 1, SupDim(2)
                cz_loop: do cz = 1, SupDim(3)
        
                    ts_loop: do tstep = 1, ts
                        md_frc_mu_alpha(num) = force(tstep, alpha, mu, cz, cy, cx)
                        num = num + 1
                    end do ts_loop
        
                end do cz_loop
            end do cy_loop
        end do cx_loop
    
        chk_num: if ( (num-1) /= tsSupN ) then
            write(*, *) "Error: All the timesteps and cells are not covered"
        end if chk_num
    
    end subroutine GetMDForce

end module FindForce

