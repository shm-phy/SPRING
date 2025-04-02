
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


module EwaldMod

    use     kinds,      only : dp
    use constants,      only : PI
    use unit_cell,      only : cell
    use mklWrap,        only : dinverse, ddet
    use CreateMesh,     only : mesh_points
    use EwaldSum,       only : RealSpaceEW, KSpaceEW, CombineZStar

    implicit none
    private

    type, public :: EwaldParam

        logical                                         :: decide_EwParam
        real(dp)                                        :: prec_Ew, Lmb
        integer, dimension(3)                           :: Gmesh, Rmesh

        real(dp), dimension(3,3)                        :: eps_inv
        real(dp)                                        :: det

        real(dp), allocatable, dimension(:, :)          :: Rpoints, Gpoints
        complex(dp), allocatable, dimension(:,:,:,:)    :: CConst
        complex(dp), allocatable, dimension(:,:,:)      :: CSumq0

        complex(dp), allocatable, dimension(:,:,:)      :: DynEwq

    contains
        procedure, public, pass         :: set_EwaldParam, set_DynEw
        procedure, private, nopass      :: DecideEwaldParameters, ConstTerm, SumRule

    end type EwaldParam

contains


    subroutine set_EwaldParam(this, sys, filename)
        
        implicit none

        type(cell), intent(in)                          :: sys
        character(len = *), intent(in)                  :: filename

        class(EwaldParam)                               :: this

        logical                                         :: decide_EwParam=.true.
        real(dp)                                        :: prec_Ew=1.0E-6_dp, Lmb = 0.80_dp
        integer, dimension(3)                           :: Gmesh, Rmesh

        namelist        /EwaldInfo/      decide_EwParam, prec_Ew, Lmb, Gmesh, Rmesh

        real(dp)                                        :: Ew_alpha
        integer                                         :: NKcut
        integer, dimension(3)                           :: Ew_Gmesh

        real(dp), dimension(3,3)                        :: eps_inv
        real(dp)                                        :: det

        real(dp), allocatable, dimension(:, :)          :: Rpoints, Gpoints
        complex(dp), allocatable, dimension(:,:,:,:)    :: CConst
        complex(dp), allocatable, dimension(:,:,:)      :: CSumq0

        integer                     :: err
        character(len=512)          :: err_msg


        open(unit=5, file=filename, status='OLD', iostat=err, iomsg=err_msg, &
             action='READ', delim='APOSTROPHE')

            open_chk: if ( err /= 0 ) then
                write(*, *) 'Input file OPEN failed: iostat = ', err
                write(*, *) 'Error message = ', err_msg
            end if open_chk

            read(unit=5, nml=EwaldInfo)

        close(unit=5, status='KEEP', iostat=err, iomsg=err_msg)

        eps_inv = sys%eps
        eps_inv = transpose(eps_inv)
        call dinverse(eps_inv)
        eps_inv = transpose(eps_inv)
        call ddet(sys%eps, det)

        chk_EwPrm: if ( decide_EwParam ) then

            call DecideEwaldParameters(sys%eps, eps_inv, sys%latvec, sys%G, sys%sup_cell, &
                                     & prec_Ew, Ew_alpha, NKcut)
            Ew_Gmesh(:) = 2*NKcut + 1

        else chk_EwPrm
            Ew_alpha = Lmb
            Ew_Gmesh = Gmesh

        end if chk_EwPrm

        this%decide_EwParam = decide_EwParam
        this%prec_Ew = prec_Ew

        this%Lmb = Ew_alpha
        this%Gmesh = Ew_Gmesh
        this%Rmesh = Rmesh

        this%eps_inv = eps_inv
        this%det = det

        call mesh_points(sys%latvec, this%Rmesh, .false., Rpoints)
        call mesh_points(sys%G, this%Gmesh, .false., Gpoints)

        !allocate( this%Rpoints(shape(Rpoints)) )
        !allocate( this%Gpoints(shape(Gpoints)) )
        this%Rpoints = Rpoints
        this%Gpoints = Gpoints

        call ConstTerm(this%eps_inv, this%det, this%Lmb, sys%natm, CConst)

        !allocate( this%CConst(shape(CConst)) )
        this%CConst = CConst

        call SumRule(sys, this%Rpoints, this%Gpoints, this%eps_inv, &
                   & this%det, this%Lmb, this%CConst, CSumq0)

        !allocate( this%CSumq0(shape(CSumq0)) )
        this%CSumq0 = CSumq0

        deallocate(Rpoints, Gpoints)
        deallocate(CConst, CSumq0)

    end subroutine set_EwaldParam


    subroutine set_DynEw(this, DynEwq)

        implicit none
        class(EwaldParam)                               :: this
        complex(dp), dimension(:,:,:), intent(in)       :: Dynewq

        alct_chk: if ( allocated(this%DynEwq) ) then

            write(*, 35)
            35 FORMAT('DynEw is already allocated, deallocating')
            deallocate(this%DynEwq)

        end if alct_chk  

        this%DynEwq = Dynewq

    end subroutine set_DynEw


    subroutine DecideEwaldParameters(eps, eps_inv, lat_vec, G, SupDim, prec_Ew, Ew_alpha, NKcut)
    
        implicit none
    
        real(dp), dimension(3, 3), intent(in)                   :: eps, eps_inv, &
                                                                 & lat_vec, G
        integer, dimension(3), intent(in)                       :: SupDim
        real(dp), intent(in)                                    :: prec_Ew
    
        real(dp), intent(out)                                   :: Ew_alpha
        integer, intent(out)                                    :: NKcut
    
        integer, dimension(3)           :: Ncut
        real(dp), dimension(3)          :: Rcut, dV, delV
        real(dp)                        :: p, D, Kcut, normG
    
        p = (-1.0_dp * dlog(prec_Ew))
    
        Ncut = (SupDim / 2) + 1
        Ncut(2:) = 0 !****!
        Rcut = matmul(lat_vec, dble(Ncut))
    
        dV = Rcut
        delV = matmul(dV, eps_inv)
        D = dot_product(delV, dV)
        D = dsqrt(D)
        Ew_alpha = dsqrt(p) / D
    
        normG = norm2(G(:,1))
        Kcut = 2.0_dp * p / (D * dsqrt( (eps(1,1) + eps(2,2) + eps(3,3))/3.0_dp ))
        NKcut = ceiling(KCut/normG)
    
    end subroutine DecideEwaldParameters

    subroutine ConstTerm(eps_inv, det_eps, Lmb, Nbasis, CConst)
    
        implicit none
    
        real(dp), dimension(3, 3), intent(in)               :: eps_inv
        real(dp), intent(in)                                :: det_eps, Lmb
        integer, intent(in)                                 :: Nbasis
        complex(dp), allocatable, intent(out)               :: CConst(:, :, :, :)
    
        real(dp)                                            :: cnst
        integer                                             :: mu, nu, alpha, beta
    
        allocate(CConst(3, 3, Nbasis, Nbasis))
        CConst = dcmplx(0.0_dp, 0.0_dp)
    
        do mu = 1, Nbasis
            do nu = 1, Nbasis
    
                if (mu == nu) then
                    do alpha = 1, 3
                        do beta = 1, 3
    
                            CConst(beta, alpha, nu, mu) = eps_inv(beta, alpha)
    
                        end do
                    end do
                end if
    
            end do
        end do
    
        cnst = 4.0_dp * (Lmb ** 3) / (3.0_dp * dsqrt(PI) * dsqrt(det_eps))
    
        CConst = CConst * cnst
    
    end subroutine ConstTerm

    subroutine SumRule(sys, Rpoints, Gpoints, eps_inv, det, Lmb, CConst, Csumq0)
    
        implicit none
    
        type(cell), intent(in)                                  :: sys
        real(dp), contiguous, dimension(:, :), intent(in)       :: Rpoints, Gpoints
        real(dp), dimension(3, 3), intent(in)                   :: eps_inv
        real(dp), intent(in)                                    :: det, Lmb
        complex(dp), contiguous, dimension(:,:,:,:), intent(in) :: CConst
        complex(dp), allocatable, intent(out)                   :: CSumq0(:, :, :)
    
        complex(dp), allocatable, dimension(:,:,:,:)            :: CRq0, CKq0, Cq0 &
                                                                &, Cq0Z
    
        real(dp), dimension(3)                                  :: q
        integer                                                 :: Nbasis, mu, nu, &
                                                                 & alpha, beta
    
        Nbasis = sys%natm
        allocate(Cq0(3, 3, Nbasis, Nbasis))
        Cq0 = dcmplx(0.0_dp, 0.0_dp)
    
        q = (/0.0_dp, 0.0_dp, 0.0_dp/)
        call RealSpaceEW(q, Rpoints, sys%tau, eps_inv, det, Lmb, CRq0)
        call KSpaceEW(q, Gpoints, sys%tau, sys%eps, sys%vol, Lmb, CKq0)
        Cq0 = CKq0 - CRq0 - CConst
        call CombineZstar(Cq0, sys%Zstar, Cq0Z)
    
        allocate(CSumq0(3, 3, Nbasis))
        CSumq0 = dcmplx(0.0_dp, 0.0_dp)
    
    mu_loop: do mu = 1, Nbasis
        a_loop: do alpha = 1, 3
            b_loop: do beta = 1, 3
    
                nu_loop: do nu = 1, Nbasis
                            CSumq0(beta, alpha, mu) = &
                          & CSumq0(beta, alpha, mu) + &
                          & Cq0Z(beta, alpha, nu, mu)
                         end do nu_loop
    
                    end do b_loop
                end do a_loop
             end do mu_loop
    
        deallocate(CRq0)
        deallocate(CKq0)
        deallocate(Cq0)
        deallocate(Cq0Z)
    
    end subroutine SumRule

end module EwaldMod

