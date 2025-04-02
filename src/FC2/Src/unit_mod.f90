
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

module unit_cell

    use     kinds,          only : dp
    use     constants,      only : PI

    implicit none
    private

    type, public :: cell

        character(len=64)               :: prefix
        logical                         :: SecondFC, ThirdFC, FourthFC, Force_dd
        integer                         :: nn_2nd, nn_3rd, nn_4th

        integer                         :: ntyp, natm
        integer                         :: sup_cell(3)
        real(dp)                        :: latvec(3, 3)
        real(dp)                        :: vol

        integer, allocatable            :: typ_indx(:)
        real(dp), allocatable           :: basis_frac(:, :)
        character(len=3), allocatable   :: typ_lbl_unq(:)
        integer, allocatable            :: atomic_no_unq(:)
        real(dp), allocatable           :: mass_unq(:) 
        real(dp), allocatable           :: mass(:)
        real(dp)                        :: eps(3, 3)
        real(dp), allocatable           :: Zstar(:, :, :)

        real(dp), allocatable           :: tau(:, :)

        real(dp), dimension(3,3)        :: G

    contains
        procedure, public, pass     :: init_data
        procedure, private, nopass  :: cross

    end type cell


contains

    subroutine init_data(this, filename)
    
        implicit none

        class(cell)                                 :: this
        character(len=*), intent(in)                :: filename

        ! ================================== Local Variables ================================== !

        character(len=64)                           :: prefix
        logical                                     :: SecondFC, ThirdFC, FourthFC, Force_dd
        integer                                     :: nn_2nd, nn_3rd, nn_4th

        namelist    /FCInfo/    prefix, SecondFC, ThirdFC, FourthFC, Force_dd, nn_2nd, nn_3rd, nn_4th

        integer                                     :: Nbasis=2, Ntyp=1
        real(dp), dimension(3,3)                    :: Latvec=0.0_dp
        integer, dimension(3)                       :: Sup_Cell=5

        namelist    /CrystalInfo/   Nbasis, Ntyp, Latvec, Sup_Cell

        real(dp), allocatable, dimension(:, :)      :: atm_pos
        integer, allocatable, dimension(:)          :: typ_indx
        character(len=3), allocatable, dimension(:) :: typ_lbl
        integer, allocatable, dimension(:)          :: atomic_no
        real(dp), allocatable, dimension(:)         :: mass

        real(dp), dimension(3, 3)                   :: dielec = 1.0_dp
        real(dp), allocatable, dimension(:,:,:)     :: BornZ

        namelist    /AtmInfo/   atm_pos, typ_indx, typ_lbl, atomic_no, mass, dielec, BornZ

        integer                                     :: aa, i, i2, i3, atm_typ_indx
        real(dp)                                    :: vl
        integer                                     :: err
        character(len=512)                          :: err_msg

        ! ================================== Local Variables ================================== !

        open(unit=5, file=filename, status='OLD', iostat=err, iomsg=err_msg, &
           & action='READ', delim='APOSTROPHE')
            
            open_chk1: if ( err /= 0 ) then
                write(*, *) 'Input file OPEN failed: iostat = ', err
                write(*, *) 'Error message = ', err_msg
            end if open_chk1

            read(unit=5, nml=FCInfo)
            read(unit=5, nml=CrystalInfo)

        close(unit=5, status='KEEP', iostat=err)

        this%prefix = prefix
        this%SecondFC = SecondFC
        this%ThirdFC = ThirdFC        
        this%FourthFC = FourthFC
        this%Force_dd = Force_dd
        this%nn_2nd = nn_2nd
        this%nn_3rd = nn_3rd
        this%nn_4th = nn_4th

        this%natm = Nbasis
        this%ntyp = Ntyp
        this%sup_cell = Sup_Cell
        this%latvec = Latvec

        allocate(atm_pos(3, Nbasis))
        allocate(typ_indx(Nbasis))
        allocate(typ_lbl(Ntyp))
        allocate(atomic_no(Ntyp))
        allocate(mass(Ntyp))
        allocate(BornZ(3, 3, Nbasis))

        open(unit=5, file=filename, status='OLD', iostat=err, iomsg=err_msg, &
           & action='READ', delim='APOSTROPHE')
            
            open_chk2: if ( err /= 0 ) then
                write(*, *) 'Input file OPEN failed: iostat = ', err
                write(*, *) 'Error message = ', err_msg
            end if open_chk2

            read(unit=5, nml=AtmInfo)

        close(unit=5, status='KEEP', iostat=err)
    
        ! call move_alloc
        allocate(this%basis_frac(3, this%natm))
        allocate(this%typ_indx(this%natm))
        allocate(this%typ_lbl_unq(this%ntyp))
        allocate(this%atomic_no_unq(this%ntyp))
        allocate(this%mass_unq(this%ntyp))
        allocate(this%Zstar(3, 3, this%natm))
    
        this%basis_frac = atm_pos
        this%typ_indx = typ_indx
        this%typ_lbl_unq = typ_lbl
        this%atomic_no_unq = atomic_no
        this%mass_unq = mass
        this%eps = dielec
        this%Zstar = BornZ

        allocate(this%mass(this%natm))
        do i = 1, this%natm
            atm_typ_indx = this%typ_indx(i)
            this%mass(i) = this%mass_unq(atm_typ_indx)
        end do

        vl = dot_product(this%latvec(:, 1), &
                       & cross(this%latvec(:, 2),this%latvec(:, 3)))
        vl = dabs(vl)
        this%vol = vl

        allocate(this%tau(3, this%natm))
        do aa = 1, this%natm 
            this%tau(:, aa) = matmul(this%latvec, this%basis_frac(:, aa))
        end do

        do i = 1, 3
            i2 = mod((i+1), 3)
            if (i2 == 0) i2 = 3

            i3 = mod((i+2), 3)
            if (i3 == 0) i3 = 3

            this%G(:, i) = cross(this%latvec(:, i2), this%latvec(:, i3)) &
                         & * 2.0_dp * PI / this%vol
        end do

        deallocate(atm_pos)
        deallocate(typ_indx, typ_lbl, atomic_no, mass)
        deallocate(BornZ)
    
    end subroutine init_data


    function cross(a, b)

        implicit none
        real(dp), intent(in) :: a(:), b(:)
        real(dp), dimension(3) :: cross

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)

    end function cross

end module unit_cell

