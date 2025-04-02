
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


module mklWrap

    use kinds,      only : dp
    implicit none

    EXTERNAL            :: dgetrf, dgetri 
    private
    public  :: dinverse, ddet

contains

    subroutine dinverse(a)
    
        implicit none
    
        real(dp), contiguous, dimension(:, :), intent(inout)    :: a
        integer, dimension(:), allocatable                      :: ipiv   ! pivot indices
    
        real(dp), dimension(:), allocatable                     :: work
        integer                                                 :: lwork
        integer                                                 :: n, lda, info
    
        n = size(a, 1)
        allocate(ipiv(n))
    
        !** dgetrf:=> Computes the LU factorization of a general m-by-n matrix: A = P*L*U **!
        lda = n
        call dgetrf(n, n, a, lda, ipiv, info)
    
        if (info /= 0) then
            stop 'eps matrix is numerically singular!'
        end if
        !** dgetrf:=> Computes the LU factorization of a general m-by-n matrix: A = P*L*U **!
    
        !** dgetri:=> Computes the inverse of an LU-factored general matrix **!
        lwork = -1
        allocate(work(max(1, lwork)))
    
        call dgetri(n, a, lda, ipiv, work, lwork, info)
    
        lwork = int(work(1))
        deallocate(work)
        allocate(work(max(1, lwork)))
    
        call dgetri(n, a, lda, ipiv, work, lwork, info)
    
        if (info /= 0) then
            stop 'Matrix inversion failed!'
        end if
        !** dgetri:=> Computes the inverse of an LU-factored general matrix **!
    
        deallocate(ipiv)
        deallocate(work)
    
    end subroutine dinverse
    
    
    subroutine ddet(a, det)
    
        implicit none
    
        real(dp), contiguous, dimension(:, :), intent(in)       :: a
        real(dp), intent(out)                                   :: det
    
        integer, dimension(:), allocatable                      :: ipiv   ! pivot indices
        integer                                                 :: n, lda, info
    
        real(dp), dimension(:, :), allocatable                  :: a_loc
        integer                                                 :: ii
    
        n = size(a, 1)
    
        allocate(a_loc(n, size(a, 2)))
        a_loc = a
    
        allocate(ipiv(n))
    
        !** dgetrf:=> Computes the LU factorization of a general m-by-n matrix: A = P*L*U **!
        lda = n
        call dgetrf(n, n, a_loc, lda, ipiv, info)
    
        if (info /= 0) then
            Write(*, *) 'WARNING: eps matrix is numerically singular!'
        end if
        !** dgetrf:=> Computes the LU factorization of a general m-by-n matrix: A = P*L*U **!
    
        !** Computes the determinant of an PLU-factored general matrix **!
        det = 1.0_dp
        do ii = 1, n
            
            det = det * a_loc(ii, ii)
            if (ipiv(ii) /= ii) det = det * (-1.0_dp)
    
        end do
        !** Computes the determinant of an PLU-factored general matrix **!
    
        deallocate(ipiv)
        deallocate(a_loc)
    
    end subroutine ddet

end module mklWrap

