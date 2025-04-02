
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

    EXTERNAL            :: dgetrf, dgetri, zgetrf, zgetri, zheevd, dgetrs

    private
    public  :: dinverse, zinverse, ddet, SolveEigen, LU_decompose, LinearSysSolve

contains

    subroutine dinverse(a)
    
        implicit none
    
        real(dp), contiguous, dimension(:, :), intent(inout)    :: a

        ! =================================== Local Variables =================================== !
        integer, dimension(:), allocatable                      :: ipiv   ! pivot indices
    
        real(dp), dimension(:), allocatable                     :: work
        integer                                                 :: lwork
        integer                                                 :: n, lda, info
        ! =================================== Local Variables =================================== !
    
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
    

    subroutine zinverse(a, not_invertible)
    
        implicit none

        real(dp), parameter                                     :: EPS = 1.0E-8_dp
    
        complex(dp), dimension(:, :), intent(inout)             :: a
        logical, intent(out)                                    :: not_invertible

        ! =================================== Local Variables =================================== !
        integer, dimension(:), allocatable                      :: ipiv   ! pivot indices
    
        !-!character(len=1)                                        :: sign_str
        complex(dp), dimension(:), allocatable                  :: work
        complex(dp)                                             :: det

        real(dp)                                                :: abs_det

        integer                                                 :: lwork
        integer                                                 :: n, lda, info, ii
        ! =================================== Local Variables =================================== !
    
        not_invertible = .false.
        n = size(a, 1)
        allocate(ipiv(n))
    
        !** zgetrf:=> Computes the LU factorization of a general m-by-n matrix: A = P*L*U **!
        lda = n
        call zgetrf(n, n, a, lda, ipiv, info)
    
        if (info /= 0) then
            STOP 'eps matrix is numerically singular! (zinverse) '
        end if
        !** zgetrf:=> Computes the LU factorization of a general m-by-n matrix: A = P*L*U **!

        !** Computes the determinant of an PLU-factored general matrix **!
        det = dcmplx( 1.0_dp, 0.0_dp )
        do ii = 1, n
            
            det = det * a(ii, ii)
            if (ipiv(ii) /= ii) det = det * (-1.0_dp)
    
        end do
        !** Computes the determinant of an PLU-factored general matrix **!

        !**            check whether 'a' matrix is singular            **!
        abs_det = abs( det )
        if( abs_det < EPS ) not_invertible = .true.

        !-! sign_str = ' '
        !-! if (dimag(det) > 0.0_dp) sign_str = '+'
        !-! write(*, "( 'Determinant = ', G20.12, A1, G20.12, 'j' )") dble(det), sign_str, dimag(det)
        !**            check whether 'a' matrix is singular            **!
    
        checkForInvertible: if ( .not. not_invertible ) then

            !** zgetri:=> Computes the inverse of an LU-factored general matrix **!
            lwork = -1
            allocate(work(max(1, lwork)))
    
            call zgetri(n, a, lda, ipiv, work, lwork, info)

            lwork = int( dble(work(1)) )
            deallocate(work)
            allocate(work(max(1, lwork)))
    
            call zgetri(n, a, lda, ipiv, work, lwork, info)
    
            if (info /= 0) then
                STOP 'Matrix inversion (zinverse) failed!'
            end if
            !** zgetri:=> Computes the inverse of an LU-factored general matrix **!
    
            deallocate(work)

        end if checkForInvertible
    
        deallocate(ipiv)

    end subroutine zinverse
    
    
    subroutine ddet(a, det)
    
        implicit none
    
        real(dp), contiguous, dimension(:, :), intent(in)       :: a
        real(dp), intent(out)                                   :: det
    
        ! =================================== Local Variables =================================== !
        integer, dimension(:), allocatable                      :: ipiv   ! pivot indices
        integer                                                 :: n, lda, info
    
        real(dp), dimension(:, :), allocatable                  :: a_loc
        integer                                                 :: ii
        ! =================================== Local Variables =================================== !
    
        n = size(a, 1)
    
        allocate(a_loc(n, size(a, 2)))
        a_loc = a
    
        allocate(ipiv(n))
    
        !** dgetrf:=> Computes the LU factorization of a general m-by-n matrix: A = P*L*U **!
        lda = n
        call dgetrf(n, n, a_loc, lda, ipiv, info)
    
        if (info /= 0) then
            write(*, *) 'WARNING: eps matrix is numerically singular!'
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
    
    
    subroutine SolveEigen(A, W)
    
        implicit none
        
        complex(dp), dimension(:, :), intent(inout)             :: A
        real(dp), dimension(:), allocatable, intent(out)        :: W
    
        ! =================================== Local Variables =================================== !
        integer                                                 :: n
        character(len=1)                                        :: jobz, uplo
        integer                                                 :: lda, lwork, liwork, lrwork, info
        complex(dp), dimension(:), allocatable                  :: work
        integer, dimension(:), allocatable                      :: iwork
        real(dp), dimension(:), allocatable                     :: rwork
        ! =================================== Local Variables =================================== !
    
        n = size(A, 1)
        allocate(W(n))
    
        lwork = -1
        allocate(work(max(1, lwork)))
        liwork = -1
        allocate(iwork(max(1, liwork)))
        lrwork = -1
        allocate(rwork(max(1, liwork)))
    
        jobz = 'V' !both eigenvalues and eigenvectors are computed
        uplo = 'U' !a stores the upper triangular part of A
        !uplo = 'L' !a stores the lower triangular part of A
    
        lda = n
    
        call zheevd(jobz, uplo, n, A, lda, W, work, lwork, rwork, lrwork, iwork, liwork, info)
    
        lwork = int(work(1))
        deallocate(work)
        allocate(work(max(1, lwork)))
    
        liwork = int(iwork(1))
        deallocate(iwork)
        allocate(iwork(max(1, liwork)))
    
        lrwork = int(rwork(1))
        deallocate(rwork)
        allocate(rwork(max(1, lrwork)))
    
        call zheevd(jobz, uplo, n, A, lda, W, work, lwork, rwork, lrwork, iwork, liwork, info)
    
        if (info /= 0) then
            write(*, 98) info
            98 FORMAT('Eigen solver routine fails! Error code: ', I3)
        end if
    
    end subroutine SolveEigen

    
    subroutine LU_decompose(A, ipiv)
    
        implicit none
    
        real(dp), dimension(:,:), intent(inout)             :: A
        integer, dimension(:), allocatable, intent(out)     :: ipiv
    
        ! =================================== Local Variables =================================== !
    
        integer                                 :: m, n, lda, info
        logical                                 :: warn
    
        ! =================================== Local Variables =================================== !
    
        warn = .true.
        !warn = .false.

        m = size(A, 1)
        n = size(A, 2)
    
        lda = m
    
        allocate( ipiv( max(1, min(m, n)) ) )
    
        call dgetrf(m, n, A, lda, ipiv, info)
    
        error_chk: if ( info > 0) then

            warning: if ( warn ) then
                write(*, *)
                write(*, 65)
                write(*, 85) info, info
                write(*, 65)
                write(*, *)
            end if warning

        else if ( info < 0) then error_chk
            write(*, 125) info
            ERROR STOP

        end if error_chk
    
        65 FORMAT(36X, " ** WARNING ** ")
        85 FORMAT(37X, "U(", I2, ",", I2, ") = 0.0.", /, &
                & 16X, "The factorization has been completed, but U is exactly singular.", /, &
                & 6X, "Division by 0 will occur if you use the factor U for solving a system of linear equations")
    
        125 FORMAT("LU decomposition routine (dgetrf) fails! Error code: ", I3)
    
    end subroutine LU_decompose


    subroutine LinearSysSolve( A, b, ipiv )

        implicit none

        real(dp), dimension(:, :), intent(in)           :: A
        real(dp), dimension(:), intent(inout)           :: b
        integer, dimension(:), intent(in)               :: ipiv

        ! ================================= Local variables ================================= !

        character(len=1)                                :: trans
        integer                                         :: n, nrhs, lda, ldb, info

        ! ================================= Local variables ================================= !

        trans = 'N'
        n = size( A, 1 )
        nrhs = 1
        lda = max( 1, n)
        ldb = max(1, n)

        call dgetrs( trans, n, nrhs, A, lda, ipiv, b, ldb, info )
        
        error_chk: if ( info /= 0) then 
            write(*, 125) info
            ERROR STOP

        end if error_chk

        125 FORMAT(" Linear equation solver routine (dgetrs) fails! Error code: ", I3)

    end subroutine LinearSysSolve


end module mklWrap

