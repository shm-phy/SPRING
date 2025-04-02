
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

    EXTERNAL            :: dgetrf, dgetri, zheevd, dgelsd, dgelsy
    private

    public  :: dinverse, ddet, SolveEigen, LU_decompose, RankDefLSqr, RankDefLSqr2

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
    
        integer                                                 :: n
        character(len=1)                                        :: jobz, uplo
        integer                                                 :: lda, lwork, liwork, lrwork, info
        complex(dp), dimension(:), allocatable                  :: work
        integer, dimension(:), allocatable                      :: iwork
        real(dp), dimension(:), allocatable                     :: rwork
    
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
    
    
    subroutine LU_decompose(A, Aout)
    
        implicit none
    
        real(dp), dimension(:,:), intent(inout)             :: A
        real(dp), dimension(:,:), allocatable, intent(out)  :: Aout
    
        ! =================================== Local Variables =================================== !
    
        real(dp), dimension(:,:), allocatable   :: Mat_tmp
        integer                                 :: m, n, lda, info
        integer, dimension(:), allocatable      :: ipiv
        logical                                 :: warn
    
        ! =================================== Local Variables =================================== !
    
        !warn = .true.
        warn = .false.

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
            STOP

        end if error_chk
    
        size_chk: if ( m < n ) then
            write(*, 150) m, n
            write(*, 160)
            STOP
    
        else size_chk
            allocate( Aout(n, n) )
            allocate( Mat_tmp(n, n) )
            Mat_tmp(:, :) = A(1:n, 1:n)
            Aout = transpose( Mat_tmp )

            deallocate( Mat_tmp )
    
        end if size_chk
    
        deallocate( ipiv )

        65 FORMAT(36X, " ** WARNING ** ")
        85 FORMAT(37X, "U(", I2, ",", I2, ") = 0.0.", /, &
                & 16X, "The factorization has been completed, but U is exactly singular.", /, &
                & 6X, "Division by 0 will occur if you use the factor U for solving a system of linear equations")
    
        125 FORMAT("LU decomposition routine (dgetrf) fails! Error code: ", I3)
        150 FORMAT("For present case num of row must be greater than number of column: [", I6, "x", I6, "]")
        160 FORMAT("Use zero padding!")
    
    end subroutine LU_decompose


    subroutine RankDefLSqr( A, x, Xout, SingValCut )

        implicit none

        real(dp), dimension(:, :), intent(inout)            :: A
        real(dp), dimension(:), intent(inout)               :: x
        real(dp), dimension(:), allocatable, intent(out)    :: Xout
        real(dp), intent(in)                                :: SingValCut

        ! ==================================== Local Variables ==================================== !

        integer                                 :: m, &     !No of Rows of A 
                                                 & n, &     !No of Cols of A
                                                 & nrhs, &  !The number of right-hand sides; the number of columns in B
                                                 & lda, &   !The leading dimension of A; at least max(1, m)
                                                 & ldb, &   !The leading dimension of b; must be at least max(1, m, n)
                                                 & lwork

        real(dp)                                :: rcond    !rcond is used to determine the effective rank of A. Singular values 
                                                            ! s(i) ≤ rcond *s(1) are treated as zero

        real(dp), dimension(:), allocatable     :: work     !work is a workspace array, its dimension max(1, lwork)
                                                            !on exit, work(1) contains the minimum value of lwork
                                                            !required for optimum performance

        integer, dimension(:), allocatable      :: iwork    !on exit, iwork(1) returns the minimum size of the workspace
                                                            !array iwork required for optimum performance.

        real(dp), dimension(:), allocatable     :: s        !Array, size at least max(1, min(m, n))
        integer                                 :: rank     !The effective rank of A
        integer                                 :: info, workSz, iworkSz

        real(dp)                                :: Res

        ! ==================================== Local Variables ==================================== !

        write(*, *)
        write(*, 20)

        ! ** !
        rcond = SingValCut
        ! ** !

        m = size(A, 1)
        n = size(A, 2)
        nrhs = 1
        lda = max(1, m)
        ldb = max(1, m, n)

        allocate( s(max(1, min(m, n))) )

        ! **** !
        lwork = -1      ! workspace query is assumed; the routine only
                        ! calculates the optimal size of the array work and the minimum sizes of the
                        ! arrays iwork , and returns these values as the first entries of the
                        ! work and iwork arrays
        ! **** !

        allocate( work(max(1, lwork)) )
        allocate( iwork(max(1, lwork)) )

        call dgelsd (m, n, nrhs, A, lda, x, ldb, s, rcond, rank, work, lwork, iwork, info)

        workSz = int( work(1) )
        iworkSz = iwork(1)

        deallocate( work )
        deallocate( iwork )

        ! **** !
        lwork = workSz
        ! **** !

        allocate( work(max(1, lwork)) )
        allocate( iwork(iworkSz) )

        call dgelsd (m, n, nrhs, A, lda, x, ldb, s, rcond, rank, work, lwork, iwork, info)

        ErrorChk: if ( info < 0 ) then
            write(*, 50) abs(info)
            STOP

        else if ( info > 0 ) then ErrorChk
            write(*, 100) info
            STOP

        else ErrorChk
            write(*, 150)
            write(*, 200) m, n, rank

            Res = dot_product( x(n+1:m), x(n+1:m) )
            write(*, 250) Res

            write(*, *)

        end if ErrorChk

        allocate( Xout(n) )
        Xout(:) = x(1:n)

        deallocate( work, iwork, s )

        20 FORMAT("Starting Linear least squares routine (RankDefLSqr)")
        50 FORMAT("ERROR(in RankDefLSqr) : The ", I2, "-th parameter had an illegal value")
        100 FORMAT("ERROR(in RankDefLSqr) : The algorithm for computing the SVD failed to converge ", I6)
        150 FORMAT("Linear least squares routine (RankDefLSqr) is succesful")
        200 FORMAT("Rank of the (", I6, " x ", I6, ") matrix is : ", I6)
        250 FORMAT("Residue ||b-Ax||^2: ", F12.4)

    end subroutine RankDefLSqr
    

    subroutine RankDefLSqr2( A, x, Xout, SingValCut )

        implicit none

        real(dp), dimension(:, :), intent(inout)            :: A
        real(dp), dimension(:), intent(inout)               :: x
        real(dp), dimension(:), allocatable, intent(out)    :: Xout
        real(dp), intent(in)                                :: SingValCut

        !*! call dgelsy (m, n, nrhs, A, lda, x, ldb, jpvt, rcond, rank, work, lwork, info)
        ! ==================================== Local Variables ==================================== !

        integer                                 :: m, &     !No of Rows of A 
                                                 & n, &     !No of Cols of A
                                                 & nrhs, &  !The number of right-hand sides; the number of columns in B
                                                 & lda, &   !The leading dimension of A; at least max(1, m)
                                                 & ldb, &   !The leading dimension of b; must be at least max(1, m, n)
                                                 & lwork

        real(dp)                                :: rcond    !rcond is used to determine the effective rank of A. Singular values 
                                                            ! s(i) ≤ rcond *s(1) are treated as zero

        real(dp), dimension(:), allocatable     :: work     !work is a workspace array, its dimension max(1, lwork)
                                                            !on exit, work(1) contains the minimum value of lwork
                                                            !required for optimum performance

        integer, dimension(:), allocatable      :: jpvt     !Array, size at least max(1, n).
                                                            !On entry, if jpvt(i)≠ 0, the i-th column of A is permuted to the front of AP,
                                                            !otherwise the i-th column of A is a free column.

        integer                                 :: rank     !The effective rank of A
        integer                                 :: info, workSz

        !-!real(dp)                                :: Res

        ! ==================================== Local Variables ==================================== !

        write(*, *)
        write(*, 20)

        ! ** !
        rcond = SingValCut
        ! ** !

        m = size(A, 1)
        n = size(A, 2)
        nrhs = 1
        lda = max(1, m)
        ldb = max(1, m, n)

        allocate( jpvt(max(1, n)) )
        jpvt = 0

        ! **** !
        lwork = -1      ! workspace query is assumed; the routine only
                        ! calculates the optimal size of the array work and the minimum sizes of the
                        ! arrays work and returns these values as the first entries of the work arrays
        ! **** !

        allocate( work(max(1, lwork)) )

        call dgelsy (m, n, nrhs, A, lda, x, ldb, jpvt, rcond, rank, work, lwork, info)

        workSz = int( work(1) )

        deallocate( work )

        ! **** !
        lwork = workSz
        ! **** !

        allocate( work(max(1, lwork)) )
        jpvt = 0

        call dgelsy (m, n, nrhs, A, lda, x, ldb, jpvt, rcond, rank, work, lwork, info)

        ErrorChk: if ( info < 0 ) then
            write(*, 50) abs(info)
            STOP

        else ErrorChk
            write(*, 150)
            write(*, 200) m, n, rank

            !-!Res = dot_product( x(n+1:m), x(n+1:m) )
            !-!write(*, 250) Res
            write(*, *)

        end if ErrorChk

        allocate( Xout(n) )
        Xout(:) = x(1:n)

        deallocate( work, jpvt )

        20 FORMAT("Starting Linear least squares routine (RankDefLSqr2)")
        50 FORMAT("ERROR(in RankDefLSqr2) : the ", I2, "-th parameter had an illegal value")
        150 FORMAT("Linear least squares routine (RankDefLSqr2) is succesful")
        200 FORMAT("Rank of the (", I6, " x ", I6, ") matrix is : ", I6)
        !-!250 FORMAT("Residue ||b-Ax||^2: ", F12.4)

    end subroutine RankDefLSqr2

end module mklWrap


