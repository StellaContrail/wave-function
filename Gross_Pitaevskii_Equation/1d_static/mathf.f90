module mathf
    implicit none
contains
    ! Integration of function f using Trapezoidal rule
    ! f   : Integrand array
    ! N   : Dimension of space excluding the first element
    ! dh  : Step distance of space
    ! sum : The result of the integration
    subroutine integrate(f, N, dh, sum)
        integer,intent(in)           :: N
        double precision,intent(in)  :: f(0:N), dh
        double precision,intent(out) :: sum
        integer i
        sum = 0d0
        do i = 0, N
            if (i == 0 .or. i == N) then
                sum = sum + 0.5*f(i)*dh
            else
                sum = sum + f(i)*dh
            end if
        end do
    end subroutine

    ! Normalize the given function f
    ! f   : Function to be normalized
    ! N   : Dimension of f
    ! dh  : step distance of space
    subroutine normalize(f, N, dh)
        integer,intent(in)                :: N
        double precision,intent(in)       :: dh
        complex(kind(0d0)),intent(inout)  :: f(0:N)
        double precision sum
        call Integrate(abs(f(:))**2d0, N, dh, sum)
        f(:) = f(:) / sqrt(sum)
    end subroutine normalize

    ! Solve Ax = lambda*x for lambda and x
    subroutine solve_eigen(A, Z, lambda, N)
        integer,intent(in)               :: N                       ! Number of division in space
        double precision,intent(in)      :: A(1:N, 1:N)             ! Input Matrix
        character,parameter              :: JOBZ  = 'V'             ! Eigenvalues and eigenvectors are computed
        character,parameter              :: RANGE = 'I'             ! The ILth through IUth eigenvectors will be found
        double precision                 :: D(1:N)                  ! The n diagonal elements of the tridiagonal matrix
        double precision                 :: E(1:N)                  ! The subdiagonal elements of the tridiagonal matrix
        double precision,parameter       :: VL    = 0d0             ! The lower bound of the interval to be searched for eigenvalues
        double precision,parameter       :: VU    = 0d0             ! The upper bound of the interval to be searched for eigenvalues
        integer,parameter                :: IL    = 1               ! The index of the smallest eigenvalues to be returned
        integer,parameter                :: IU    = 1               ! The index of the largest eigenvalues to be returned
        double precision                 :: ABSTOL                  ! Obsolete feature of LAPACK
        integer                          :: M                       ! Total number of eigenvalues found
        double precision                 :: W(1:N)                  ! Eigenvalues in ascending order
        complex(kind(0d0))               :: Z_(1:N, 1:N)            ! The ith column of Z holding the eigenvector associated with W(i)
        complex(kind(0d0)),intent(out)   :: Z(1:N)                  ! The array holding the eigenvector associated with lambda
        integer                          :: LDZ                     ! The first dimension of the array Z
        integer                          :: ISUPPZ(1:2*N)           ! The indices indicating the nonzero elements in Z
        double precision,allocatable     :: WORK(:)                 ! Workspace
        integer                          :: LWORK                   ! The dimension of the array work
        integer,allocatable              :: IWORK(:)                ! Workspace
        integer                          :: LIWORK                  ! The dimension of the array iwork
        integer                          :: INFO                    ! Success/Error indicator
        integer                          :: i                       ! Loop variable
        double precision,intent(out)     :: lambda                  ! Eigenvalue

        ! Substitute diagonal/non-diagonal elements of input matrix A
        do i = 1, N
            D(i) = A(i, i)
            if (i < N) then
                E(i) = A(i, i+1)
            end if
        end do
        
        LDZ = N

        allocate (WORK(1), IWORK(1))
        call zstegr(JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z_, N, ISUPPZ, WORK, -1, IWORK, -1, INFO)
        LWORK  = int(WORK(1))
        LIWORK = int(WORK(1))
        
        deallocate(WORK, IWORK)
        allocate(WORK(LWORK), IWORK(LIWORK))
        call zstegr(JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z_, N, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO)
        deallocate(WORK, IWORK)

        if (INFO < 0) then
            write (*, *) "ERROR : Argument ", -INFO, " has illegal value"
        else if (INFO == 1) then
            write (*, *) "ERROR : The dqds algorithm failed to converge"
            stop
        else if (INFO == 2) then
            write (*, *) "ERROR : Inverse iteration failed to converge"
            stop
        end if

        lambda = W(1)
        Z(:) = Z_(:, 1)
    end subroutine solve_eigen

    ! Solve Ax = lambda*x for lambda and x
    subroutine solve_eigen2(A, Z, lambda, N)
        character,parameter              :: JOBZ  = 'V'             ! Eigenvalues and eigenvectors are computed
        character,parameter              :: UPLO  = 'U'             ! The upper triangular part of A is stored
        integer,intent(in)               :: N                       ! Number of division in space
        double precision,intent(in)      :: A(1:N, 1:N)             ! Input Matrix
        integer                          :: LDA                     ! The first dimension of the array A
        double precision                 :: W(1:N)                  ! Eigenvalues in ascending order
        complex(kind(0d0)),allocatable   :: WORK(:)                 ! Workspace
        integer                          :: LWORK                   ! The dimension of the array work
        double precision                 :: RWORK(1:3*N-2)          ! Workspace
        integer                          :: INFO                    ! Success/Error indicator
        complex(kind(0d0))               :: A_(1:N, 1:N)            ! Workspace
        complex(kind(0d0))               :: Z(1:N)                  ! Eigenvector
        double precision,intent(out)     :: lambda                  ! Eigenvalue

        ! As the input matrix A_ would be partly overwritten, we don't want A to be changed.
        A_(:, :) = dcmplx(A(:, :), 0d0)
        ! LDZ is the first dimension of the array Z
        LDA = N

        ! Set parameters so the routine only calculates the optimal size of the WORK array.
        allocate (WORK(1))
        LWORK = -1
        ! Call CHEEV subroutine to calculate the optimal size of the WORK array.
        call zheev(JOBZ, UPLO, N, A_, LDA, W, WORK, LWORK, RWORK, INFO)
        LWORK  = int(WORK(1))
        ! Re-allocate the size of the array WORK
        deallocate(WORK)
        allocate(WORK(LWORK))
        ! Actual calculation of eigenvalue equation
        call zheev(JOBZ, UPLO, N, A_, LDA, W, WORK, LWORK, RWORK, INFO)
        deallocate(WORK)

        if (INFO < 0) then
            write (*, *) "ERROR : Argument ", -INFO, " has illegal value"
        else if (INFO == 1) then
            write (*, *) "ERROR : The dqds algorithm failed to converge"
            stop
        else if (INFO == 2) then
            write (*, *) "ERROR : Inverse iteration failed to converge"
            stop
        end if

        lambda = W(1)
        Z(:)   = A_(:, 1)
    end subroutine solve_eigen2

    ! Shift the phase of input complex vector f by Phase
    ! f        : COMPLEX array having dimension of N
    ! N        : Integer dimension of F
    ! iu       : Imaginary unit. sqrt(-1)
    ! Phase    : The phase of array f would be shifted by Phase
    ! f_result : The phase-shifted input vector would be substituted into this variable
    subroutine apply_phase_shift(f, N, iu, Phase, f_result)
        integer,intent(in)             :: N
        double precision,intent(in)    :: Phase
        complex(kind(0d0)),intent(in)  :: f(0:N), iu
        complex(kind(0d0)),intent(out) :: f_result(0:N)
        
        f_result(:) = exp(-iu*Phase)*f(:)
    end subroutine
end module mathf