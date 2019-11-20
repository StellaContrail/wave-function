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
    subroutine solve_eigen(A, Z, W, N)
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
        double precision,intent(out)     :: W(1:N)                  ! Eigenvalues in ascending order
        complex(kind(0d0)),intent(out)   :: Z(1:N, 1:N)             ! The ith column of Z holding the eigenvector associated with W(i)
        integer                          :: LDZ                     ! The first dimension of the array Z
        integer                          :: ISUPPZ(1:2*N)           ! The indices indicating the nonzero elements in Z
        double precision,allocatable     :: WORK(:)                 ! Workspace
        integer                          :: LWORK                   ! The dimension of the array work
        integer,allocatable              :: IWORK(:)                ! Workspace
        integer                          :: LIWORK                  ! The dimension of the array iwork
        integer                          :: INFO                    ! Success/Error indicator
        integer                          :: i                       ! Loop variable

        ! Substitute diagonal/non-diagonal elements of input matrix A
        do i = 1, N
            D(i) = A(i, i)
            if (i < N) then
                E(i) = A(i, i+1)
            end if
        end do
        
        LDZ = N

        allocate (WORK(1), IWORK(1))
        call zstegr(JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z, N, ISUPPZ, WORK, -1, IWORK, -1, INFO)
        LWORK  = int(WORK(1))
        LIWORK = int(WORK(1))
        
        deallocate(WORK, IWORK)
        allocate(WORK(LWORK), IWORK(LIWORK))
        call zstegr(JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z, N, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO)
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
    end subroutine solve_eigen

    ! Solve Ax = lambda*x for lambda and x
    subroutine solve_eigen2(A, Z, W, N)
        character,parameter              :: JOBZ  = 'V'             ! Eigenvalues and eigenvectors are computed
        character,parameter              :: RANGE = 'I'             ! The ILth through IUth eigenvectors will be found
        character,parameter              :: UPLO  = 'U'             ! The upper triangular part of A is stored
        integer,intent(in)               :: N                       ! Number of division in space
        double precision,intent(in)      :: A(1:N, 1:N)             ! Input Matrix
        integer                          :: LDA                     ! The first dimension of the array A
        double precision,parameter       :: VL    = 0d0             ! The lower bound of the interval to be searched for eigenvalues
        double precision,parameter       :: VU    = 0d0             ! The upper bound of the interval to be searched for eigenvalues
        integer,parameter                :: IL    = 1               ! The index of the smallest eigenvalues to be returned
        integer,parameter                :: IU    = 1               ! The index of the largest eigenvalues to be returned
        double precision                 :: ABSTOL                  ! Obsolete feature of LAPACK
        integer,parameter                :: M = IU - IL + 1         ! Total number of eigenvalues found
        double precision,intent(out)     :: W(1:N)                  ! Eigenvalues in ascending order
        complex(kind(0d0)),intent(out)   :: Z(1:N, 1:M)             ! The ith column of Z holding the eigenvector associated with W(i)
        integer                          :: LDZ                     ! The first dimension of the array Z
        double precision,allocatable     :: WORK(:)                 ! Workspace
        integer                          :: LWORK                   ! The dimension of the array work
        double precision,allocatable     :: RWORK(:)                ! Workspace
        integer,allocatable              :: IWORK(:)                ! Workspace
        integer,allocatable              :: JFAIL(:)                ! The indices of the eigenvectors that failed to converge
        integer                          :: INFO                    ! Success/Error indicator
        integer                          :: i                       ! Loop variable
        complex(kind(0d0))               :: A_(1:N, 1:N)            ! Workspace

        ! As the input matrix A_ would be partly overwritten, we don't want A to be changed.
        A_(:, :) = dcmplx(A(:, :), 0d0)
        ! LDZ is the first dimension of the array Z
        LDZ = N
        ! ABSTOL is the absolute error tolerance for the eigenvalues.
        ABSTOL = 2d0 * dlamch()
        ! Allocate the arrays
        allocate (RWORK(1:7*N), IWORK(1:5*N), JFAIL(1:N))

        ! Set parameters so the routine only calculates the optimal size of the WORK array.
        allocate (WORK(1))
        LWORK = -1
        ! Call ZHEEVX subroutine to calculate the optimal size of the WORK array.
        call zheevx(JOBZ, RANGE, UPLO, N, A_, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, JFAIL, INFO)
        LWORK  = int(WORK(1))
        ! Re-allocate the size of the array WORK
        deallocate(WORK)
        allocate(WORK(LWORK))
        call zheevx(JOBZ, RANGE, UPLO, N, A_, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, JFAIL, INFO)
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

    double precision function dlamch()
        double precision,parameter :: one  = 1d0
        double precision,parameter :: zero = 0d0
        double precision           :: rnd, eps, sfmin, small

        ! epsilon : The smallest number of the same kind as zero such that 1+eps > 1
        eps = epsilon(zero)

        ! tiny : The smallest positive number in the model of the type of zero
        sfmin = tiny(zero)

        ! huge : The largest number that is not an infinity in the model of the type of zero
        small = one / huge(zero)

        if ( small >= sfmin ) then
            sfmin = small*( one+eps )
        end if

        dlamch = sfmin
    end function dlamch
end module mathf