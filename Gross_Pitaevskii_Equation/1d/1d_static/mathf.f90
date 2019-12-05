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
        double precision,intent(inout)  :: f(0:N)
        double precision sum
        call Integrate(abs(f(:))**2d0, N, dh, sum)
        f(:) = f(:) / sqrt(sum)
    end subroutine normalize

    ! Solve Ax = lambda*x for lambda and x
    subroutine solve_eigen(A, Z, lambda, N)
        character,parameter              :: JOBZ  = 'V'             ! Eigenvalues and eigenvectors are computed
        character,parameter              :: UPLO  = 'U'             ! The upper triangular part of A is stored
        integer,intent(in)               :: N                       ! Number of division in space
        double precision,intent(in)      :: A(1:N, 1:N)             ! Input Matrix
        integer                          :: LDA                     ! The first dimension of the array A
        double precision                 :: W(1:N)                  ! Eigenvalues in ascending order
        double precision,allocatable     :: WORK(:)                 ! Workspace
        integer                          :: LWORK                   ! The dimension of the array work
        integer                          :: INFO                    ! Success/Error indicator
        double precision                 :: A_(1:N, 1:N)            ! Workspace
        double precision,intent(out)     :: Z(1:N)                  ! Eigenvector
        double precision,intent(out)     :: lambda                  ! Eigenvalue

        ! As the input matrix A_ would be partly overwritten, we don't want A to be changed.
        A_(:, :) = A(:, :)
        ! LDZ is the first dimension of the array Z
        LDA = N

        ! Set parameters so the routine only calculates the optimal size of the WORK array.
        allocate (WORK(1))
        LWORK = -1
        ! Call CHEEV subroutine to calculate the optimal size of the WORK array.
        call dsyev(JOBZ, UPLO, N, A_, LDA, W, WORK, LWORK, INFO)
        LWORK  = int(WORK(1))
        ! Re-allocate the size of the array WORK
        deallocate(WORK)
        allocate(WORK(LWORK))
        ! Actual calculation of eigenvalue equation
        call dsyev(JOBZ, UPLO, N, A_, LDA, W, WORK, LWORK, INFO)
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
    end subroutine solve_eigen

    ! Shift the phase of input complex vector f by Phase
    ! f        : COMPLEX array having dimension of N
    ! N        : Integer dimension of F
    ! iu       : Imaginary unit. sqrt(-1)
    ! Phase    : The phase of array f would be shifted by Phase
    ! f_result : The phase-shifted input vector would be substituted into this variable
    subroutine apply_phase_shift(f, N, iu, Phase, f_result)
        integer,intent(in)             :: N
        double precision,intent(in)    :: Phase
        double precision,intent(in)    :: f(0:N)
        complex(kind(0d0)),intent(in)  :: iu
        complex(kind(0d0)),intent(out) :: f_result(0:N)
        
        f_result(:) = exp(-iu*Phase)*f(:)
    end subroutine


    subroutine calc_asymmetry_degree(Phi, N, Degree)
        integer,intent(in)           :: N
        double precision,intent(in)  :: Phi(0:N)
        double precision,intent(out) :: Degree
        integer                      :: i
        Degree = 0d0
        do i = 0, floor(N/2d0)
            Degree = Degree + abs(Phi(ceiling(N/2d0) + i) - Phi(floor(N/2d0) - i))
        end do
    end subroutine
end module mathf