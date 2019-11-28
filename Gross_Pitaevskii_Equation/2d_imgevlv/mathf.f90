! Mathematical Procedures
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
        double precision,intent(in)  :: f(0:N, 0:N), dh
        double precision,intent(out) :: sum
        double precision             :: sum_temp(0:N)
        integer                      :: i
        sum = 0d0
        sum_temp(:) = 0d0

        do i = 0, N
            if (i == 0 .or. i == N) then
                sum_temp(:) = sum_temp(:) + 0.5d0*f(i,:)*dh
            else
                sum_temp(:) = sum_temp(:) + f(i,:)*dh
            end if
        end do

        do i = 0, N
            if (i == 0 .or. i == N) then
                sum = sum + 0.5d0*sum_temp(i)*dh
            else
                sum = sum + sum_temp(i)*dh
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
        double precision,intent(inout)    :: f(0:N, 0:N)
        double precision                  :: sum
        call integrate(abs(f(:, :))**2d0, N, dh, sum)
        f(:, :) = f(:, :) / sqrt(sum)

        write (*, *) "Integration result : ", sum
    end subroutine normalize

    ! Calculate C := exp(Ax)
    ! A : REAL array having dimension of NxN
    ! x : Complex array having dimension of N
    ! C : Complex array having dimension of N
    subroutine exp_mat(A, f, N, dt, epsilon, iu, ans)
        integer,intent(in)             :: N
        double precision,intent(in)    :: A(0:N, 0:N), dt, epsilon
        complex(kind(0d0)),intent(in)  :: f(0:N), iu
        complex(kind(0d0)),intent(out) :: ans(0:N)
        integer                        :: i
        complex(kind(0d0))             :: temp(0:N), Atemp(0:N)

        ! First term of Taylor expansion
        temp(:)   = f(:)
        ans(:) = temp(:)

        ! Other terms of Taylor expansion
        do i = 1, 4
            call multiply_symmetry(A, temp, N, Atemp)
            temp(:)   = -Atemp*dt/(epsilon*i)
            ans(:) = ans(:) + temp(:)
        end do
    end subroutine exp_mat

    ! Calculate C := AB
    ! A : REAL array having dimension of NxN
    ! B : COMPLEX array having dimension of N
    ! C : COMPLEX array having dimension of N
    subroutine multiply_symmetry(A, B, N, C)
        integer,intent(in)             :: N
        double precision,intent(in)    :: A(0:N, 0:N)
        complex(kind(0d0)),intent(in)  :: B(0:N)
        complex(kind(0d0)),intent(out) :: C(0:N)
        integer                        :: i

        do i = 0, N
            if (i == 0) then
                C(i) = A(i,i)*B(i)+A(i,i+1)*B(i+1)+A(i,i+2)*B(i+2)+A(i,i+3)*B(i+3)+A(i,i+4)*B(i+4)
            else if (i == 1) then
                C(i) = A(i,i-1)*B(i-1)+A(i,i)*B(i)+A(i,i+1)*B(i+1)+A(i,i+2)*B(i+2)+A(i,i+3)*B(i+3)+A(i,i+4)*B(i+4)
            else if (i == 2) then
                C(i) = A(i,i-2)*B(i-2)+A(i,i-1)*B(i-1)+A(i,i)*B(i)+A(i,i+1)*B(i+1)+A(i,i+2)*B(i+2)+A(i,i+3)*B(i+3)+A(i,i+4)*B(i+4)
            else if (i == 3) then
                C(i) = A(i,i-3)*B(i-3)+A(i,i-2)*B(i-2)+A(i,i-1)*B(i-1)+A(i,i)*B(i)+A(i,i+1)*B(i+1)+A(i,i+2)*B(i+2)+A(i,i+3)*B(i+3)&
                +A(i,i+4)*B(i+4)
            else if (i == N-3) then
                C(i) = A(i,i-4)*B(i-4)+A(i,i-3)*B(i-3)+A(i,i-2)*B(i-2)+A(i,i-1)*B(i-1)+A(i,i)*B(i)+A(i,i+1)*B(i+1)+A(i,i+2)*B(i+2)&
                +A(i,i+3)*B(i+3)
            else if (i == N-2) then
                C(i) = A(i,i-4)*B(i-4)+A(i,i-3)*B(i-3)+A(i,i-2)*B(i-2)+A(i,i-1)*B(i-1)+A(i,i)*B(i)+A(i,i+1)*B(i+1)+A(i,i+2)*B(i+2)
            else if (i == N-1) then
                C(i) = A(i,i-4)*B(i-4)+A(i,i-3)*B(i-3)+A(i,i-2)*B(i-2)+A(i,i-1)*B(i-1)+A(i,i)*B(i)+A(i,i+1)*B(i+1)
            else if (i == N) then
                C(i) = A(i,i-4)*B(i-4)+A(i,i-3)*B(i-3)+A(i,i-2)*B(i-2)+A(i,i-1)*B(i-1)+A(i,i)*B(i)
            else
                C(i) = A(i,i-4)*B(i-4)+A(i,i-3)*B(i-3)+A(i,i-2)*B(i-2)+A(i,i-1)*B(i-1)+A(i,i)*B(i)+A(i,i+1)*B(i+1)+A(i,i+2)*B(i+2)&
                +A(i,i+3)*B(i+3)+A(i,i+4)*B(i+4)
            end if
        end do
    end subroutine multiply_symmetry

    ! Calculate expected value of a symmetric matrix A
    ! i.e. calculate ans := <f|A|f>
    ! A   : REAL array having dimension of NxN
    ! f   : COMPLEX array having dimension of N
    ! ans : REAL value
    subroutine expected_value_symm(f, A, N, ans, dh)
        integer,intent(in)             :: N
        double precision,intent(in)    :: A(0:N, 0:N), dh
        complex(kind(0d0)),intent(in)  :: f(0:N)
        double precision,intent(out)   :: ans
        complex(kind(0d0))             :: Af(0:N), temp
        call multiply_symmetry(A, f, N, Af)
        ! Note that f is coordinate representation (that means 1 = ∫<f|x><x|f>*dh is correct)
        temp = dot_product(f, Af)*dh
        if (aimag(temp) > 1d-5) then
            print *, "ERROR : Possible calculation error at expected_value_symm"
        end if
        ans = dble(temp)
    end subroutine

    ! Calculate probability current at given time
    subroutine calc_current(f, N, dh, hbar, mass, j)
        integer,intent(in)            :: N
        double precision,intent(in)   :: dh, hbar, mass
        complex(kind(0d0)),intent(in) :: f(0:N)
        double precision,intent(out)  :: j(0:N)
        integer                       :: i
        do i = 0, N
            if (i == 0) then
                j(0) = hbar * aimag( conjg(f(0))*(f(1)-f(0))/dh ) / mass
            else if (i == N) then
                j(N) = hbar * aimag( conjg(f(N))*(f(N)-f(N-1))/dh ) / mass
            else
                j(i) = hbar * aimag( conjg(f(i))*(f(i+1)-f(i-1))/(2d0*dh) ) / mass
            end if
        end do
    end subroutine

    ! Solve Ax = lambda*x for lambda and x
    subroutine solve_eigen(A, Z, lambda, N)
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
        complex(kind(0d0)),intent(in)  :: f(0:N), iu
        complex(kind(0d0)),intent(out) :: f_result(0:N)
        
        f_result(:) = exp(-iu*Phase)*f(:)
    end subroutine
end module mathf