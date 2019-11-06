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
        double precision,intent(in)  :: f(1:N), dh
        double precision,intent(out) :: sum
        integer i
        sum = 0d0
        do i = 1, N
            if (i == 1 .or. i == N) then
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
        complex(kind(0d0)),intent(inout)  :: f(1:N)
        double precision sum
        call Integrate(abs(f(:))**2d0, N, dh, sum)
        f(:) = f(:) / sqrt(sum)
    end subroutine normalize

    ! Calculate C := exp(Ax)
    ! A : REAL array having dimension of NxN
    ! x : Complex array having dimension of N
    ! C : Complex array having dimension of N
    subroutine exp_mat(A, f, N, dt, epsilon, iu, ans)
        integer,intent(in)             :: N
        double precision,intent(in)    :: A(1:N, 1:N), dt, epsilon
        complex(kind(0d0)),intent(in)  :: f(1:N), iu
        complex(kind(0d0)),intent(out) :: ans(1:N)
        integer i
        complex(kind(0d0))             :: temp(1:N), Atemp(1:N)

        ! First term of Taylor expansion
        temp(:)   = f(:)
        ans(:) = temp(:)

        ! Other terms of Taylor expansion
        do i = 1, 4
            call multiply_symmetry(A, temp, N, Atemp)
            temp(:)   = -iu*Atemp*dt/(epsilon*i)
            ans(:) = ans(:) + temp(:)
        end do
    end subroutine exp_mat

    ! Calculate C := AB
    ! A : REAL array having dimension of NxN
    ! B : COMPLEX array having dimension of N
    ! C : COMPLEX array having dimension of N
    subroutine multiply_symmetry(A, B, N, C)
        integer,intent(in)             :: N
        double precision,intent(in)    :: A(1:N, 1:N)
        complex(kind(0d0)),intent(in)  :: B(1:N)
        complex(kind(0d0)),intent(out) :: C(1:N)
        integer                        :: i

        do i = 1, N
            if (i == 1) then
                C(i) = A(i,i)*B(i)+A(i,i+1)*B(i+1)+A(i,i+2)*B(i+2)+A(i,i+3)*B(i+3)+A(i,i+4)*B(i+4)
            else if (i == 2) then
                C(i) = A(i,i-1)*B(i-1)+A(i,i)*B(i)+A(i,i+1)*B(i+1)+A(i,i+2)*B(i+2)+A(i,i+3)*B(i+3)+A(i,i+4)*B(i+4)
            else if (i == 3) then
                C(i) = A(i,i-2)*B(i-2)+A(i,i-1)*B(i-1)+A(i,i)*B(i)+A(i,i+1)*B(i+1)+A(i,i+2)*B(i+2)+A(i,i+3)*B(i+3)+A(i,i+4)*B(i+4)
            else if (i == 4) then
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
    subroutine expected_value_symm(f, A, N, ans)
        integer,intent(in)             :: N
        double precision,intent(in)    :: A(1:N, 1:N)
        complex(kind(0d0)),intent(in)  :: f(1:N)
        double precision,intent(out)   :: ans
        complex(kind(0d0))             :: Af(1:N), temp
        call multiply_symmetry(A, f, N, Af)
        temp = dot_product(f, Af)
        if (aimag(temp) > 1d-5) then
            print *, "ERROR : Possible calculation error at expected_value_symm"
        end if
        ans = dble(temp)
    end subroutine
end module mathf