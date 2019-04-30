! Used package : LAPACK95 (https://www.nag.co.uk/numeric/fl/nagdoc_fl26/pdf/f08/f08jbf.pdf)
!                         (http://www.nag-j.co.jp/lapack/dstevx.htm)
module extensions
    implicit none
contains
    ! solve eigenvalues, and eigenvectors of a real symmetric tridiagonal matrix A.
    subroutine solve_eigen(A, n, E, f, vl, vu, m)
        double precision,intent(in) :: A(:, :), vl, vu
        integer,intent(in) :: n
        integer,intent(out) :: m
        double precision,intent(out) :: E(:), f(:, :)
        double precision,allocatable :: A_diag(:), A_nondiag(:), work(:), iwork(:)
        ! il : smallest eigenvalue
        ! iu : largest eigenvalue
        ! m : the total number of eigenvalues found.
        ! info : INFO = 0 when the routine successfully executed.
        integer i, il, iu, info
        integer,allocatable :: jfail(:)
        double precision :: abstol = 0d0
        allocate(A_diag(1:n), A_nondiag(1:n-1), work(5*n), iwork(5*n), jfail(n))
        do i = 1, n-1
            A_diag(i) = A(i, i)
            A_nondiag(i) = A(i, i+1)
        end do
        A_diag(n) = A(n, n)
        call dstevx("V","V",n,A_diag,A_nondiag,vl,vu,il,iu,abstol,m,E,f,n,work,iwork,jfail,info)
    end subroutine
end module

program main
    use extensions
    implicit none
    integer,parameter :: n = 1000
    double precision,parameter :: PI = acos(-1d0), a = 10d0
    double precision :: H(1:n, 1:n) = 0d0, E(1:n), f(1:n,1:n), dh = a / dble(n)
    integer i, m
    
    do i = 1, n
        H(i, i) = -2d0
        if (i /= 1) then
            H(i, i-1) = 1d0
        end if
        if (i /= n) then
            H(i, i+1) = 1d0
        end if
    end do
    H = H * (-0.5d0 / (dh * dh))
    
    call solve_eigen(H, n, E, f, 0d0, 10d0, m)
    do i = 1, m
        write (*, '(a, i2, a, f6.4)') "Theoretical E(n=", i, ") = ", E(i)
        write (*, '(a, i2, a, f6.4)') "Solved      E(n=", i, ") = ", PI * PI * 0.5d0 * i * i / (a * a)
        write (*, *) "-----------------------------"
    end do
end program