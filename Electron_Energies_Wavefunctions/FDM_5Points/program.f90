module extension
    implicit none
    double precision,parameter :: COULOMB_K = 1d0
    double precision,parameter :: ANG_L = 1d0
    integer,parameter :: MAX_N = 3
contains
    subroutine solve_schrodinger(n, dr)
        integer,intent(in) :: n
        double precision,intent(in) :: dr
        double precision H(1:n, 1:n), r, values(1:n)
        double precision,allocatable :: dummy(:)
        integer i, j, lwork, info
        H = 0d0
        do i = 1, n
            H(i, i) = -30d0
            if (1 < i) then
                H(i, i-1) = 16d0
            end if
            if (2 < i) then
                H(i, i-2) = -1d0
            end if
            if (i < n) then
                H(i, i+1) = 16d0
            end if
            if (i < n-1) then
                H(i, i+2) = -1d0
            end if
        end do
        H = - 0.5d0 * H / (12d0*dr*dr)
        ! debug----------------------------------
        !do i = 1, n
        !    write (*, '(100(f6.2, x))') H(i, :)
        !end do
        !stop
        ! ---------------------------------------
        
        r = 0d0
        do i = 1, n
            r = r + dr
            H(i, i) = H(i, i) + 0.5d0*ANG_L*(ANG_L+1d0)/(r*r) + V(r)
        end do
        
        allocate(dummy(1))
        call dsyev("V", "U", n, H, n, values, dummy, -1, info)
        lwork = int(dummy(1))
        deallocate(dummy)
        allocate(dummy(lwork))
        call dsyev("V", "U", n, H, n, values, dummy, lwork, info)
        if (info < 0) then
            write (*, *) "Error when executing dsyev() manipulation"
        else
            write (*, *)
            write (*, '(a, f2.0)') "Orbital Angular Momentum : L = ", ANG_L
            do i = 1, MAX_N
                write (*, '(a, i0, a)') "n = ", i, " ----------------------"
                write (*, '(a, f13.8)') "E_NUMERICAL  =", values(i)
                write (*, '(a, f13.8)') "E_ANALYTICAL =", -0.5/(i+ANG_L)**2d0
                write (*, *)
            end do
        end if

        call normalize(H, n, dr)
        open(10, file="data.txt")
        do j = 1, MAX_N
            do i = 1, n
                write(10, '(2(f10.5, x))') dr*i, H(i, j)
            end do
            write (10, *)
        end do
    end subroutine

    double precision function V(r)
        double precision,intent(in) :: r
        V = - COULOMB_K / r
    end function

    subroutine normalize(u, n, dr)
        integer,intent(in) :: n
        double precision,intent(inout) :: u(1:n, 1:n)
        double precision,intent(in) :: dr
        integer i, j
        double precision sum, r
        r = 0d0
        do j = 1, MAX_N
            sum = 0d0
            do i = 1, n
                if (i == 1 .or. i == n) then
                    sum = sum + 0.5d0 * u(i, j)**2d0 * dr
                else
                    sum = sum + u(i, j)**2d0 * dr
                end if
            end do
            u(:, j) = u(:, j) / sqrt(sum)
        end do
    end subroutine
end module

program main
    use extension
    implicit none
    double precision :: R = 50d0, dr
    integer :: n = 500
    dr = R / n
    call solve_schrodinger(n, dr)
end program