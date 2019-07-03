module extension
    implicit none
    double precision,parameter :: COULOMB_K = 1d0
    double precision,parameter :: ANG_L = 0d0
contains
    subroutine solve_schrodinger(n, dr)
        integer,intent(in) :: n
        double precision,intent(in) :: dr
        double precision H(1:n, 1:n), r, values(1:n)
        double precision,allocatable :: dummy(:)
        integer i, lwork, info
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
        ! debug----------------------------------
        do i = 1, n
            write (*, '(100(f6.2, x))') H(i, :)
        end do
        ! ---------------------------------------
        write (*, *) -0.5d0 * (H(1,1)/12d0*dr*dr)
        write (*, *) -0.5d0 * H(1,1)/(12d0*dr*dr)
        H = -0.5d0 * (H/12d0*dr*dr)
        ! debug----------------------------------
        do i = 1, n
            write (*, '(100(f6.2, x))') H(i, :)
        end do
        stop
        ! ---------------------------------------
        
        r = 0d0
        do i = 1, n
            r = r + dr
            H(i, i) = H(i, i) + 0.5d0*ANG_L*(ANG_L+1d0)/(r*r) - COULOMB_K / r
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
            write (*, '(f10.5)') values(1)
            write (*, '(f10.5)') -0.5/(1d0)**2d0
        end if

        call normalize(H, n, dr)
        open(10, file="data.txt")
        do i = 1, n
            write(10, '(2(f10.5, x))') dr*i, H(i, 1)
        end do
    end subroutine

    subroutine normalize(u, n, dr)
        integer,intent(in) :: n
        double precision,intent(inout) :: u(1:n, 1:n)
        double precision,intent(in) :: dr
        integer i
        double precision sum, r
        sum = 0d0
        r = 0d0
        do i = 1, n
            if (i == 1 .or. i == n) then
                sum = sum + 0.5d0 * u(i, 1)**2d0 * dr
            else
                sum = sum + u(i, 1)**2d0 * dr
            end if
        end do
        do i = 1, n
            u(i, 1) = u(i, 1) / sqrt(sum)
        end do
    end subroutine
end module

program main
    use extension
    implicit none
    double precision :: R = 10d0, dr
    integer :: n = 500
    dr = R / n
    call solve_schrodinger(5, dr)
    call solve_schrodinger(n, dr)
end program