module extensions
    implicit none
    double precision,parameter :: ALPHA = 1d0
    double precision,parameter :: COULOMB_K = 1d0
contains
    subroutine solve_schroedinger(u, v, n, dr, l, E, startIndex)
        integer,intent(in) :: n, l, startIndex
        double precision,intent(inout) :: u(startIndex:n), v(startIndex:n)
        double precision,intent(in) :: dr, E
        double precision r
        integer i
        do i = startIndex, n-1
            r = dr * i
            call rungekutta(u(i), v(i), u(i+1), v(i+1), r, dr, l, E)
        end do
    end subroutine
    subroutine initialize(u, v, n, startIndex)
        integer,intent(in) :: n, startIndex
        double precision,intent(out) :: u(startIndex:n), v(startIndex:n)
        if (startIndex == 0) then
            u(0) = 0d0
            v(0) = 1d0
        else
            u(1) = 0d0
            v(1) = 1d0
        end if
    end subroutine

    subroutine rungekutta(u, v, u_new, v_new, r, dr, l, E)
        double precision,intent(in) :: u, v, dr, r, E
        double precision,intent(out) :: u_new, v_new
        double precision ku(1:4), kv(1:4)
        integer,intent(in) :: l
        kv(1) = dr*dv_dr(u,l,r,E)
        ku(1) = dr*du_dr(v)

        kv(2) = dr*dv_dr(u+0.5d0*ku(1),l,r+0.5d0*dr,E)
        ku(2) = dr*du_dr(v+0.5d0*kv(1))

        kv(3) = dr*dv_dr(u+0.5d0*ku(2),l,r+0.5d0*dr,E)
        ku(3) = dr*du_dr(v+0.5d0*kv(2))

        kv(4) = dr*dv_dr(u+ku(3),l,r+dr,E)
        ku(4) = dr*du_dr(v+kv(3))

        v_new = v + (kv(1)+2d0*kv(2)+2d0*kv(3)+kv(4))/6d0
        u_new = u + (ku(1)+2d0*ku(2)+2d0*ku(3)+ku(4))/6d0
    end subroutine
    double precision function du_dr(v)
        double precision,intent(in) :: v
        du_dr = v
    end function
    double precision function dv_dr(u, l, r, E)
        double precision,intent(in) :: u, r, E
        integer,intent(in) :: l
        double precision,parameter :: R_zero = 0.1d0
        if (r == 0d0 .and. l == 0d0) then
            dv_dr = -2d0*COULOMB_K*ALPHA*R_zero
        else
            dv_dr = (l*(l+1)/(r*r)-2d0*COULOMB_K*ALPHA/r-2d0*E*ALPHA)*u
        end if
    end function
    subroutine output_to_file(u, n, dr, startIndex)
        integer,intent(in) :: n, startIndex
        double precision,intent(in) :: u(startIndex:n), dr
        integer i
        open(10, file="data.txt")
        do i = startIndex, n
            write (10, *) dr*i, u(i)
        end do
        close(10)
    end subroutine
    integer function count_nodes(u, n, startIndex)
        integer,intent(in) :: n, startIndex
        double precision,intent(in) :: u(startIndex:n)
        integer i
        count_nodes = 0
        do i = startIndex, n-1
            if (u(i)*u(i+1) < 0d0) then
                count_nodes = count_nodes + 1
            end if
        end do
    end function
    subroutine normalize(u, n, dr, startIndex)
        integer,intent(in) :: n, startIndex
        double precision,intent(inout) :: u(startIndex:n)
        double precision,intent(in) :: dr
        double precision sum
        integer i
        sum = 0d0
        do i = startIndex, n
            if (i == startIndex .or. i == n)  then
                sum = sum + 0.5d0*u(i)*u(i)*dr
            else
                sum = sum + u(i)*u(i)*dr
            end if
        end do
        u(:) = u(:) / sqrt(sum)
    end subroutine
end module

program main
    use extensions
    implicit none
    double precision :: R = 20d0
    integer,parameter :: n = 5000
    double precision,allocatable :: u(:), v(:)
    double precision EL, EM, EH, dr
    integer i, l, k, startIndex
    l = 0
    k = 0
    if (l == 0) then
        startIndex = 0
    else
        startIndex = 1
    end if
    allocate(u(startIndex:n), v(startIndex:n))
    call initialize(u, v, n, startIndex)
    dr = R / n

    EL = -1d0
    EH = 0d0
    do i = 1, 1000
        EM = 0.5d0*(EL+EH)
        call solve_schroedinger(u, v, n, dr, l, EM, startIndex)
        if (count_nodes(u,n,l) <= k) then
            EL = EM
        else if (count_nodes(u,n,l) > k) then
            EH = EM
        end if
    end do
    write (*, *) "numerical  :", EM
    write (*, *) "analytical :", -0.5d0/((l+k+1)*(l+k+1))
    write (*, '(x, a, f6.3, a)') "error      :", (1d0-abs((EM/(0.5d0/((l+k+1)*(l+k+1))))))*100d0, " %"
    call normalize(u, n, dr, l)
    call output_to_file(u, n, dr, l)
end program