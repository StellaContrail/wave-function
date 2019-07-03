module extensions
    implicit none
    double precision,parameter :: ALPHA = 1d0
    double precision,parameter :: COULOMB_K = 1d0
contains
    subroutine solve_schroedinger(u, v, n, dr, l, E)
        integer,intent(in) :: n, l
        double precision,intent(inout) :: u(0:n), v(0:n)
        double precision,intent(in) :: dr, E
        double precision r
        integer i
        do i = 0, n-1
            r = dr * i
            call rungekutta(u(i), v(i), u(i+1), v(i+1), r, dr, l, E)
        end do
    end subroutine
    subroutine initialize(u, v, n)
        integer,intent(in) :: n
        double precision,intent(out) :: u(0:n), v(0:n)
        u(0) = 0d0
        v(0) = 1d0
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
    subroutine output_to_file(u, n, dr)
        integer,intent(in) :: n
        double precision,intent(in) :: u(0:n), dr
        integer i
        open(10, file="data.txt")
        do i = 0, n
            write (10, *) dr*i, u(i)
        end do
        close(10)
    end subroutine
    integer function count_nodes(u, n)
        integer,intent(in) :: n
        double precision,intent(in) :: u(0:n)
        integer i
        count_nodes = 0
        do i = 1, n-1
            if (u(i)*u(i+1) < 0d0) then
                count_nodes = count_nodes + 1
            end if
        end do
    end function
    subroutine normalize(u, n, dr)
        integer,intent(in) :: n
        double precision,intent(inout) :: u(0:n)
        double precision,intent(in) :: dr
        double precision sum
        integer i
        sum = 0d0
        do i = 1, n
            if (i == 1 .or. i == n)  then
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
    double precision u(0:n), v(0:n), dr
    double precision EL, EM, EH
    integer i, l, k
    call initialize(u, v, n)
    l = 0
    k = 0
    dr = R / n

    EL = -1d0
    EH = 0d0
    do i = 1, 1000
        EM = 0.5d0*(EL+EH)
        call solve_schroedinger(u, v, n, dr, l, EM)
        if (count_nodes(u,n) <= k) then
            EL = EM
        else if (count_nodes(u,n) > k) then
            EH = EM
        end if
    end do
    write (*, *) "numerical  :", EM
    write (*, *) "analytical :", -0.5d0/((k+1)*(k+1))
    write (*, '(x, a, f6.3, a)') "error      :", (1d0-abs((EM/(0.5d0/((k+1)*(k+1))))))*100d0, " %"
    call normalize(u, n, dr)
    call output_to_file(u, n, dr)
end program