module methods
    implicit none
contains
    ! Solve Schrodinger equation using given variables
    ! E : Energy eigenvalue
    ! phi : Wave function array
    ! x : Position samples used for phi
    ! a : Width of zero potential
    subroutine solve_schrodinger(E, phi, x, a)
        double precision,intent(in) :: E, a
        double precision,intent(inout) :: phi(0:), x(0:)
        double precision :: h, u
        integer i, n
        n = ubound(x, 1)
        h = a / dble(n)
        u = (phi(1) - phi(0)) / h

        do i = 0, n - 1
            call rk4(x(i), phi(i), x(i+1), phi(i+1), u, h, E)
        end do
    end subroutine
    ! Rungekutta
    ! x0 : position for last calculated phi
    ! f0 : last calculated phi
    ! x : position for calculating phi
    ! f : calculating phi
    ! v : defined as df/dx, last differentiated phi
    ! h : step of x
    ! E : energy eigenvalue
    subroutine rk4(x0, f0, x, f, v, h, E)
        double precision,intent(in) :: x0, f0, h, E
        double precision,intent(out) :: x, f
        double precision,intent(inout) :: v
        double precision kf(4), kv(4)

        kf(1) = h * df_dx(v)
        kv(1) = h * dv_dx(f0, E)

        kf(2) = h * df_dx(v + 0.5d0*kv(1))
        kv(2) = h * dv_dx(f0 + 0.5d0*kf(1), E)

        kf(3) = h * df_dx(v + 0.5d0*kv(2))
        kv(3) = h * dv_dx(f0 + 0.5d0*kf(2), E)

        kf(4) = h * df_dx(v + kv(3))
        kv(4) = h * dv_dx(f0 + kf(3), E)

        x = x0 + h
        v = v + (kv(1) + 2d0 * kv(2) + 2d0 * kv(3) + kv(4)) / 6d0
        f = f0 + (kf(1) + 2d0 * kf(2) + 2d0 * kf(3) + kf(4)) / 6d0
    end subroutine
    ! definition of df/dx
    double precision function df_dx(u)
        double precision,intent(in) :: u
        df_dx = u 
    end function
    ! definition of dv/dx
    double precision function dv_dx(phi, E)
        double precision,intent(in) :: phi, E
        dv_dx = - 2d0 * E * phi
    end function
    ! Count nodes for phi(x)=0
    integer function count_nodes(phi)
        double precision,intent(in) :: phi(0:)
        integer n, i
        n = ubound(phi, 1)
        count_nodes = 0

        do i = 0, n-1
            if (phi(i)*phi(i+1) < 0) then
                count_nodes = count_nodes + 1
            end if
        end do  
    end function
    ! Initialize phi and x with given conditions
    subroutine initialize(phi, x)
        double precision,intent(out) :: phi(0:), x(0:)
        phi = 0d0
        x = 0d0
        ! initial condition
        phi(0) = 0d0
        phi(1) = 1d-2
    end subroutine
end module

program main
    use methods
    implicit none
    ! n : wave function is devided by this number to get sample points
    integer,parameter :: n = 100
    ! phi : Wave function
    ! x : Sample positions of the wave function
    ! El : Minimum Energy Eigenvalue set initially
    ! Eh : Maximum Energy Eigenvalue set initially
    ! Em, Em0 : Average of El and Eh. Only used for calculation
    double precision :: phi(0:n), x(0:n), El = 0d0, Eh = 100d0, Em, Em0
    ! a : Width of zero potential
    ! PI : 3.141592.....
    ! epsilon : Extremly Small Number.
    double precision,parameter :: a = 10d0, PI=acos(-1d0), epsilon=1d-10
    integer i, Nm
    !integer  Nl, Nh

    ! Increase range to determine E more precisely
    ! Find an Energy Eigenvalue with condition of phi(a) = 0 by bisection method
    do i = 0, 10000
        !call initialize(phi, x)
        !call solve_schrodinger(El, phi, x, a)
        !Nl = count_nodes(phi)

        !call initialize(phi, x)
        !call solve_schrodinger(Eh, phi, x, a)
        !Nh = count_nodes(phi)

        call initialize(phi, x)
        call solve_schrodinger(Em, phi, x, a)
        Nm = count_nodes(phi)
        if (Nm == 0) then
            El = Em
        else
            Eh = Em
        end if
        Em0 = Em
        Em = 0.5d0 * (El + Eh)

        if (abs(Em0 - Em) < epsilon) then
            write (*, '(a, i0, a)') "Loop ended (n=", i, ")"
            exit
        end if
    end do
    write (*, '(a, f0.15)') "E     =", Em
    write (*, '(a, f0.15)') "Theory=", PI*PI/200d0
    write (*, '(a, f0.15)') "Error =", abs(Em - PI*PI/200d0)
end program