module methods
    implicit none
    public calculate_energy, output
    private solve_schrodinger, rk4, df_dx, dv_dx, count_nodes, initialize, normalize
contains
    ! Solve Schrodinger equation using given variables
    ! E : Energy eigenvalue
    ! phi : Wave function array
    ! x : Position samples used for phi
    ! a : Width of zero potential
    subroutine solve_schrodinger(E, phi, x, a)
        double precision,intent(in) :: E, a
        double precision,intent(inout) :: phi(0:), x(0:)
        double precision :: h, v
        integer i, n
        n = ubound(x, 1)
        ! initialization
        h = (a - x(0)) / dble(n)
        v = (phi(1) - phi(0)) / h

        do i = 0, n - 1
            call rk4(x(i), phi(i), x(i+1), phi(i+1), v, h, E)
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
        kv(1) = h * dv_dx(f0, E, x0)

        kf(2) = h * df_dx(v + 0.5d0*kv(1))
        kv(2) = h * dv_dx(f0 + 0.5d0*kf(1), E, x0 + 0.5d0 * h)

        kf(3) = h * df_dx(v + 0.5d0*kv(2))
        kv(3) = h * dv_dx(f0 + 0.5d0*kf(2), E, x0 + 0.5d0 * h)

        kf(4) = h * df_dx(v + kv(3))
        kv(4) = h * dv_dx(f0 + kf(3), E, x0 + h)

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
    double precision function dv_dx(phi, E, x)
        double precision,intent(in) :: phi, E, x
        dv_dx = - (2d0 * E - x*x) * phi
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
        x(0) = -5d0
        phi(1) = 1d-2
    end subroutine
    ! Calculate Energy Eigenvalue with arbitrary nodes count
    double precision function calculate_energy(nodes, n, Eh0, phi0, x0)
        integer,intent(in) :: nodes
        ! n : wave function is devided by this number to get sample points
        integer,intent(in) :: n
        double precision,intent(in) :: Eh0
        ! phi0 : [Optional parameter] Normalized Wave function will be returned
        ! x0 : [Optional parameter] Sample positions of Wave function will be returned
        double precision,intent(out),optional :: phi0(0:), x0(0:)
        ! phi : Wave function
        ! x : Sample positions of the wave function
        ! El : Minimum Energy Eigenvalue set initially
        ! Eh : Maximum Energy Eigenvalue set initially
        ! Em, Em0 : Average of El and Eh. Only used for calculation
        double precision phi(0:n), x(0:n), El, Eh, Em, Em0
        ! a : Width of zero potential
        ! PI : 3.141592.....
        ! epsilon : Extremly Small Number.
        double precision,parameter :: a = 5d0, epsilon=1d-10
        integer i, Nm
        !integer  Nl, Nh

        ! Initialization
        El = 0d0
        Eh = Eh0 ! Eh must be bigger than E

        ! Increase range to determine E more precisely
        ! Find an Energy Eigenvalue with condition of phi(a) = 0 by bisection method
        Em = 0.5d0 * (El + Eh)
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

            !if (i == 0) then
            !    call output(phi, x)
            !    write (*, *) "Nm=", Nm
            !    write (*, *) "Eh=", Eh, " Em=", Em, " El=", El
            !    read (*, *)
            !end if

            ! Nm<=nodesのとき、解はElとEmの間にあるはずなので、捜索範囲をEl ~ Emに変更する
            if (Nm <= nodes) then
                El = Em
            else if (Nm > nodes) then
                Eh = Em
            end if
            Em0 = Em
            Em = 0.5d0 * (El + Eh)

            if (abs(Em0 - Em) < epsilon) then
                write (*, *) Nm
                write (*, '(a, i0, a)') "Loop ended (n=", i, ")"
                exit
            end if
        end do
        ! nが小さいとき、サンプル数が少なすぎるあまりに最初の根の数が的確に取られない可能性がある。
        ! 例：n=100 -> 根の数:12
        ! そのため、指定された根の数よりも大きい根の数を持つ場合のnを使うのが適切である

        calculate_energy = Em

        if (present(phi0)) then
            call normalize(phi, phi(0), a)
            phi0 = phi
        end if
        if (present(x0)) then
            x0 = x
        end if
    end function
    ! Output phi and positions into a file
    subroutine output(phi, x)
        double precision,intent(in) :: phi(0:), x(0:)
        integer n, i
        integer,parameter :: fo = 10
        n = ubound(phi, 1)

        open(fo, file="data.txt")
        do i = 0, n
            write (fo, *) x(i), phi(i)
        end do
        close(fo)
    end subroutine
    ! Normalize wave function
    ! phi : Wave function
    ! start : Start value of integration
    ! end : End value of integration
    subroutine normalize(phi, start, end)
        double precision,intent(in) :: start, end
        double precision,intent(inout) :: phi(0:)
        double precision sum, h
        integer i, n
        sum = 0d0
        n = ubound(phi, 1)
        h = (end - start) / n

        do i = 0, n
            if (i == 0 .or. i == n) then
                sum = sum + 0.5d0 * abs(phi(i) * phi(i)) * h
            else
                sum = sum + abs(phi(i) * phi(i)) * h
            end if
        end do

        do i = 0, n
            phi(i) = phi(i) / sqrt(sum)
        end do
    end subroutine
end module

program main
    use methods
    implicit none
    ! Number of samples of wave function
    integer,parameter :: samples = 1000
    integer n
    ! PI : 3.141592.....
    double precision,parameter :: PI=acos(-1d0)
    double precision :: E, theory_value, phi(0:samples), x(0:samples)
    
    write (*, *) "SOLVE Nth ENERGY EIGENVALUE"
    write (*, '(a)', advance="no") "Input N :"
    read (*, *) n

    E = calculate_energy(n, samples, 1d4, phi, x)
    write (*, '(a, f0.15)') "E     =", E
    theory_value = (n + 0.5d0)
    write (*, '(a, f0.15)') "Theory=", theory_value
    write (*, '(a, f0.15)') "Error =", abs(E - theory_value)
    call output(phi, x)
end program