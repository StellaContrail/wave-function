module methods
    implicit none
    public calculate_energy, output
    private solve_schrodinger, fdm3, count_nodes, initialize, normalize
contains
    ! Solve Schrodinger equation using given variables
    ! E : Energy eigenvalue
    ! phi : Wave function array
    ! x : Position samples used for phi
    ! a : Width of zero potential
    subroutine solve_schrodinger(E, phi, x, a)
        double precision,intent(in) :: E, a
        double precision,intent(inout) :: phi(0:), x(0:)
        double precision :: h
        integer n
        n = ubound(x, 1)
        ! initialization
        h = (a - x(0)) / dble(n)

        call fdm5(x, phi, E, h, n)
    end subroutine
    ! Finite-Difference Methods
    subroutine fdm3(x, phi, E, h, n)
        double precision,intent(out) :: x(0:), phi(0:)
        double precision,intent(in) :: E, h
        integer,intent(in) :: n
        integer i
        ! initial condition
        phi(1) = 1d0

        do i = 0, n
            x(i) = x(0) + i * h
            if (1 <= i .and. i < n) then
                phi(i+1) = 2d0*phi(i) - phi(i-1) + h*h*f(x(i), phi(i), E)
            end if
        end do
    end subroutine
    subroutine fdm5(x, phi, E, h, n)
        double precision,intent(out) :: x(0:), phi(0:)
        double precision,intent(in) :: E, h
        integer,intent(in) :: n
        integer i

        do i = 0, n
            x(i) = x(0) + i * h
        end do
        ! initial condition
        phi(1) = 1d0

        !phi(2) = -12d0*h*h*f(x(0), phi(0), E)+16d0*phi(1)-30d0*phi(0)+16d0*(-phi(1))-(-phi(1)*2d0)
        !phi(3) = -12d0*h*h*f(x(1), phi(1), E)+16d0*phi(2)-30d0*phi(1)+16d0*phi(0)-(-phi(1))
        phi(2) = 2d0
        phi(3) = 3d0
        do i = 0, n
            if (2 <= i .and. i <= n - 2) then
                phi(i+2) = -12d0*h*h*f(x(i), phi(i), E)+16d0*phi(i+1)-30d0*phi(i)+16d0*phi(i-1)-phi(i-2)
            end if
        end do
    end subroutine
    double precision function f(x, phi, E)
        double precision,intent(in) :: phi, x, E
        f = (x*x - 2d0*E)*phi
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
        x = 0d0
        phi = 0d0
        ! initial condition
        x(0) = -5d0
        phi(0) = 0d0
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
            call initialize(phi, x)
            call solve_schrodinger(Em, phi, x, a)
            if (.false. .and. mod(i, 10) == 0) then
                write (*, *) "i=", i
                call output(x, phi)
                read (*, *)
            end if
            Nm = count_nodes(phi)

            if (Nm <= nodes) then
                El = Em
            else if (Nm > nodes) then
                Eh = Em
            end if
            Em0 = Em
            Em = 0.5d0 * (El + Eh)

            if (abs(Em0 - Em) < epsilon) then
                write (*, '(a, i0, a)') "Loop ended (n=", i, ")"
                exit
            end if

        end do

        calculate_energy = Em

        if (present(phi0)) then
            call normalize(phi, x(0), a)
            phi0 = phi
        end if
        if (present(x0)) then
            x0 = x
        end if
    end function
    ! Output phi and positions into a file
    subroutine output(x, phi)
        double precision,intent(in) :: x(0:), phi(0:)
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
    integer,parameter :: samples = 100
    integer n
    ! PI : 3.141592.....
    double precision,parameter :: PI=acos(-1d0)
    double precision :: E, theory_value, phi(0:samples), x(0:samples)
    
    write (*, *) "SOLVE Nth ENERGY EIGENVALUE"
    write (*, '(a)', advance="no") "Input N :"
    read (*, *) n

    E = calculate_energy(n, samples, 10d0, phi, x)
    write (*, '(a, f0.15)') "E     =", E
    theory_value = (n + 0.5d0)
    write (*, '(a, f0.15)') "Theory=", theory_value
    write (*, '(a, f0.15)') "Error =", abs(E - theory_value)
    call output(x, phi)
end program