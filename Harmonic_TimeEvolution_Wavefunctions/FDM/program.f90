module extensions
    implicit none
    double precision,parameter :: ALPHA = 1d0 ! (HBAR*C)^2/(2*MC^2)
    double precision,parameter :: BETA = 1d0 ! (MASS*OMEGA)^2
    double precision,parameter :: HBAR = 1d0
    double precision,parameter :: K = -10d0 ! wave number
    double precision,parameter :: a = 1d0
contains
    ! i*(a+ib)=-b+ia : (REAL=-b), (IMAG=a)
    function ix(A, n)
        integer,intent(in) :: n
        double complex,intent(in) :: A(1:n)
        double complex ix(1:n)
        integer i
        do i = 1, n
            ix(i) = complex(-aimag(A(i)), dble(A(i)))
        end do
    end function
    function ix_scaler(z)
        double complex,intent(in) :: z
        double complex ix_scaler
        ix_scaler = complex(-aimag(z), dble(z))
    end function

    subroutine solve_schroedinger(phi, n, m, dh, dt, xl)
        integer,intent(in) :: n, m
        double complex,intent(inout) :: phi(1:n, 1:m)
        double precision,intent(in) :: dh, dt, xl
        double precision H(1:n, 1:n), x
        double precision t1, t2, speed
        integer i
        H = 0d0
        do i = 1, n
            x = xl + dh * i
            call construct_hamiltonian(H, x, n, dh)
            !H(i, i) = 0.5d0
        end do
        
        phi = complex(0d0, 0d0)
        call initialize(phi(:, 1), n, dh, xl)
        call normalize(phi(:, 1), n, dh)

        call cpu_time(t1)
        do i = 1, m-1
            phi(:, i+1) = calc_future_wavefunction(phi(:,i), H, n, dt)
            if (mod(i, 100) == 0) then
                call cpu_time(t2)
                speed = 100d0/(t2-t1)
                write (*, *) (100d0*i)/m, " % ", speed, " items/sec", " ETA : ", (m-1-i)/speed, " sec"
                call cpu_time(t1)
            endif
        end do
    end subroutine

    subroutine normalize(phi, n, dh)
        integer,intent(in) :: n
        double complex,intent(inout) :: phi(1:n)
        double precision,intent(in) :: dh

        phi = phi / sqrt(calc_probability(phi, n, dh))
    end subroutine

    double precision function calc_probability(phi, n, dh)
        integer,intent(in) :: n
        double complex,intent(in) :: phi(1:n)
        double precision,intent(in) :: dh
        double precision sum
        integer i

        sum = 0d0
        do i = 1, n
            if (i == 1 .or. i == n) then
                sum = sum + 0.5d0*abs(phi(i))**2d0*dh
            else
                sum = sum + abs(phi(i))**2d0*dh
            end if
        end do
        calc_probability = sum
    end function

    subroutine initialize(phi, n, dh, xl)
        integer,intent(in) :: n
        double complex,intent(out) :: phi(1:n)
        double precision,intent(in) :: dh, xl
        integer i
        double precision x

        do i = 1, n
            x = xl + dh * i
            phi(i) = complex(cos(K*x)*exp(-0.5d0*x**2d0/a), sin(K*x)*exp(-0.5d0*x**2d0/a))
        end do  
    end subroutine  

    function calc_future_wavefunction(phi, H, n, dt)
        integer,intent(in) :: n
        double complex,intent(in) :: phi(1:n)
        double precision,intent(in) :: dt, H(1:n, 1:n)
        double complex :: calc_future_wavefunction(1:n)
        double complex :: phi_next(1:n), phi_temp(1:n)
        integer i

        phi_next = complex(0d0, 0d0)

        ! n = 0
        phi_temp = phi
        phi_next = phi_next + phi_temp

        i = 1
        do
            phi_temp = -ix(matmul(H, phi_temp), n)*dt / dble(i)
            phi_next = phi_next + phi_temp
            if (maxval(abs(dble(phi_temp))) < 1d-10 .and. maxval(abs(aimag(phi_temp))) < 1d-10) then
                !write (*, *) i, "Convergence"
                exit
            end if
            i = i + 1
        end do

        calc_future_wavefunction = phi_next
    end function

    subroutine construct_hamiltonian(H, x, n, dh)
        integer,intent(in) :: n
        double precision,intent(out) :: H(1:n, 1:n)
        double precision,intent(in)  :: x, dh
        integer i

        H = 0d0
        do i = 1, n
            H(i, i) = -30d0
            if (i > 1) then
                H(i, i-1) = 16d0
            end if
            if (i > 2) then
                H(i, i-2) = -1d0
            end if
            if (i < n-1) then
                H(i, i+1) = 16d0
            end if
            if (i < n-2) then
                H(i, i+2) = -1d0
            end if
        end do

        do i = 1, n
            H(i, i) = H(i, i) / (12d0*dh*dh)
            H(i, i) = - ALPHA * H(i, i)
            H(i, i) = H(i, i) + V(x)
        end do
    end subroutine
    double precision function V(x)
        double precision,optional,intent(in) :: x
        V = 0.5d0*x*x*BETA
    end function
    
    subroutine plot(phi, n, m, dh, xl)
        integer,intent(in) :: n, m
        double precision,intent(in) :: dh, xl
        double complex,intent(in) :: phi(1:n, 1:m)
        integer i, j
        double precision x
        open(10, file="data.txt")

        do j = 1, m
            do i = 1, n
                x = xl + dh * i
                write (10, *) x, dble(phi(i, j)), aimag(phi(i, j)), abs(phi(i,j))**2d0
            end do
            write (10, *)
        end do

        close(10)
    end subroutine
    
    function solve_x(phi, n, m, dh, xl)
        integer,intent(in) :: n, m
        double precision,intent(in) :: dh, xl
        double complex,intent(in) :: phi(1:n, 1:m)
        double precision solve_x(1:m)
        double precision x
        integer i, j
        solve_x = 0d0

        do j = 1, m
            do i = 1, n
                x = xl + dh * i
                if (i == 1 .or. i == n) then
                    solve_x(j) = solve_x(j) + 0.5d0*x*abs(phi(i, j))**2d0*dh
                else
                    solve_x(j) = solve_x(j) + x*abs(phi(i, j))**2d0*dh
                end if
            end do
        end do

    end function

    function solve_p(phi, n, m)
        integer,intent(in) :: n, m
        double complex,intent(in) :: phi(1:n, 1:m)
        double complex solve_p(1:m)
        integer i, j
        solve_p = complex(0d0, 0d0)

        do j = 1, m
            do i = 1, n
                if (i == 1) then
                    solve_p(j) = solve_p(j) - 0.5d0*HBAR*conjg(phi(1,j))*ix_scaler(phi(2,j)-phi(1,j))
                else if (i == n) then
                    solve_p(j) = solve_p(j) - 0.5d0*HBAR*conjg(phi(n,j))*ix_scaler(phi(n,j)-phi(n-1,j))
                else
                    solve_p(j) = solve_p(j) - HBAR*conjg(phi(i,j))*ix_scaler(phi(i+1,j)-phi(i,j))
                end if
            end do
        end do

    end function
end module

program main
    use extensions
    implicit none
    integer,parameter :: n = 1000, m = 4000+1
    double precision,parameter :: dh = 0.05d0, dt = 0.001d0, xl = -dh*(n/2)
    double precision t1, t2
    double complex phi(1:n, 1:m), p(1:m)
    double precision x(1:m)
    integer j
    phi = complex(0d0, 0d0)

    write (*, *) "Start Calculation..."
    call cpu_time(t1)
    call solve_schroedinger(phi, n, m, dh, dt, xl)
    call cpu_time(t2)
    write (*, *) t2 - t1, " seconds"

    write (*, *) "Start Plotting..."
    call plot(phi, n, m, dh, xl)

    write (*, *) "Time dependance of probability"
    do j = 1, m, 1000
        write (*, *) "j=", j, calc_probability(phi(:,j), n, dh)
    end do

    write (*, *) "Time dependance of expected value of x"
    x = solve_x(phi, n, m, dh, xl)
    do j = 1, m, 1000
        write (*, *) "j=", j, x(j)
    end do

    write (*, *) "Time dependance of expected value of p"
    p = solve_p(phi, n, m)
    do j = 1, m, 1000
        write (*, *) "j=", j, p(j)
    end do
end program 