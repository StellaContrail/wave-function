module extensions
    implicit none
    double precision,parameter :: HBAR = 1d0
    double precision,parameter :: MASS = 0.5d0
    double precision,parameter :: OMEGA = 1d0
    double precision,parameter :: K = 0d0 ! wave number
    double precision,parameter :: a = HBAR / (MASS * OMEGA)
    double precision,parameter :: x0 = 0d0 ! initial center point of Gaussian wave packet
contains
    ! i*(a+ib)=-b+ia : (REAL=-b), (IMAG=a)
    function ix(A_matrix, n)
        integer,intent(in) :: n
        double complex,intent(in) :: A_matrix(1:n)
        double complex ix(1:n)
        integer i
        do i = 1, n
            ix(i) = dcmplx(-aimag(A_matrix(i)), dble(A_matrix(i)))
        end do
    end function
    function ix_scaler(z)
        double complex,intent(in) :: z
        double complex ix_scaler
        ix_scaler = dcmplx(-aimag(z), dble(z))
    end function

    subroutine solve_schroedinger(phi, n, m, dh, dt, xl, H)
        integer,intent(in) :: n, m
        double complex,intent(inout) :: phi(1:n, 1:m)
        double precision,intent(in) :: dh, dt, xl
        double precision,intent(out) ::  H(1:n, 1:n)
        double precision t1, t2, speed
        integer j
        H = 0d0
        open(10, file="data.txt")
        call construct_hamiltonian(H, n, xl, dh)

        phi = dcmplx(0d0, 0d0)
        call initialize(phi(:, 1), n, dh, xl)
        call normalize(phi(:, 1), n)

        call cpu_time(t1)
        do j = 1, m-1
            phi(:, j+1) = calc_future_wavefunction(phi(:,j), H, n, dt)

            if (mod(j, 100) == 0) then
                call cpu_time(t2)
                speed = 100d0/(t2-t1)
                write (*, 100, advance='no') (100d0*j)/m, " % ", speed, " items/sec", " ETA : ", (m-1-j)/speed, " sec"
                100 format(F6.2, A, F15.7, 2A, F10.5, A)
                write (*, *)
                call cpu_time(t1)
            endif
        end do
        close(10)
    end subroutine

    subroutine normalize(phi, n)
        integer,intent(in) :: n
        double complex,intent(inout) :: phi(1:n)
        double complex prob
        prob = calc_probability(phi, n)
        if (aimag(prob) > 1d-15) then
            stop "Probability calculation failed"
        end if
        phi = phi / sqrt(prob)
    end subroutine

    ! Calculate total probability of the system
    double complex function calc_probability(phi, n)
        integer,intent(in) :: n
        double complex,intent(in) :: phi(1:n)
        calc_probability = dot_product(phi, phi)
    end function

    subroutine initialize(phi, n, dh, xl)
        integer,intent(in) :: n
        double complex,intent(out) :: phi(1:n)
        double precision,intent(in) :: dh, xl
        integer i
        double precision x

        do i = 1, n
            x = xl + dh * i
            phi(i) = dcmplx(cos(K*x)*f(x), sin(K*x)*f(x))
        end do  
    end subroutine

    double precision function f(x)
        double precision,intent(in) :: x
        f = exp(-0.5d0*(x-x0)**2d0/a)
    end function

    function calc_future_wavefunction(phi, H, n, dt)
        integer,intent(in) :: n
        double complex,intent(in) :: phi(1:n)
        double precision,intent(in) :: dt, H(1:n, 1:n)
        double complex :: calc_future_wavefunction(1:n)
        double complex :: phi_next(1:n), phi_temp(1:n)
        integer i

        phi_next = dcmplx(0d0, 0d0)

        ! n = 0
        phi_temp = phi
        phi_next = phi_next + phi_temp

        i = 1
        do
            phi_temp = -ix(symmatmul(H, phi_temp, n), n)*dt / (HBAR*dble(i))
            
            phi_next = phi_next + phi_temp
            ! Condition to tell if calculation converges
            if (sum(abs(dble(phi_temp))) < 1d-10 .and. sum(abs(aimag(phi_temp))) < 1d-10) then
                exit
            end if
            i = i + 1
        end do

        calc_future_wavefunction = phi_next
    end function

    ! H : Double precision matrix
    ! B : Double complex matrix
    ! C : Double complex matrix
    function symmatmul(H, B, n) result(C)
        integer,intent(in) :: n
        double precision,intent(in) :: H(1:n, 1:n)
        double complex,intent(in) :: B(1:n)
        double complex C(1:n)
        integer i
        C = (0d0, 0d0)
        
        do i = 1, n
            C(i) = H(i,i)*B(i)
            if (i > 1) then
                C(i) = C(i) + H(i,i-1)*B(i-1)
            end if
            if (i > 2) then
                C(i) = C(i) + H(i,i-2)*B(i-2)
            end if
            if (i > 3) then
                C(i) = C(i) + H(i,i-3)*B(i-3)
            end if
            if (i > 4) then
                C(i) = C(i) + H(i,i-4)*B(i-4)
            end if
            if (i > 5) then
                C(i) = C(i) + H(i,i-5)*B(i-5)
            end if
            if (i < n-1) then
                C(i) = C(i) + H(i, i+1)*B(i+1)
            end if
            if (i < n-2) then
                C(i) = C(i) + H(i, i+2)*B(i+2)
            end if
            if (i < n-3) then
                C(i) = C(i) + H(i, i+3)*B(i+3)
            end if
            if (i < n-4) then
                C(i) = C(i) + H(i, i+4)*B(i+4)
            end if
            if (i < n-5) then
                C(i) = C(i) + H(i, i+5)*B(i+5)
            end if
        end do
    end function

    subroutine construct_hamiltonian(H, n, xl, dh)
        double precision,parameter :: ALPHA = 0.5d0*HBAR*HBAR/MASS ! Coefficient of Kinetic Energy
        integer,intent(in) :: n
        double precision,intent(out) :: H(1:n, 1:n)
        double precision,intent(in)  :: xl, dh
        integer i
        double precision x

        H = 0d0
        do i = 1, n
            H(i, i) = -73766d0
            if (i > 1) then
                H(i, i-1) = 42000d0
            end if
            if (i > 2) then
                H(i, i-2) = -6000d0
            end if
            if (i > 3) then
                H(i, i-3) = 1000d0
            end if
            if (i > 4) then
                H(i, i-4) = -125d0
            end if
            if (i > 5) then
                H(i, i-5) = 8d0
            end if

            if (i < n-1) then
                H(i, i+1) = 42000d0
            end if
            if (i < n-2) then
                H(i, i+2) = -6000d0
            end if
            if (i < n-3) then
                H(i, i+3) = 1000d0
            end if
            if (i < n-4) then
                H(i, i+4) = -125d0
            end if
            if (i < n-5) then
                H(i, i+5) = 8d0
            end if

        end do
        H = H / (25200d0*dh*dh)
        H = - ALPHA * H

        do i = 1, n
            x = xl + dh * i
            H(i, i) = H(i, i) + V(x)
        end do
    end subroutine
    double precision function V(x)
        double precision,parameter :: BETA = MASS * OMEGA * OMEGA
        double precision,optional,intent(in) :: x
        V = 0.5d0*x*x*BETA
        
        ! Well potential
        !if (abs(x) < 5d0) then
        !    V = -5d0
        !else
        !    V = 0d0
        !end if
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
                write (10, '(4(F15.10, X))', advance='no') x, dble(phi(i, j)), aimag(phi(i, j)), abs(phi(i,j))
                write (10, *)
            end do
            write (10, *)
        end do

        close(10)
    end subroutine
    
    function solve_x(phi, n, m, dh, xl)
        integer,intent(in) :: n, m
        double precision,intent(in) :: dh, xl
        double complex,intent(in) :: phi(1:n, 1:m)
        double complex solve_x(1:m)
        double precision x
        integer i, j
        solve_x = 0d0

        do j = 1, m
            do i = 1, n
                x = xl + dh * i
                solve_x(j) = solve_x(j) + x * abs(phi(i, j))**2d0
            end do
        end do

    end function

    function solve_p(phi, n, m, dh)
        integer,intent(in) :: n, m
        double complex,intent(in) :: phi(1:n, 1:m)
        double complex solve_p(1:m), temp(1:n)
        double precision,intent(in) :: dh
        double precision diff_op(1:n, 1:n)
        integer i, j
        solve_p = dcmplx(0d0, 0d0)

        diff_op = 0d0
        diff_op(1,1) = -1d0
        diff_op(1,2) = 1d0
        diff_op(n,n-1) = -1d0
        diff_op(n,n) = 1d0
        do i = 2, n-1
            diff_op(i, i+1) = 0.5d0
            diff_op(i, i-1) = -0.5d0
        end do
        diff_op = diff_op / dh
        
        do j = 1, m
            temp(1) = diff_op(1,1)*phi(1,j)+diff_op(1,2)*phi(2,j)
            temp(n) = diff_op(n,n)*phi(n,j)+diff_op(n,n-1)*phi(n-1,j)
            do i = 2, n-1
                temp(i) = diff_op(i,i-1)*phi(i-1,j)+diff_op(i,i+1)*phi(i+1,j)
            end do
            solve_p(j) = -HBAR*ix_scaler(dot_product(phi(:, j), temp))
        end do
    end function

    subroutine solve_energy(phi, n, m, H, energies)
        integer,intent(in) :: n, m
        double complex,intent(in) :: phi(1:n, 1:m)
        double precision,intent(in) :: H(1:n, 1:n)
        double complex,intent(out) :: energies(1:m)
        double complex :: temp(1:n, 1:m)
        integer j
        do j = 1, m
            temp(:, j) = symmatmul(H, phi(:, j), n)
            energies(j) = dot_product(phi(:, j), temp(:, j))
        end do
    end subroutine
end module

program main
    use extensions
    implicit none
    integer,parameter :: n = 40, m = 630+1, m_display_step = 100
    double precision,parameter :: dh = 0.8d0, dt = 0.02d0, xl = -dh*(n/2)
    double precision t1, t2
    double complex phi(1:n, 1:m), p(1:m), E(1:m), x(1:m)
    double precision H(1:n, 1:n)
    integer j
    phi = dcmplx(0d0, 0d0)

    write (*, *) "Start Calculation..."
    call cpu_time(t1)
    call solve_schroedinger(phi, n, m, dh, dt, xl, H)
    call cpu_time(t2)
    write (*, '(A, F15.10, A)') "Calculation Time : ", t2 - t1, " seconds"

    write (*, *) "Start Plotting..."
    call plot(phi, n, m, dh, xl)
    write (*, *) "Plotting end"

    write (*, *) "Time dependance of probability"
    do j = 1, m, m_display_step
        write (*, *) "j=", j, dble(calc_probability(phi(:,j), n))
    end do

    write (*, *) "Time dependance of expected value of x"
    x = solve_x(phi, n, m, dh, xl)
    do j = 1, m, m_display_step
        write (*, '(A, I10, 2(A, F15.10))') "j=", j, " REAL=", dble(x(j)), " IMAG=", aimag(x(j))
    end do

    write (*, *) "Time dependance of expected value of p"
    p = solve_p(phi, n, m, dh)
    do j = 1, m, m_display_step
        write (*, '(A, I10, 2(A, F15.10))') "j=", j, " REAL=", dble(p(j)), " IMAG=", aimag(p(j))
    end do

    write (*, *) "Time dependance of expected energy E"
    call solve_energy(phi, n, m, H, E)
    do j = 1, m, m_display_step
        write (*, '(A, I10, 2(A, F15.10))') "j=", j, " REAL=", dble(E(j)), " IMAG=", aimag(E(j))
    end do
end program 
