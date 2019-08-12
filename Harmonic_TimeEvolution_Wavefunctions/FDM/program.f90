module extensions
    implicit none
    double precision,parameter :: ALPHA = 1d0 ! (HBAR*C)^2/(2*MC^2)
    double precision,parameter :: BETA = 1d0 ! (MASS*OMEGA)^2
    double precision,parameter :: HBAR = 1d0
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

    subroutine solve_schroedinger(phi, n, m, dh, dt, xl)
        integer,intent(in) :: n, m
        double complex,intent(inout) :: phi(1:n, 1:m)
        double precision,intent(in) :: dh, dt, xl
        double precision H(1:n, 1:n), x
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

        do i = 1, m-1
            phi(:, i+1) = calc_future_wavefunction(phi(:,i), H, n, dt)
        end do
    end subroutine

    subroutine normalize(phi, n, dh)
        integer,intent(in) :: n
        double complex,intent(inout) :: phi(1:n)
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
            phi(i) = complex(exp(-0.5d0*(x/2d0)**2d0), 0d0)
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
                write (10, *) x, dble(phi(i, j)), aimag(phi(i, j))
            end do
            write (10, *)
        end do

        close(10)
    end subroutine
end module

program main
    use extensions
    implicit none
    integer,parameter :: n = 400, m = 2000
    double precision,parameter :: dh = 0.1d0, dt = 0.0001d0, xl = -dh*(n/2)
    double complex phi(1:n, 1:m)
    phi = complex(0d0, 0d0)

    call solve_schroedinger(phi, n, m, dh, dt, xl)

    call plot(phi, n, m, dh, xl)

    write (*, *) "j=1", calc_probability(phi(:,1), n, dh)
    write (*, *) "j=10", calc_probability(phi(:,10), n, dh)
    write (*, *) "j=100", calc_probability(phi(:,100), n, dh)
    write (*, *) "j=200", calc_probability(phi(:,200), n, dh)
    write (*, *) "j=500", calc_probability(phi(:,500), n, dh)
end program 