! Three dimensional harmonic oscillator with time development
module extension
    implicit none
    double precision,parameter :: HBAR = 1d0    ! Plank constant
    double precision,parameter :: MASS = 0.5d0  ! mass of the harmonic oscillator
    double precision,parameter :: OMEGA = 1d0   ! angular frequency of the potential
    ! squared covariance of initial Gaussian wavepacket
    double precision,parameter :: a = HBAR / (MASS * OMEGA)
contains
    subroutine solve(n, m, dh, dt, xl)
        integer,intent(in) :: n, m
        double precision,intent(in) :: dh, dt, xl
        double complex phi(1:n, 1:n, 1:n)
        integer i
        call initialize(phi, n, dh, xl)
        !call normalize(phi, n, dh)
        open(10, file="data.txt")
        call plot(phi, n, dh, xl, 10)
        write (*, '(A, I5)') "i=", 0
        write (*, '(A, F15.10)') "P=", probability(phi, n, dh)
        write (*, '(A, 2F15.10)') "E=", energy(phi, n, dh, xl)
        write (*, '(A, 2F15.10)') "x=", expect_x(phi, n, dh, xl)
        write (*, *)

        do i = 1, m
            phi = evolve(phi, n, dh, dt, xl)
            if (mod(i, 25) == 0d0) then
                call plot(phi, n, dh, xl, 10)
                write (*, '(A, I5)') "i=", i
                write (*, '(A, F15.10)') "P=", probability(phi, n, dh)
                write (*, '(A, 2F15.10)') "E=", energy(phi, n, dh, xl)
                write (*, '(A, 2F15.10)') "x=", expect_x(phi, n, dh, xl)
                write (*, *)
            end if
        end do
        close(10)
    end subroutine
    subroutine initialize(phi, n, dh, xl)
        integer,intent(in) :: n
        double complex,intent(out) :: phi(1:n, 1:n, 1:n)
        double precision,intent(in) :: dh, xl
        integer i, j, k
        double precision x, y, z

        do k = 1, n
            z = xl + dh * k
            z = z - 2.3d0
            do j = 1, n
                y = xl + dh * j
                y = y - 2.3d0
                do i = 1, n
                    x = xl + dh * i
                    x = x - 2.3d0
                    phi(i, j, k) = (a*acos(-1d0))**(-0.75d0)*exp(-0.5d0*(x*x+y*y+z*z)/a)
                end do
            end do
        end do
    end subroutine
    function probability(phi, n, dh) result(prob)
        integer,intent(in) :: n
        double precision,intent(in) :: dh
        double complex,intent(inout) :: phi(1:n, 1:n, 1:n)
        double precision temp_1(1:n, 1:n), temp_2(1:n), sum, prob
        integer i, j, k
        sum = 0d0
        temp_1 = 0d0
        temp_2 = 0d0
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    if (i == 1 .or. i == n) then
                        temp_1(j, k) = temp_1(j, k) + 0.5d0*abs(phi(i, j, k))**2d0
                    else
                        temp_1(j, k) = temp_1(j, k) + abs(phi(i, j, k))**2d0
                    end if
                end do
            end do
        end do
        temp_1 = temp_1 * dh
        do k = 1, n
            do j = 1, n
                if (j == 1 .or. j == n) then
                    temp_2(k) = temp_2(k) + 0.5d0*temp_1(j, k)
                else
                    temp_2(k) = temp_2(k) + temp_1(j, k)
                end if
            end do
        end do
        temp_2 = temp_2 * dh
        do k = 1, n
            if (k == 1 .or. k == n) then
                sum = sum + 0.5d0*temp_2(k)
            else
                sum = sum + temp_2(k)
            end if
        end do
        prob  = sum * dh
    end function
    subroutine plot(phi, n, dh, xl, unit)
        integer,intent(in) :: n, unit
        double precision,intent(in) :: dh, xl
        double complex,intent(in) :: phi(1:n, 1:n, 1:n)
        integer i, j, k
        double precision x, y, z
        double complex density(1:n, 1:n)
        do j = 1, n
            do k = 1, n
                do i = 1, n
                    if (k == 1 .or. k == n) then
                        density(i,j) = density(i,j) + 0.5d0*abs(phi(i,j,k))**2d0 
                    else
                        density(i,j) = density(i,j) + abs(phi(i,j,k))**2d0
                    end if
                end do
            end do
        end do
        density = density * dh

        do j = 1, n
            y = xl + dh * j
            do i = 1, n
                x = xl + dh * i
                write (unit,'(3F15.10)',advance='no') x, y, density(i,j)
                write (unit, *)
            end do
            write (unit, *)
        end do
    end subroutine
    function evolve(phi, n, dh, dt, xl) result(phi_new)
        integer,intent(in) :: n
        double complex,intent(in) :: phi(1:n, 1:n, 1:n)
        double precision,intent(in) :: dh, dt, xl
        double complex phi_new(1:n, 1:n, 1:n), phi_temp(1:n, 1:n, 1:n)
        integer i

        ! n = 0
        phi_temp = phi
        phi_new = phi_temp

        ! n = 1 ~ 4
        do i = 1, 10
            phi_temp = H(phi_temp, n, dh, xl)*dt/(HBAR*i)
            phi_temp = -ix(phi_temp, n)
            phi_new = phi_new + phi_temp
        end do
    end function
    function H(phi, n, dh, xl) result(HPhi)
        integer,intent(in) :: n
        double complex,intent(in) :: phi(1:n, 1:n, 1:n)
        double precision,intent(in) :: dh, xl
        double complex HPhi(1:n, 1:n, 1:n)
        integer i, j, k
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    HPhi(i, j, k) = -30d0*3d0*phi(i, j, k)
                    if (i > 1) then
                        HPhi(i, j, k) = HPhi(i, j, k) + 16d0*phi(i-1, j, k)
                    end if
                    if (i > 2) then
                        HPhi(i, j, k) = HPhi(i, j, k) - phi(i-2, j, k)
                    end if
                    if (i < n-1) then
                        HPhi(i, j, k) = HPhi(i, j, k) - phi(i+2, j, k)
                    end if
                    if (i < n) then
                        HPhi(i, j, k) = HPhi(i, j, k) + 16d0*phi(i+1, j, k)
                    end if

                    if (j > 1) then
                        HPhi(i, j, k) = HPhi(i, j, k) + 16d0*phi(i, j-1, k)
                    end if
                    if (j > 2) then
                        HPhi(i, j, k) = HPhi(i, j, k) - phi(i, j-2, k)
                    end if
                    if (j < n-1) then
                        HPhi(i, j, k) = HPhi(i, j, k) - phi(i, j+2, k)
                    end if
                    if (j < n) then
                        HPhi(i, j, k) = HPhi(i, j, k) + 16d0*phi(i, j+1, k)
                    end if

                    if (k > 1) then
                        HPhi(i, j, k) = HPhi(i, j, k) + 16d0*phi(i, j, k-1)
                    end if
                    if (k > 2) then
                        HPhi(i, j, k) = HPhi(i, j, k) - phi(i, j, k-2)
                    end if
                    if (k < n-1) then
                        HPhi(i, j, k) = HPhi(i, j, k) - phi(i, j, k+2)
                    end if
                    if (k < n) then
                        HPhi(i, j, k) = HPhi(i, j, k) + 16d0*phi(i, j, k+1)
                    end if
                end do
            end do
        end do
        HPhi = -0.5d0*HBAR*HBAR*HPhi/(12d0*MASS*dh*dh)

        do k = 1, n
            do j = 1, n
                do i = 1, n
                    HPhi(i, j, k) = HPhi(i, j, k) + V(dh*i+xl, dh*j+xl, dh*k+xl)*phi(i, j, k)
                end do
            end do
        end do
    end function
    function V(x, y, z) result(potential)
        double precision,intent(in) :: x, y, z
        double precision potential
        double precision,parameter :: BETA = MASS * OMEGA * OMEGA
        double precision r2, r
        r2 = x*x+y*y+z*z
        r = sqrt(r2)
        !potential = 0.5d0*r2*BETA
        !potential = -5d0*exp(-r2)
        !if (abs(x) < 5d0 .and. abs(y) < 5d0 .and. abs(z) < 5d0) then
        !    potential = -50d0
        !else
        !    potential = 0d0
        !end if
        if (2d0 < r .and. r < 6d0) then
            potential = -50d0
        else 
            potential = 0d0
        end if
        !potential = 0d0
    end function
    function ix(matrix, n) result(output_matrix)
        integer,intent(in) :: n
        double complex,intent(in) :: matrix(1:n, 1:n, 1:n)
        double complex output_matrix(1:n, 1:n, 1:n)
        integer i, j, k

        do k = 1, n
            do j = 1, n
                do i = 1, n
                    output_matrix(i,j,k) = dcmplx(-aimag(matrix(i,j,k)), dble(matrix(i,j,k)))
                end do
            end do
        end do
    end function
    function energy(phi, n, dh, xl) result(E)
        integer,intent(in) :: n
        double precision,intent(in) :: dh, xl
        double complex,intent(in) :: phi(1:n, 1:n, 1:n)
        double complex E, HPhi(1:n, 1:n, 1:n)
        double complex temp_1(1:n, 1:n), temp_2(1:n), sum
        integer i,j,k
        sum = 0d0
        temp_1 = 0d0
        temp_2 = 0d0
        HPhi = H(phi, n, dh, xl)
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    if (i == 1 .or. i == n) then
                        temp_1(j,k) = temp_1(j,k) + 0.5d0*conjg(phi(i,j,k))*HPhi(i,j,k)
                    else 
                        temp_1(j,k) = temp_1(j,k) + conjg(phi(i,j,k))*HPhi(i,j,k)
                    end if
                end do
            end do
        end do
        temp_1 = temp_1 * dh
        do k = 1, n
            do j = 1, n
                if (j == 1 .or. j == n) then
                    temp_2(k) = temp_2(k) + 0.5d0*temp_1(j, k)
                else
                    temp_2(k) = temp_2(k) + temp_1(j, k)
                end if
            end do
        end do
        temp_2 = temp_2 * dh
        do k = 1, n
            if (k == 1 .or. k == n) then
                sum = sum + 0.5d0*temp_2(k)
            else
                sum = sum + temp_2(k)
            end if
        end do
        E  = sum * dh
    end function
    function expect_x(phi, n, dh, xl) result(x_value)
        integer,intent(in) :: n
        double precision,intent(in) :: dh, xl
        double complex,intent(in) :: phi(1:n, 1:n, 1:n)
        double complex x_value
        double complex temp_1(1:n, 1:n), temp_2(1:n), sum
        integer i,j,k
        double precision x
        sum = 0d0
        temp_1 = 0d0
        temp_2 = 0d0
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    x = xl + dh * i
                    if (i == 1 .or. i == n) then
                        temp_1(j,k) = temp_1(j,k) + 0.5d0*conjg(phi(i,j,k))*x*phi(i,j,k)
                    else 
                        temp_1(j,k) = temp_1(j,k) + conjg(phi(i,j,k))*x*phi(i,j,k)
                    end if
                end do
            end do
        end do
        temp_1 = temp_1 * dh
        do k = 1, n
            do j = 1, n
                if (j == 1 .or. j == n) then
                    temp_2(k) = temp_2(k) + 0.5d0*temp_1(j, k)
                else
                    temp_2(k) = temp_2(k) + temp_1(j, k)
                end if
            end do
        end do
        temp_2 = temp_2 * dh
        do k = 1, n
            if (k == 1 .or. k == n) then
                sum = sum + 0.5d0*temp_2(k)
            else
                sum = sum + temp_2(k)
            end if
        end do
        x_value  = sum * dh
    end function
end module

program main
    use extension
    implicit none
    integer,parameter :: n = 40, m = 1200
    double precision,parameter :: dh = 0.8d0, dt = 0.02d0
    ! x = xl ~ xl+dh*n
    double precision,parameter :: xl = -dh*(n/2)
    double precision t1, t2

    write (*, *) "Start Calculation..."
    call cpu_time(t1)
    call solve(n, m, dh, dt, xl)
    call cpu_time(t2)
    write (*, '(A, F15.10, A)') "Calculation Time : ", t2 - t1, " seconds"
end program
