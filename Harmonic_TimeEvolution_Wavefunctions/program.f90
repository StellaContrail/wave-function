module extensions
    implicit none
contains
    subroutine initialize(phi, n, xl, xu, x0, a, dh)
        double complex, intent(out) :: phi(0:, 0:)
        double precision, intent(in) :: xl, xu, x0, a, dh
        integer, intent(in) :: n
        integer i
        double precision x
        x = 0d0
        ! Initialize wave functions with zero
        phi = (0d0, 0d0)
        ! Initialize at t=0 wave function as Gaussian Wave packet
        do i = 0, n
            x = xl + dh * i
            phi(i, 0) = dcmplx(exp(-0.5d0*a*a*(x-x0)**2d0), 0d0)
        end do
        call normalize(phi(:, 0), n, dh)
        ! Set Boundary Condition
        phi(0, :) = (0d0, 0d0)
        phi(n, :) = (0d0, 0d0)
    end subroutine

    subroutine normalize(phi, n, dh)
        integer,intent(in) :: n
        double complex,intent(inout) :: phi(:)
        double precision,intent(in) :: dh
        double precision sum
        integer i
        sum = 0d0
        do i = 0, n
            if (i /= 0 .or. i /= n) then
                sum = sum + norm_2(phi(i))*dh
            else
                sum = sum + 0.5d0 * norm_2(phi(i))*dh
            end if
        end do
        phi(:) = phi(:) / sqrt(sum)
    end subroutine  

    subroutine solve_schroedinger(phi, n, m, xl, xu, dh, dt)
        double complex, intent(inout) :: phi(0:, 0:)
        double precision, intent(in) :: xl, xu, dh, dt
        integer, intent(in) :: n, m
        double complex k1(0:n), k2(0:n), k3(0:n), k4(0:n)
        double precision x
        integer i, j

        ! j : Time
        do j = 0, m-1
            ! i : Coordinates
            do i = 1, n-1
                x = xl + i * dh
                k1(i) = -conjg(H(phi(i-1, j), phi(i, j), phi(i+1,j), x, dh))*dt
            end do
            k1(0) = (0d0, 0d0)
            k1(n) = (0d0, 0d0)

            do i = 1, n-1
                x = xl + i * dh
                k2(i) = -conjg(H(phi(i-1, j)+0.5d0*k1(i-1)*dt, phi(i, j)+0.5d0*k1(i)*dt, phi(i+1,j)+0.5d0*k1(i+1)*dt, x, dh))*dt
            end do
            k2(0) = (0d0, 0d0)
            k2(n) = (0d0, 0d0)

            do i = 1, n-1
                x = xl + i * dh
                k3(i) = -conjg(H(phi(i-1, j)+0.5d0*k2(i-1)*dt, phi(i, j)+0.5d0*k2(i)*dt, phi(i+1,j)+0.5d0*k2(i+1)*dt, x, dh))*dt
            end do
            k3(0) = (0d0, 0d0)
            k3(n) = (0d0, 0d0)
            
            do i = 1, n-1
                x = xl + i * dh
                k4(i) = -conjg(H(phi(i-1, j)+k3(i-1)*dt, phi(i, j)+k3(i)*dt, phi(i+1,j)+k3(i+1)*dt, x, dh))*dt
            end do

            do i = 1, n-1
                phi(i, j+1) = phi(i, j) + (k1(i) + 2d0*k2(i) + 2d0*k3(i) + k4(i)) / 6d0
            end do
        end do
    end subroutine

    double complex function H(phi_back, phi_center, phi_forward, x, dh)
        double complex,intent(in) :: phi_back, phi_center, phi_forward
        double precision,intent(in) :: dh, x
        double precision re, im

        re = -0.5d0*(dble(phi_forward)-2d0*dble(phi_center)+dble(phi_back))/(dh*dh) + 0.5d0*x*x*dble(phi_center)
        im = -0.5d0*(imag(phi_forward)-2d0*imag(phi_center)+imag(phi_back))/(dh*dh) + 0.5d0*x*x*imag(phi_center)
        H = dcmplx(re, im)
    end function

    subroutine total_probability(phi, prob, dh, n, m)
        double complex,intent(in) :: phi(0:, 0:)
        double precision,intent(in) :: dh
        double precision,intent(out) :: prob(0:)
        integer,intent(in) ::n,  m
        double precision sum
        integer i, j

        do j = 0, m
            sum = 0d0
            do i = 0, n
                if (i /= 0 .or. i /= n) then
                    sum = sum + norm_2(phi(i, j)) * dh
                else
                    sum = sum + 0.5d0 * norm_2(phi(i, j)) * dh
                end if
            end do
            prob(j) = sum
        end do
    end subroutine

    double precision function norm_2(cmp)
        double complex,intent(in) :: cmp
        norm_2 = dble(cmp) * dble(cmp) + imag(cmp) * imag(cmp)
    end function    

    subroutine output(phi, mode, i, j, xl, xu, n, m, dh, dt)
        double complex,intent(in) :: phi(0:, 0:)
        double precision,intent(in) :: xl, xu, dt, dh
        character, intent(in) :: mode
        integer,intent(in) :: i, j, n, m
        integer k, l
        integer,parameter :: fo = 10
        double precision time
        open(fo, file="data.txt")

        if (mode == "X") then
            ! Specify one time point, output wave amplitudes at all positions
            do l = 0, n
                write (fo, *) xl + dh*l, dble(phi(l, j)), imag(phi(l, j))
            end do
        else if (mode == "T") then
            do k = 0, m
                write (fo, *) dt*k, dble(phi(i, k)), imag(phi(i, k))
            end do
        else if (mode == "A") then
            do k = 0, m
                time = dt * k
                do l = 0, n
                    write (fo, *) time, xl + dh*l, dble(phi(l, k)), imag(phi(l, k))
                end do
                write (fo, *)
            end do
        else
            close(fo)
            stop "Specify one type of modes"
        end if

        close(fo)
    end subroutine  
end module 

program main
    use extensions
    implicit none
    integer,parameter          :: n  = 200, m  = 2000
    double precision,parameter :: dh = 0.2d0, dt = 0.01d0
    double precision :: xl = -20d0, xu = 20d0, x0 = 0d0, a = 0.5d0, prob(0:m)
    double complex :: phi(0:n, 0:m) = (0d0, 0d0)

    call initialize(phi, n, xl, xu, x0, a, dh)
    call solve_schroedinger(phi, n, m, xl, xu, dh, dt)
    call output(phi, "A", 0, 0, xl, xu, n, m, dh, dt)
    call total_probability(phi, prob, dh, n, m)
    write (*, *) prob(0), prob(1), prob(10), prob(100), prob(1000), prob(2000)
end program