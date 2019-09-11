! Three dimensional harmonic oscillator with time development
module extension
    implicit none
    double precision,parameter :: HBAR = 1d0    ! Plank constant
    double precision,parameter :: MASS = 0.5d0  ! mass of the harmonic oscillator
    double precision,parameter :: OMEGA = 1d0   ! angular frequency of the potential
    ! wave number which wave function initially has
    double precision,parameter :: wavenum_x = 0d0
    double precision,parameter :: wavenum_y = 0d0
    double precision,parameter :: wavenum_z = 0d0
    ! initial center point of Gaussian wave packet
    double precision,parameter :: x0 = 0d0
    double precision,parameter :: y0 = 0d0
    double precision,parameter :: z0 = 0d0
    ! squared covariance of initial Gaussian wavepacket
    double precision,parameter :: a = HBAR / (MASS * OMEGA)
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

    ! Solve all space and time steps of harmonic oscillator
    subroutine solve(n, m, dh, dt, xl, H)
        integer,intent(in) :: n, m
        double precision,intent(in) :: dh, dt, xl
        double precision,intent(out) ::  H(1:n, 1:n)
        double complex :: phi_next(1:n, 1:n, 1:n), phi_old(1:n, 1:n, 1:n)
        double precision t1, t2, speed, x, y, z
        integer i, j, k, l
        H = 0d0
        phi_old = dcmplx(0d0, 0d0)
        open(10, file="data.txt")

        call construct_hamiltonian(H, n, xl, dh)
        call initialize(phi_old, n, dh, xl)
        call normalize(phi_old, n, dh)

        do l = 1, 10
            phi_next = (0d0, 0d0)
            do k = 1, n
                do j = 1, n
                    phi_next(:,j,k) = phi_next(:,j,k) + calc_next(phi_old(:,j,k), H, n, dt)
                end do
            end do

            do k = 1, n
                do i = 1, n
                    phi_next(i,:,k) = phi_next(:,j,k) + calc_next(phi_old(i,:,k), H, n, dt)
                end do
            end do

            do j = 1, n
                do i = 1, n
                    phi_next(i,j,:) = phi_next(i,j,:) + calc_next(phi_old(i,j,:), H, n, dt)
                end do
            end do

            do j = 1, n
                z = xl + dh*j
                do k = 1, n
                    y = xl + dh*k
                    do i = 1, n
                        x = xl + dh*i
                        write (10, '(4F15.10)', advance='no') x, y, z, abs(phi_next(i,j,k))
                        write (10, *)
                    end do
                end do
            end do
            write (10, *)

            phi_old = phi_next
            write (*, *) l, " end"
            write (*, *) "Prob = ", calc_probability(phi_old, n, dh)
        end do
        
        close(10)
    end subroutine

    ! Construct the Hamiltonian matrix with 11-points stencil
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

    ! Potential
    double precision function V(x)
        double precision,parameter :: BETA = MASS * OMEGA * OMEGA
        double precision,optional,intent(in) :: x
        V = 0.5d0*x*x*BETA
    end function

    ! Initialize wave function with given function f(r)
    subroutine initialize(phi, n, dh, xl)
        integer,intent(in) :: n
        double complex,intent(out) :: phi(1:n, 1:n, 1:n)
        double precision,intent(in) :: dh, xl
        integer i, j, k
        double precision x, y, z, k_r

        do k = 1, n
            z = xl + dh * k
            do j = 1, n
                y = xl + dh * j
                do i = 1, n
                    x = xl + dh * i
                    k_r = wavenum_x*x+wavenum_y*y+wavenum_z*z
                    phi(i, j, k) = dcmplx(cos(k_r)*f(x,y,z), sin(k_r)*f(x,y,z))
                end do
            end do
        end do  
    end subroutine
    double precision function f(x, y, z)
        double precision,intent(in) :: x, y, z
        double precision r_r0_2
        r_r0_2 = (x-x0)**2d0+(y-y0)**2d0+(z-z0)**2d0
        f = exp(-0.5d0*r_r0_2/a)
    end function

    ! Calculate total probability of the system
    ! TODO: Try summation version
    double complex function calc_probability(phi, n, dh)
        integer,intent(in) :: n
        double complex,intent(in) :: phi(1:n, 1:n, 1:n)
        double precision,intent(in) :: dh
        integer i, j, k
        calc_probability = (0d0, 0d0)

        do k = 1, n
            do j = 1, n
                do i = 1, n
                    if (i == 1 .or. i == n .or. j == 1 .or. j == n .or. k == 1 .or. k == n) then
                        calc_probability = calc_probability + 0.5d0*abs(phi(i,j,k))**2d0*dh
                    else
                        calc_probability = calc_probability + abs(phi(i,j,k))**2d0*dh
                    end if
                end do
            end do
        end do
    end function

    ! Normalization
    subroutine normalize(phi, n, dh)
        integer,intent(in) :: n
        double complex,intent(inout) :: phi(1:n, 1:n, 1:n)
        double precision,intent(in) :: dh
        double complex prob
        prob = calc_probability(phi, n, dh)
        if (aimag(prob) > 1d-15) then
            stop "Probability calculation failed"
        end if
        phi = phi / sqrt(prob)
    end subroutine

    
    function calc_next(phi, H, n, dt) result(phi_next)
        integer,intent(in) :: n
        double complex,intent(in) :: phi(1:n)
        double precision,intent(in) :: dt, H(1:n, 1:n)
        double complex :: phi_next(1:n), phi_temp(1:n)
        integer i
        phi_next = dcmplx(0d0, 0d0)

        ! n = 0
        phi_temp = phi
        phi_next = phi_next + phi_temp

        ! We expect the series converges well enough at 4th order
        do i = 1, 4
            phi_temp = -ix(symmatmul(H, phi_temp, n), n)*dt / (HBAR*dble(i))
            phi_next = phi_next + phi_temp
        end do
    end function

    ! Multiplication of a Symmetric matix and a vector
    ! C = HB
    ! H : Double precision matrix
    ! B : Double complex vector
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
end module

program main
    use extension
    implicit none
    ! Space size (Assuming regular cube), and Time size
    integer,parameter :: n = 40, m = 630
    ! Step amount of space, and time respectively
    double precision,parameter :: dh = 0.8d0, dt = 0.02d0
    ! Lowest bound of space (Represent each coordinates variables as x)
    double precision,parameter :: xl = -dh*(n/2)
    ! Three dimensional wave function
    ! Each index indicates x, y, z respectively
    double complex phi(1:n, 1:n, 1:n)
    ! Hamiltonian
    double precision H(1:n, 1:n)
    ! Other variables
    double precision t1, t2
    phi = dcmplx(0d0, 0d0)

    write (*, *) "Start Calculation..."
    call cpu_time(t1)
    call solve(n, m, dh, dt, xl, H)
    call cpu_time(t2)
    write (*, '(A, F15.10, A)') "Calculation Time : ", t2 - t1, " seconds"
end program