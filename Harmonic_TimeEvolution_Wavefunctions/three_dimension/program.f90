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
    function ix_3d(A_matrix, n)
      integer,intent(in) :: n
      double complex,intent(in) :: A_matrix(1:n, 1:n, 1:n)
      double complex ix(1:n, 1:n, 1:n)
      integer i, j, k
      do k = 1, n
         do j = 1, n
            do i = 1, n
               ix(i,j,k) = dcmplx(-aimag(A_matrix(i,j,k)), dble(A_matrix(i,j,k)))
            end do
         end do
      end do
    end function ix_3d
    function ix_scaler(z)
        double complex,intent(in) :: z
        double complex ix_scaler
        ix_scaler = dcmplx(-aimag(z), dble(z))
    end function

    ! Solve all space and time steps of harmonic oscillator
    subroutine solve(n, m, dh, dt, xl)
        integer,intent(in) :: n, m
        double precision,intent(in) :: dh, dt, xl
        double complex :: phi_next(1:n, 1:n, 1:n), phi_old(1:n, 1:n, 1:n)
        double precision t1, t2, speed, x, y, z
        integer i, j, k, waste_loop_var, time_loop_var
        phi_old = dcmplx(0d0, 0d0)
        open(10, file="data.txt")

        ! Construction of Hamiltonian
        ! Temporary disabled: I dunno the proper Hamiltonian for this case...
        !call construct_hamiltonian(H, n, xl, dh)

        ! Initialize the wave function
        call initialize(phi_old, n, dh, xl)
        ! Normalize the initial wave function
        call normalize(phi_old, n, dh)
        
        ! Calculate future wave function of all time steps
        do time_loop_var = 1, m/10
           ! Plot data in every 10 calculations
           do waste_loop_var = 1, 10
              phi_next = (0d0, 0d0)
              ! Calculate the future wave function of next time step
              phi_next = calc_next(phi_old, n, dt, dh)
           end do
           
           ! Plot data into a text file
           do k = 1, n
              z = xl + dh*k
              do j = 1, n
                 y = xl + dh*j
                 do i = 1, n
                    x = xl + dh*i
                    write (10, '(4F15.10)', advance='no') x, y, z, abs(phi_next(i,j,k))
                    write (10, *)
                 end do
              end do
           end do
           write (*, *)
        end do
    
        close(10)
    end subroutine

    function apply_H(phi_old, n, dh) result(phi_next)
      integer,intent(in) :: n
      double precision,intent(in) :: dh
      double complex,intent(in) :: phi_old(1:n, 1:n, 1:n)
      double complex phi_next(1:n, 1:n, 1:n)
      integer i, j, k
      
      do k = 1, n
         do j = 1, n
            do i = 1, n
               phi_next(i,j,k) = -6d0*phi_old(i,j,k)
               if (i /= 1) then
                  phi_next(i,j,k) = phi_next(i,j,k) + phi_old(i-1,j,k)
               end if
               if (i /= n) then
                  phi_next(i,j,k) = phi_next(i,j,k) + phi_old(i+1,j,k)
               end if
               if (j /= 1) then
                  phi_next(i,j,k) = phi_next(i,j,k) + phi_old(i,j-1,k)
               end if
               if (j /= n) then
                  phi_next(i,j,k) = phi_next(i,j,k) + phi_old(i,j+1,k)
               end if
               if (k /= 1) then
                  phi_next(i,j,k) = phi_next(i,j,k) + phi_old(i,j,k-1)
               end if
               if (k /= n) then
                  phi_next(i,j,k) = phi_next(i,j,k) + phi_old(i,j,k+1)
               end if
               phi_next(i,j,k) = phi_next(i,j,k) / (dh*dh)

               phi_next(i,j,k) = -(HBAR*HBAR/(2d0*MASS))*phi_next(i,j,k)
               ! I dunno if this is correct, maybe potential can only affect the diagonal elements of the wave function
               phi_next(i,j,k) = phi_next(i,j,k) + V(dh*i,dh*j,dh*k)*phi_next(i,j,k)
            end do
         end do
      end do
    end function apply_H

    ! Construct the Hamiltonian matrix with 5-point stencil
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
    double precision function V(x,y,z)
        double precision,parameter :: BETA = MASS * OMEGA * OMEGA
        double precision,optional,intent(in) :: x,y,z
        double precision r2
        r2 = x*x+y*y+z*z
        V = 0.5d0*r2*BETA
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

    
    function calc_next(phi, n, dh, dt) result(phi_next)
        integer,intent(in) :: n
        double complex,intent(in) :: phi(1:n, 1:n, 1:n)
        double precision,intent(in) :: dt, dh
        double complex :: phi_next(1:n, 1:n, 1:n), phi_temp(1:n, 1:n, 1:n)
        integer i
        phi_next = dcmplx(0d0, 0d0)

        ! n = 0
        phi_temp = phi
        phi_next = phi_next + phi_temp

        ! We expect the series converges well enough at 4th order
        do i = 1, 4
            phi_temp = -ix_3d(apply_H(phi_temp,n,dh),n)*dt / (HBAR*dble(i)))
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
    call solve(n, m, dh, dt, xl)
    call cpu_time(t2)
    write (*, '(A, F15.10, A)') "Calculation Time : ", t2 - t1, " seconds"
end program
