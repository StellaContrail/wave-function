module extensions
    implicit none
    double precision,parameter :: HBAR = 1d0  ! Reduced Plank Constant
    double precision,parameter :: ALPHA = 1d0 ! (HBAR*C)^2/MC^2
  contains
    ! Multiply imaginary unit to an arbitrary complex z
    ! i(a+ib) = -b + ia = -(imag) + i(real)
    double complex function ix(z)
      double complex, intent(in) :: z
      ix = dcmplx(-imag(z), dble(z))
    end function ix

    ! Initialize given wave function array
    ! phi : wave function
    ! n : dimension of phi
    ! xl, xu : lower and upper bound for range of x
    ! x0 : expected value (arbitrary)
    ! a : inverse of standard deviation
    ! dh : step of x
    subroutine initialize(phi, n, xl, x0, a, dh)
      double complex, intent(out) :: phi(0:, 0:)
      double precision, intent(in) :: xl, x0, a, dh
      integer, intent(in) :: n
      integer i
      double precision x
      x = 0d0
      ! Reset all wave function data
      phi = (0d0, 0d0)

      ! Initialize the wave function at t=0 with initialize condition
      do i = 0, n
         x = xl + dh * i
         phi(i, 0) = dcmplx(exp(-0.5d0*a*a*(x-x0)**2d0), 0d0)
      end do

      ! Apply boundary condition at x=xl, x=xu
      phi(0, :) = (0d0, 0d0)
      phi(n, :) = (0d0, 0d0)
    end subroutine initialize

    ! Solve Schroedinger equation
    ! phi : wave function
    ! n : dimension of the wave function
    ! xl, xu : lower and upper bound of range x
    ! dh : step of space
    ! dt : step of time
    subroutine solve_schroedinger(phi, n, m, xl, dh, dt)
      double complex, intent(inout) :: phi(0:, 0:)
        double precision, intent(in) :: xl, dh, dt
        integer, intent(in) :: n, m
        double complex k1(0:n), k2(0:n), k3(0:n), k4(0:n)
        double precision x
        integer i, j
        
        do j = 0, m-1

           ! when i=0, i=n (which are the indices of both ends'), wave functions are always zero
           ! so no need to implementation for the calculation
           ! however, when solving for phi(1), we need information of k(1) which means we are short of k(0)
           ! this causes a problem because when we solve for k(0), we need k(-1) and so on...
           do i = 1, n-1
              x = xl + i * dh
              k1(i) = -ix(H(phi(i-1, j), phi(i, j), phi(i+1,j), x, dh)/HBAR)*dt
           end do
           
           do i = 1, n-1
              x = xl + i * dh
              k2(1) = -ix(H(phi(0, j), phi(1, j)+0.5d0*k1(1), phi(2, j)+0.5d0*k1(2), x, dh)/HBAR)*dt
              k2(n-1) = -ix(H(phi(n-2, j)+0.5d0*k1(n-2), phi(n-1, j)+0.5d0*k1(n-1), phi(n,j), x, dh)/HBAR)*dt
              if (1 < i .and. i < n-1) then
                 k2(i) = -ix(H(phi(i-1, j)+0.5d0*k1(i-1), phi(i, j)+0.5d0*k1(i), phi(i+1,j)+0.5d0*k1(i+1), x, dh)/HBAR)*dt
              end if
           end do
            
           do i = 1, n-1
              x = xl + i * dh
              k3(1) = -ix(H(phi(0, j), phi(1, j)+0.5d0*k2(1), phi(2, j)+0.5d0*k2(2), x, dh)/HBAR)*dt
              k3(n-1) = -ix(H(phi(n-2, j)+0.5d0*k2(n-2), phi(n-1, j)+0.5d0*k2(n-1), phi(n,j), x, dh)/HBAR)*dt
              if (1 < i .and. i < n-1) then
                 k3(i) = -ix(H(phi(i-1, j)+0.5d0*k2(i-1), phi(i, j)+0.5d0*k2(i), phi(i+1,j)+0.5d0*k2(i+1), x, dh)/HBAR)*dt
              end if
           end do
           
           do i = 1, n-1
              x = xl + i * dh
              k4(1) = -ix(H(phi(0, j), phi(1, j)+k3(1), phi(2, j)+k3(2), x, dh)/HBAR)*dt
              k4(n-1) = -ix(H(phi(n-2, j)+k3(n-2), phi(n-1, j)+k3(n-1), phi(n,j), x, dh)/HBAR)*dt
              if (1 < i .and. i < n-1) then
                 k4(i) = -ix(H(phi(i-1, j)+k3(i-1), phi(i, j)+k3(i), phi(i+1,j)+k3(i+1), x, dh)/HBAR)*dt
              end if             
            end do
           
            do i = 1, n-1
               phi(i, j+1) = phi(i, j) + (k1(i) + 2d0*k2(i) + 2d0*k3(i) + k4(i)) / 6d0
            end do
            
        end do
    end subroutine

    ! Define hamiltonian when V = 0.5*x*x
    ! phi_back, phi_center, phi_forward : wave function placed on different positions (back, center, forward)
    ! x : position of center wave function
    ! dh : step of space
   double complex function H(phi_back, phi_center, phi_forward, x, dh)
      double complex,intent(in) :: phi_back, phi_center, phi_forward
      double precision,intent(in) :: dh, x
      double precision re, im

      re = -0.5d0 * ALPHA *( ((dble(phi_forward)-2d0*dble(phi_center)+dble(phi_back))/(dh*dh)) + dble(V(phi_center, x)) )
      im = -0.5d0 * ALPHA *( ((imag(phi_forward)-2d0*imag(phi_center)+imag(phi_back))/(dh*dh)) + imag(V(phi_center, x)) )
      H = dcmplx(re, im)
   end function

   double complex function V(phi, x)
      double complex,intent(in) :: phi
      double precision,intent(in) :: x
      
      V = 0.5d0*x*x*phi   
   end function

    ! Normalization of wave function
    ! phi : wave function
    ! dh : step of space
    ! n : dimension of space
    subroutine normalize(phi, dh, n)
      double complex,intent(inout) :: phi(0:, 0:)
      double precision,intent(in) :: dh
      integer,intent(in) :: n

      phi(:, :) = phi(:, :) / sqrt(summation(phi, dh, n, 0))
    end subroutine normalize

    double precision function summation(phi, dh, n, j)
      double complex, intent(in) :: phi(0:, 0:)
      double precision,intent(in) :: dh
      integer,intent(in) :: n, j
      double precision sum
      integer i

      sum = 0d0
      do i = 0, n
         if (i /= 0 .or. i /= n) then
            sum = sum + abs(phi(i, j))**2d0 * dh
         else
            sum = sum + 0.5d0 * abs(phi(i, j))**2d0 * dh
         end if
      end do
      summation = sum
    end function summation

    subroutine average_x(phi, dh, n, m, xl, avg_x)
      double complex, intent(in) :: phi(0:, 0:)
      double precision, intent(in) :: dh, xl
      double precision, intent(out) :: avg_x(0:)
      integer,intent(in) :: n, m
      double precision sum, x
      integer i, j

      do j = 0, m
         sum = 0d0
         do i = 0, n
            x = xl + dh * i
            if (i /= 0 .or. i /= n) then
               sum = sum + abs(phi(i, j))**2d0 * x * dh
            else
               sum = sum + 0.5d0 * abs(phi(i, j))**2d0 * x * dh
            end if
         end do
         avg_x(j) = sum
      end do
    end subroutine average_x

    subroutine average_p(phi, dh, n, m, xl, avg_p)
      double complex, intent(in) :: phi(0:, 0:)
      double complex sum
      double precision, intent(in) :: dh, xl
      double complex, intent(out) :: avg_p(0:)
      double precision x
      integer, intent(in) :: n, m
      integer i, j

      do j = 0, m
         sum = dcmplx(0d0, 0d0)
         do i = 0, n
            x = xl + dh * i
            if (i == 0) then
               sum = sum - 0.5d0*ix( conjg(phi(0, j))*(phi(1,j)-phi(0,j)) )
            else if (i == n) then
               sum = sum - 0.5d0*ix( conjg(phi(n, j))*(phi(n,j)-phi(n-1,j)) )
            else
               sum = sum - ix( conjg(phi(i, j))*(phi(i+1,j)-phi(i,j)) )
            end if
         end do
         avg_p(j) = sum
      end do
    end subroutine average_p
    
    ! Output wave function
    ! mode : output mode
    ! i, j : index of space and time respectively
    ! xl, xu : lower and upper bound of range x
    ! n, m : dimensions of space and time respectively
    ! dh, dt : step of space and time respectively
    subroutine output(phi, mode, i, j, xl, n, m, dh, dt)
        double complex,intent(in) :: phi(0:, 0:)
        double precision,intent(in) :: xl, dt, dh
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
                   write (fo, *) time, xl + dh*l, dble(phi(l, k)), imag(phi(l, k)), abs(phi(l, k))**2d0
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
    integer,parameter :: n = 100, m = 2500
    double precision :: xl = -10d0, x0 = 0d0, a = 1d0
    double precision :: avg_x(0:m)
    double complex :: phi(0:n, 0:m) = (0d0, 0d0), avg_p(0:m) = (0d0, 0d0)
    double precision,parameter :: dh = 0.2d0, dt = 0.005d0

    call initialize(phi, n, xl, x0, a, dh)
    call solve_schroedinger(phi, n, m, xl, dh, dt)
    call normalize(phi, dh, n)
    call output(phi, "A", 0, 0, xl, n, m, dh, dt)
    write (*, *)
    write (*, *) "Probability"
    write (*, '(x, a, 2x, f10.7)') "j=0    : ", sqrt(summation(phi, dh, n, 0))
    write (*, '(x, a, 2x, f10.7)') "j=100  : ", sqrt(summation(phi, dh, n, 100))
    write (*, '(x, a, 2x, f10.7)') "j=1000 : ", sqrt(summation(phi, dh, n, 1000))
    write (*, '(x, a, 2x, f10.7)') "j=2000 : ", sqrt(summation(phi, dh, n, 2000))
    write (*, *)

    call average_x(phi, dh, n, m, xl, avg_x)
    write (*, *) "Expected values of x"
    write (*, '(x, a, 2x, f10.7)') "j=0    : ", avg_x(0)
    write (*, '(x, a, 2x, f10.7)') "j=100  : ", avg_x(100)
    write (*, '(x, a, 2x, f10.7)') "j=1000 : ", avg_x(1000)
    write (*, '(x, a, 2x, f10.7)') "j=2000 : ", avg_x(2000)
    write (*, *)

    call average_p(phi, dh, n, m, xl, avg_p)
    write (*, *) "Expected value of p"
    write (*, *) "REAL IMAG"
    write (*, '(x, a, 2x, f10.7,2x,f12.7)') "j=0    : ", dble(avg_p(0)), imag(avg_p(0))
    write (*, '(x, a, 2x, f10.7,2x,f12.7)') "j=100  : ", dble(avg_p(0)), imag(avg_p(0))
    write (*, '(x, a, 2x, f10.7,2x,f12.7)') "j=1000 : ", dble(avg_p(0)), imag(avg_p(0))
    write (*, '(x, a, 2x, f10.7,2x,f12.7)') "j=2000 : ", dble(avg_p(0)), imag(avg_p(0))
    write (*, *)
end program
