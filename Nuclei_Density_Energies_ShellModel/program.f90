! 4.1 Potential of nucleon

module extension
   implicit none
   double precision,parameter :: ALPHA = 939d0/197.323d0**2d0 ! mc^2/(hbar*c)^2
   double precision,parameter :: a = 0.67
   double precision,parameter :: r_zero = 1.27
   double precision,parameter :: COULOMB_CONSTANT = 197.326/137 ! e^2 in Gaussian unit
   double precision,parameter :: PI = acos(-1d0)
contains
   double precision function get_pos(i, dr)
      integer,intent(in) :: i
      double precision,intent(in) :: dr
      get_pos = dr*i - 0.5d0*dr
      !get_pos = dr*i
   end function

   ! solve schroedinger equation
   subroutine solve_schrodinger(n, dr, isNeutron, A_, N_, Z_, l, j, values, H)
      integer,intent(in) :: n, N_, A_, Z_, l
      logical,intent(in) :: isNeutron
      double precision,intent(in) :: dr, j
      double precision,intent(out) :: values(1:n)
      double precision,intent(inout) :: H(1:n, 1:n)
      integer num
    
      ! Construct Hamiltonian
      call construct_hamiltonian(H, n, dr, l, j, isNeutron, N_, A_, Z_)

      ! Solve eigen equation (Schroedinger equation)
      call solve_eigen(H, n, values)
      values = values / ALPHA ! values = (Energy*ALPHA) / ALPHA

      ! num : number of bounding states of the given J
      num = count(values < 0)
      if (num > 0) then
         ! Normalize wave functions
         call normalize(H, n, dr, num)
      end if
   end subroutine solve_schrodinger

   ! Construct Hamiltonian   
   subroutine construct_hamiltonian(H, n, dr, l, j, isNeutron, N_, A_, Z_)
      integer,intent(in) :: n, N_, A_, Z_, l
      logical,intent(in) :: isNeutron
      double precision,intent(in) :: dr, j
      double precision,intent(out) :: H(1:n, 1:n)
      integer i
      double precision r
    
      H = 0d0
      do i = 1, n
         H(i, i) = -30d0
         if (1 < i) then
            H(i, i-1) = 16d0
         end if
         if (2 < i) then
            H(i, i-2) = -1d0
         end if
         if (i < n) then
            H(i, i+1) = 16d0
         end if
         if (i < n-1) then
            H(i, i+2) = -1d0
         end if
      end do
      H(1, 1) = -46d0
      H(1, 2) =  17d0
      H(2, 1) =  17d0
      H = -0.5d0*H/(12d0*dr*dr)
      do i = 1, n
         r = get_pos(i, dr)
         H(i, i) = H(i, i) + ALPHA*V_prime(r, l, j, isNeutron, N_, A_, Z_, dr) + (0.5d0*l*(l+1))/(r*r)
      end do
   end subroutine
   
   ! calculate potential
   double precision function V_prime(r, l, j, isNeutron, N_, A_, Z, dr)
      double precision,intent(in) :: r, j, dr
      integer,intent(in) :: N_, A_, Z, l
      logical,intent(in) :: isNeutron
      double precision ls, V, W, R_, beta
      integer sign_
      if (isNeutron) then
         sign_ = -1
      else
         sign_ = 1
      end if
      ls = 0.5d0*(j*(j+1d0)-l*(l+1)-0.75d0)
      R_ = r_zero*A_**(1d0/3d0)
      beta = (r-R_)/a
    
      ! Woods-Saxon Potential
      V = -(51d0+sign_*33d0*dble(N_-Z)/A_)*f(r, R_)
      
      ! Coulomb Potential
      if (.not. isNeutron) then
          if (r < R_) then
             V = V + 0.5d0*((Z-1)*COULOMB_CONSTANT*(3d0-(r/R_)**2d0))/R_
         else
             V = V + 0.5d0*((Z-1)*COULOMB_CONSTANT)/r
         end if
      end if
      
      ! Spin-Orbital Interaction Potential
      W = (22d0+sign_*14d0*dble(N_-Z)/A_)*r_zero*r_zero/r
      W = W * (f(r+0.5d0*dr, R_)-f(r-0.5d0*dr,R_))/dr
      V_prime = V + ls*W
    end function V_prime
    ! woods-saxon potential
    double precision function f(r, R_)
      double precision,intent(in) :: r, R_
      f = (1d0/(1d0+exp((r-R_)/a)))*(1d0/(1d0+exp((-r-R_)/a)))
      !f = (1d0/(1d0+exp((r-R_)/a)))
    end function f
   
   ! solve eigen equation of schroedinger equation
   subroutine solve_eigen(H, n, values)
      integer,intent(in) :: n
      double precision,intent(inout) :: H(1:n, 1:n), values(1:n)
      double precision, allocatable :: dummy(:)
      integer lwork, info
    
      allocate(dummy(1))
      call dsyev("V", "U", n, H, n, values, dummy, -1, info)
      lwork = int(dummy(1))
      deallocate(dummy)
      allocate(dummy(lwork))
      call dsyev("V", "U", n, H, n, values, dummy, lwork, info)
      if (info < 0) then
         stop "There was an error when executing dsyev()"
      end if
      deallocate(dummy)
   end subroutine solve_eigen

   ! normalize wave functions
   subroutine normalize(phi, n, dr, output_max)
      integer,intent(in) :: n, output_max
      double precision,intent(inout) :: phi(1:n, 1:n)
      double precision,intent(in) :: dr
      double precision sum
      integer i, j
      do j = 1, output_max
         sum = 0d0
         do i = 1, n
            if (i == 1 .or. i == n) then
               sum = sum + 0.5d0*phi(i,j)*phi(i,j)*dr
            else
               sum = sum + phi(i,j)*phi(i,j)*dr
            end if
         end do
         phi(:,j) = phi(:,j) / sqrt(sum)
      end do
   end subroutine normalize

   ! solve density of nucleon
   subroutine solve_density(u, dr, n, count, j, rho, num_nucleon)
      double precision,intent(in) :: u(1:n, 1:count), j(1:), dr
      double precision,intent(out) :: rho(1:)
      integer,intent(in) :: count, n, num_nucleon
      integer i, k, num
      double precision r
      double precision sum
      
      num = num_nucleon
      rho = 0d0

      do k = 1, count
         if (num - int(2d0*j(k)+1d0) > 0) then
            num = num - int(2d0*j(k)+1d0)
            do i = 1, n
               r = get_pos(i, dr)
               rho(i) = rho(i) + (2d0*j(k)+1d0)*u(i, k)*u(i, k)
            end do
         else if (num > 0) then
            do i = 1, n
               r = get_pos(i, dr)
               rho(i) = rho(i) + num*u(i, k)*u(i, k)
            end do
            exit
         end if
      end do
      do i = 1, n
         r = get_pos(i, dr)
         rho(i) = rho(i) / (4d0*PI*r*r)
      end do

      ! check number of nucleon
      sum = 0d0
      do i = 1, n
         r = get_pos(i, dr)
         if (i == 1 .or. i == n) then
            sum = sum + 0.5d0*rho(i)*(4d0*PI*r*r)*dr
         else
            sum = sum + rho(i)*(4d0*PI*r*r)*dr
         end if
      end do
      if (abs(sum - num_nucleon) > 0.005d0) then
         write (*, *) "[WARNING] Calculated number of nucleon does not match the input data. Result may not be correct"
      end if
      write (*, '(a, f10.5)') "# of nuclei : ", sum
      write (*, '(a, I3)')    "   Expected : ", num_nucleon
    end subroutine solve_density

    ! output wave function to a file
    subroutine output_wavefunc_to_file(H, n, dr, k_max)
      integer,intent(in) :: n, k_max
      double precision,intent(in) :: H(1:, 1:), dr
      integer k, i
      open(10, file="wave.txt")
      do k = 1, k_max
         do i = 1, n
            write (10, *) get_pos(i, dr), H(i,k)
         end do
         write (10, *)
      end do
      close(10)
      write (*, *)
      write (*, '(2a)') "Wave function data -> ", "wave.txt"
      write (*, *) "{FORMAT} ", "Position[fm]  Amplitude"
   end subroutine

   ! output density to a file
   subroutine output_density_to_file(rho, n, dr)
      integer,intent(in) :: n
      double precision,intent(in) :: rho(1:), dr
      integer i
      open(10, file="density.txt")
      do i = 1, n
         write (10, *) get_pos(i, dr), rho(i)
      end do
      close(10)
      write (*, *)
      write (*, '(2a)') "Nucleon density data -> ", "density.txt"
      write (*, *) "{FORMAT} ", "Position[fm]  Density[fm^-3]"
   end subroutine

   ! sort eigenvalues in ascending order
   subroutine sort_data(phi, phi_sorted, n, energies, energies_sorted, count, k, l, j)
      integer,intent(in) :: count, n
      double precision,intent(in) :: phi(1:, 0:, 1:), energies(1:, 0:)! Row:N, Column:(l-0.5, l+0.5)*l_max
      double precision,intent(out) :: j(1:count), energies_sorted(1:count)
      double precision,intent(inout) :: phi_sorted(1:n, 1:count)
      integer,intent(out) :: k(1:count), l(1:count)
      double precision lowest_energy
      integer i, position(1:2)

      lowest_energy = -100d0
      do i = 1, count
         ! minloc() returns position of minimum value satisfying MASK
         ! note: minloc() ignores lowest bound of array, always regards arrays to start at 1
         position = minloc(energies, energies > lowest_energy) ! 1st:k, 2nd:indexOfJ
         position(2) = position(2) - 1
         lowest_energy = minval(energies, energies > lowest_energy)
         k(i) = position(1)
         l(i) = position(2)/2
         j(i) = position(2)-position(2)/2-0.5d0
         energies_sorted(i) = lowest_energy
         phi_sorted(:, i) = phi(k(i), position(2), :)
      end do
   end subroutine

   ! make plt file for plotting magic numbers
   subroutine output_to_plt(energies, k, l, j, count)
      integer,intent(in) :: count, k(1:), l(1:)
      double precision,intent(in) :: j(1:), energies(1:count)! Row:N, Column:(l-0.5, l+0.5)*l_max
      integer,parameter :: magic_numbers(1:7) = (/2, 8, 20, 28, 50, 82, 126/)
      integer i, total
      character(len=2) l_name
      integer,parameter :: BAR_WIDTH = 750, ARROW_WIDTH = 170, SPACE = 100
      integer :: WIDTH_SUM = BAR_WIDTH + ARROW_WIDTH
      total = 0
      open(11, file="plot.plt")

      write (11, *) "set ylabel 'E [MeV]"
      write (11, *) "unset xtics"
      write (11, *) "unset key"
      write (11, *) "y = 0"
      write (11, *) "plot y"
      write (11, *) "unset arrow"
      write (11, *) "unset label"
      write (11, *) "set yrange[-41:1]"
      write (11, *) "set xrange[-350:", WIDTH_SUM+SPACE+350, "]"
      
      do i = 1, count
         total = total + int(2d0*j(i)+1d0)
         ! energy levels drawing
         write (11, *) "set arrow ", i, " from 0,", energies(i), " to ", BAR_WIDTH, ",", energies(i), " nohead"
          ! magic numbers drawing
         if (any(total == magic_numbers)) then
            write (11, *) "set arrow ",i+count," from ", WIDTH_SUM+SPACE, ",", energies(i), " to ", BAR_WIDTH+SPACE,",", energies(i)
            write (11, *) "set label ", i, '"', total,'"', "at ", WIDTH_SUM+SPACE, ",", energies(i), " center"
         end if
         ! N:L:J drawing
         if (l(i) == 0) then
            l_name = "s"
         else if (l(i) == 1) then
            l_name = "p"
         else if (l(i) == 2) then
            l_name = "d"
         else if (l(i) == 3) then
            l_name = "f"
         else if (l(i) == 4) then
            l_name = "g"
         else if (l(i) == 5) then
            l_name = "h"
         else if (l(i) == 6) then
            l_name = "i"
         else if (l(i) == 7) then
            l_name = "j"
         else if (l(i) == 8) then
            l_name = "k"
         else if (l(i) == 9) then
            l_name = "l"
         else if (l(i) == 10) then
            l_name = "m"
         end if
         write (11, '(a, i3, a, i1, a, i2)', advance="no") "set label ", i+count, '"', k(i), l_name, int(2d0*j(i))
         write (11, '(3a, f15.8, a)') "/2", '"', "at -80,", energies(i), " right"
      end do
      write (11, *) "replot"
      close(11)
      
      write (*, *)
      write (*, '(2a)') "Magic number data -> ", "plot.plt"
      write (*, *) "This is for gnuplot use. Type load ", '"plot.plt"', " in gnuplot and hit enter."
   end subroutine
end module extension

program main
   use extension
   implicit none
   integer,parameter :: n = 500, l_max = 10
   double precision,parameter :: R = 20d0
   double precision :: energies(1:4, 0:l_max*2)! Row:N, Column:(l-0.5, l+0.5)*l_max
   double precision dr, j, values(1:n), u(1:4, 0:l_max*2, 1:n), H(1:n, 1:n)
   integer A_, Z_, l, k, spin, count, i
   integer,allocatable :: k_sorted(:), l_sorted(:)
   double precision,allocatable :: j_sorted(:), energies_sorted(:), u_sorted(:, :), rho(:)
   logical isNeutron
   character(len=2) yesorno 
   dr = R / n
   count = 0

   write (*, '(a)', advance="no") "Input A:"
   read (*, *) A_
   write (*, '(a)', advance="no") "Input Z:"
   read (*, *) Z_
   write (*, '(a)', advance="no") "Calculate neutron?(n->proton) [y/n]:"
   read (*, *) yesorno
   write (*, *)
   isNeutron = yesorno == "y" .or. yesorno == "yes"

   write (*, '(a, f10.5)') "Radius of the nucleus [fm]:", r_zero*A_**(1d0/3d0) 
   ! Angular orbital momentum qunatum number l=[0,l_max(user-defined)]
   do l = 0, l_max
      ! s = spin * 0.5d0 = -0.5d0 or 0.5d0
      do spin = -1, 1, 2
         ! Total momentum : J = L + S
         j = l + spin * 0.5d0
         ! check if j doesnt surpass below zero
         if (j < 0) then
            cycle
         end if
         ! solve schroedinger under the given condition
         call solve_schrodinger(n, dr, isNeutron, A_, A_-Z_, Z_, l, j, values, H)

         ! find states where neucleon is trapped inside the potential
         do k = 1, 4
            ! end this loop when there is no more bound states
            if (values(k) > 0) then
               exit
            end if
            count = count + 1
            energies(k, 2*l+floor(spin/2d0)+1) = values(k)
            u(k, 2*l+floor(spin/2d0)+1, :) = H(:, k)
         end do
      end do
   end do

   allocate(k_sorted(1:count), l_sorted(1:count), j_sorted(1:count), energies_sorted(1:count), u_sorted(1:n, 1:count), rho(1:n))
  
   call sort_data(u, u_sorted, n, energies, energies_sorted, count, k_sorted, l_sorted, j_sorted)
   write (*, *)
   write (*, '(a)') "Energy levels in ascending order :"
   do i = 1, count
      write (*, '(x, 2(a,i1),a,f4.1,a,f8.3)') "N=",k_sorted(i)," L=",l_sorted(i)," J=",j_sorted(i)," : ",energies_sorted(i)
   end do

   write (*, *)
   write (*, '(a)', advance="no") "OUTPUTs :"
   call output_to_plt(energies_sorted, k_sorted, l_sorted, j_sorted, count)
   call output_wavefunc_to_file(u_sorted, n, dr, count)
   if (isNeutron) then
      call solve_density(u_sorted, dr, n, count, j_sorted, rho, A_-Z_)
   else
      call solve_density(u_sorted, dr, n, count, j_sorted, rho, Z_)
   end if
   call output_density_to_file(rho, n, dr)
   deallocate(k_sorted, l_sorted, j_sorted, energies_sorted, u_sorted, rho)
end program main
