module setting
   implicit none
contains
   ! Initialize the wave functions
   ! Phi_next : Wave function of next space
   ! Phi_prev : Wave function of previous space
   ! Pot      : Potential function
   ! N        : Dimension of space exclusing the first element
   ! dh       : Step distance of space
   ! xmax     : Max x position
   subroutine initialize(Phi_next, Phi_prev, Pot, N, dh, xmax)
      integer,intent(in)              :: N
      complex(kind(0d0)),intent(out)  :: Phi_next(0:N), Phi_prev(0:N)
      double precision,intent(out)    :: Pot(0:N)
      double precision,intent(in)     :: dh, xmax
      integer i
      double precision x, dummy, real_part, imag_part
      Phi_next(:) = dcmplx(0d0, 0d0)

      !open(11, file="data_shifted.txt")
      do i = 0, N
         x = -xmax + dh*i
         !Pot(i) = 0.5d0*x*x
         Pot(i) = 0d0
         if (0d0 < x) then
            Pot(i) = 0.8d0
         end if

         ! Read wave function data from a file
         !read (11, *) dummy, real_part, imag_part, dummy
         !Phi_prev(i) = dcmplx(real_part, imag_part)
         ! Assume the form of the initial wave function
         !Phi_prev(i) = exp(-0.5d0*x*x/(Azero*Azero))+x*exp(-0.5d0*x*x/(Azero*Azero))
         !Phi_prev(i) = exp(-0.5d0*(x+0.5d0)**2d0)
         Phi_prev(i) = exp(-0.5d0*((x+2d0)/(0.5d0))**2d0)*exp(dcmplx(0d0,1d0)*1d0*x)
      end do
      !close(11)
   end subroutine initialize

   ! Construct the hamiltonian
   subroutine hamiltonian(H, Pot, density, N, dh, epsilon, kappa)
      integer,intent(in)            :: N
      double precision,intent(in)   :: dh, Pot(0:N), epsilon, kappa, density(0:N)
      double precision,intent(out)  :: H(0:N, 0:N)
      integer                       :: i
      double precision              :: coe
      coe = -0.5d0 * epsilon * epsilon / (5040d0 * dh * dh) ! Coefficient of the laplacian part
      H(:, :) = 0d0

      ! Laplacian part
      do i = 0, N
         H(i, i)      = -14350d0 * coe
         if (i > 0) then
            H(i, i-1) = 8064d0 * coe
         end if
         if (i < N) then
            H(i, i+1) = 8064d0 * coe
         end if
         if (i > 1) then
         H(i, i-2)  = -1008d0 * coe
      end if
      if (i < N-1) then
         H(i, i+2)  = -1008d0 * coe
      end if
      if (i > 2) then
         H(i, i-3)  = 128d0 * coe
      end if
      if (i < N-2) then
         H(i, i+3)  = 128d0 * coe
      end if
      if (i > 3) then
         H(i, i-4)  = -9d0 * coe
      end if
      if (i < N-3) then
         H(i, i+4)  = -9d0 * coe
      end if
      end do

      ! Potential and Nonlinear part
      do i = 0, N
         H(i, i) = H(i, i) + Pot(i) + kappa*density(i)
      end do
   end subroutine hamiltonian

   subroutine change_potential(Pot_Init, Pot_TD, N, iter) 
      integer,intent(in)             :: N, iter
      double precision,intent(in)    :: Pot_Init(0:N)
      double precision,intent(inout) :: Pot_TD(0:N)
      double precision               :: scale
      integer                        :: i
      
      scale = 1d0 - 0.00002d0 * iter
      if (scale > 0.02d0) then
         do i = 0, N
            Pot_TD(i) = Pot_Init(i) * scale
         end do
      else
         Pot_TD(i) = Pot_TD(i)
      end if
   end subroutine
end module