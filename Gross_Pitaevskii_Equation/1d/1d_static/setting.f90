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
      double precision,intent(out)    :: Phi_next(0:N), Phi_prev(0:N)
      double precision,intent(out)    :: Pot(0:N)
      double precision,intent(in)     :: dh, xmax
      integer                         :: i
      double precision                :: x
      double precision,parameter      :: sigma = 0.8d0
      Phi_next(:) = 0d0
      do i = 0, N
         x = dh*i -0.5d0*dh*N
         ! Potential
         Pot(i) = 0.5d0*x*x
         ! Initial Trial Wave function
         Phi_prev(i) = exp(-0.5d0*x*x/9d0)
      end do
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
            H(i, i-2) = -1008d0 * coe
         end if
         if (i < N-1) then
            H(i, i+2) = -1008d0 * coe
         end if
         if (i > 2) then
            H(i, i-3) = 128d0 * coe
         end if
         if (i < N-2) then
            H(i, i+3) = 128d0 * coe
         end if
         if (i > 3) then
            H(i, i-4) = -9d0 * coe
         end if
         if (i < N-3) then
            H(i, i+4) = -9d0 * coe
         end if
      end do

      ! Potential and Nonlinear part
      do i = 0, N
         H(i, i) = H(i, i) + Pot(i) + kappa*density(i)
      end do
  end subroutine hamiltonian
end module
