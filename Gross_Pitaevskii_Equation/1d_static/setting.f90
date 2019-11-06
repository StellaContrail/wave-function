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
  subroutine initialize(Phi_next, Phi_prev, Pot, N, dh, xmax, Azero)
    integer,intent(in)              :: N
    complex(kind(0d0)),intent(out)  :: Phi_next(1:N), Phi_prev(1:N)
    double precision,intent(out)    :: Pot(1:N)
    double precision,intent(in)     :: dh, xmax, Azero
    integer i
    double precision x
    Phi_next(:) = dcmplx(0d0, 0d0)
    do i = 1, N
         x = -xmax + dh*i
         Pot(i) = 0.5d0*x*x
         ! Assume the form of the initial wave function
         Phi_prev(i) = exp(-0.5d0*x*x/(Azero**2d0))
    end do
  end subroutine initialize

  ! Construct the hamiltonian
  subroutine hamiltonian(H, Pot, N, dh, epsilon, kappa, Phi)
    integer,intent(in)            :: N
    double precision,intent(in)   :: dh, Pot(1:N), epsilon, kappa
    double precision,intent(out)  :: H(1:N, 1:N)
    complex(kind(0d0)),intent(in) :: Phi(1:N)
    integer                       :: i
    double precision              :: coe
    coe = -0.5d0 * epsilon * epsilon / (dh * dh) ! Coefficient of the laplacian part
    H(:, :) = 0d0

    ! Laplacian part
    do i = 1, N
       H(i, i) = -2d0 * coe
       if (i > 1) then
          H(i, i-1) = 1d0 * coe
       end if
       if (i < N) then
          H(i, i+1) = 1d0 * coe
       end if
    end do

    ! Potential and Nonlinear part
    do i = 1, N
       H(i, i) = H(i, i) + Pot(i) + kappa*abs(Phi(i))**2d0
    end do
  end subroutine hamiltonian
end module