module setting
    implicit none
contains
    ! Initialize the wave functions and potential
  subroutine initialize(Phi_next, Phi_prev, Pot, N, dh, xmax, gamma_y, gamma_z)
      integer,intent(in)            :: N
      double precision,intent(out)  :: Phi_next(0:N, 0:N, 0:N), Phi_prev(0:N, 0:N, 0:N)
      double precision,intent(out)  :: Pot(0:N, 0:N, 0:N)
      double precision,intent(in)   :: dh, xmax, gamma_y, gamma_z
      double precision              :: x, y, z
      integer                       :: i, j, k
      double precision,parameter    :: sigma = 0.5d0
      integer,parameter             :: mode = 0
      Phi_next(:, :, :) = 0d0

      do k = 0, N
        z = -xmax + dh*k
        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i
                ! Initial trial wave function
                Phi_prev(i,j,k) = exp(-0.5d0*(x*x+y*y+z*z))

                ! External potential
                select case (mode)
                case (0)
                    ! Harmonic Oscillator Trap
                    Pot(i,j,k) = 0.5d0*(x*x+gamma_y*gamma_y*y*y+gamma_z*gamma_z*z*z)
                case (1)
                    ! Harmonic Oscillator Trap and Very Narrow Gaussian-shaped Wall at the center
                    Pot(i,j,k) = 0.5d0*(x*x+gamma_y*gamma_y*y*y+gamma_z*gamma_z*z*z)
                    Pot(i,j,k) = Pot(i,j,k) + 100d0*exp(-0.5d0*(x*x+y*y+z*z)/(0.5d0*sigma**2d0))
                case (2)
                    ! Box Trap
                    if (abs(x) < 7d0 .and. abs(y) < 7d0 .and. abs(z) < 7d0) then
                        Pot(i,j,k) = -5d0
                    else
                        Pot(i,j,k) = 0d0
                    end if
                case (3)
                    ! Sphere Trap
                    if (sqrt(x**2d0+y**2d0+z**2d0) < 7d0) then
                        Pot(i,j,k) = -5d0
                    else
                        Pot(i,j,k) = 0d0
                    end if
                case (4)
                    ! Axially-symmetry Harmonic Oscillator Potential and Split Wall
                    Pot(i,j,k) = 0.5d0*(x*x*2d0+gamma_y*gamma_y*y*y*0.06d0+gamma_z*gamma_z*z*z*2d0)
                    if (abs(y) < 1d0) then
                        Pot(i,j,k) = 30d0
                    end if
                case (5)
                    ! Cylinder Trap
                    if (sqrt(x**2d0+y**2d0) < 5d0 .and. abs(z) < 8d0) then
                        Pot(i,j,k) = -5d0
                    else
                        Pot(i,j,k) = 0d0
                    end if
                case default
                    stop "Invalid mode of external potential"
                end select
            end do
        end do
        end do
  end subroutine initialize
end module