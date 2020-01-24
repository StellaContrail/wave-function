! Set up wave functions and potential
module setting
    use constants
    use mathf
    implicit none
contains
    ! Initialize wave functions and potential
    subroutine initialize(Pot, mode, Phi)
        complex(kind(0d0)),intent(out),optional  :: Phi(0:N, 0:N)
        double precision,intent(out)             :: Pot(0:N, 0:N)
        integer,intent(in)                       :: mode
        double precision                         :: x, y
        integer                                  :: i, j

        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i
                if (present(Phi)) then
                    ! Initial trial wave function
                    Phi(i, j) = exp(-0.5d0*(x*x+y*y))
                end if

                ! External potential
                select case (mode)
                case (0)
                    ! Harmonic Oscillator Trap
                    Pot(i, j) = 0.5d0*(x*x+gamma*gamma*y*y)
                case (1)
                    ! Harmonic Oscillator Trap and Very Narrow Gaussian-shaped Wall at the center
                    Pot(i, j) = 0.5d0*(x*x+gamma*gamma*y*y) + 100d0*exp(-0.5d0*(x*x+y*y)/(0.5d0*(0.5d0)**2d0))
                case (2)
                    ! Box Trap
                    if (abs(x) < 3d0 .and. abs(y) < 3d0) then
                        Pot(i, j) = -50d0
                    else
                        Pot(i, j) = 0d0
                    end if
                case (3)
                    ! Circle Trap
                    if (sqrt(x**2d0+y**2d0) < 3d0) then
                        Pot(i, j) = -50d0
                    else
                        Pot(i, j) = 0d0
                    end if
                case (4)
                    ! Axially-symmetry Harmonic Oscillator Potential and Split Wall
                    Pot(i, j) = 0.5d0*(x*x*2d0+gamma*gamma*y*y*0.06d0)*0.1d0
                    if (abs(y) < 1d0) then
                        Pot(i, j) = 30d0
                    end if
                case (5)
                    ! Pinning Grid with circulary symmetric trap
                    Pot(i, j) = pinning_potential(x, y)
                case default
                    stop "Invalid mode of external potential"
                end select
            end do
        end do
    end subroutine initialize

    function pinning_potential(x, y) result(V)
        double precision,intent(in) :: x, y
        double precision            :: V, r, rdiff
        ! Radius from the origin point
        r = sqrt(x*x + y*y)
        ! Radius from the point (x0, y0)
        rdiff = sqrt((x-x0)**2d0 + (y-y0)**2d0)

        ! Circularly symmetric trap
        V = 0.5d0*Vmax*(tanh(2*(r-R0))+1d0)
        ! Pinning trap located at (x0, y0)
        V = V + V0*(tanh(delta*(rdiff-alpha))-1d0)
    end function
end module