! Set up wave functions and potential
module setting
    use constants
    use mathf
    implicit none
contains
    ! Initialize wave functions and potential
    subroutine initialize(Phi, Pot)
        complex(kind(0d0)),intent(out)  :: Phi(0:N, 0:N)
        double precision,intent(out)  :: Pot(0:N, 0:N)
        double precision                :: x, y, dummy, real_part, imag_part
        integer                         :: i, j
        ! sigma : Width of Gaussian's wave packet formed potential
        double precision,parameter      :: sigma = 0.5d0
        ! mode  : Specify type of potential forms
        integer,parameter               :: mode = 5
        ! R_0   : Radius of circle or box half width
        double precision,parameter      :: R_0 = 3d0
        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i
                ! Initial trial wave function
                Phi(i, j) = exp(-0.5d0*(x*x+y*y))

                ! External potential
                select case (mode)
                case (0)
                    ! Harmonic Oscillator Trap
                    Pot(i, j) = 0.5d0*(x*x+gamma*gamma*y*y)
                case (1)
                    ! Harmonic Oscillator Trap and Very Narrow Gaussian-shaped Wall at the center
                    Pot(i, j) = 0.5d0*(x*x+gamma*gamma*y*y) + 100d0*exp(-0.5d0*(x*x+y*y)/(0.5d0*sigma**2d0))
                case (2)
                    ! Box Trap
                    if (abs(x) < R_0 .and. abs(y) < R_0) then
                        Pot(i, j) = -50d0
                    else
                        Pot(i, j) = 0d0
                    end if
                case (3)
                    ! Circle Trap
                    if (sqrt(x**2d0+y**2d0) < R_0) then
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
                    Pot(i, j) = pinning_potential(x, y, 200d0, 120d0, 6.5d0, -1.5d0, 1.5d0, 4d0)
                case default
                    stop "Invalid mode of external potential"
                end select
            end do
        end do
        
        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i
    
                if ((x+1.5d0)**2d0 + (y-1.5d0)**2d0 < 1.5d0**2d0) then
                    Phi(i,j) = Phi(i,j) * exp(-iu*phase((y-1.5d0), (x+1.5d0)))
                end if
            end do 
        end do
    end subroutine initialize

    function pinning_potential(x, y, Vmax, V0, R0, x0, y0, delta) result(V)
        double precision,intent(in) :: x, y, Vmax, V0, R0, x0, y0, delta
        double precision            :: V, r, rdiff
        ! Radius from the origin point
        r = sqrt(x*x + y*y)
        ! Radius from the point (x0, y0)
        rdiff = sqrt((x-x0)**2d0 + (y-y0)**2d0)

        ! Circularly symmetric trap
        V = 0.5d0*Vmax*(tanh(2*(r-R0))+1)
        ! Pinning trap located at (x0, y0)
        V = V + V0*(tanh(delta*rdiff)-1)
    end function
end module