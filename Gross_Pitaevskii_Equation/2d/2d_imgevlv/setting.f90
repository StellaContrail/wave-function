module setting
    use constants
    implicit none
contains
    ! Initialize wave functions and potential
    subroutine initialize(Phi_next, Phi_prev, Pot, N, dh, xmax, gamma)
        integer,intent(in)              :: N
        complex(kind(0d0)),intent(out)  :: Phi_next(0:N, 0:N), Phi_prev(0:N, 0:N)
        complex(kind(0d0)),intent(out)  :: Pot(0:N, 0:N)
        double precision,intent(in)     :: dh, xmax, gamma
        double precision                :: x, y
        integer                         :: i,j
        ! sigma : Width of Gaussian's wave packet formed potential
        double precision,parameter      :: sigma = 0.5d0
        ! mode  : Specify type of potential forms
        integer,parameter               :: mode = 5
        ! R_0   : Radius of circle or box half width
        double precision,parameter      :: R_0 = 4.5d0
        Phi_next(:, :) = dcmplx(0d0, 0d0)

        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i
                ! Initial trial wave function
                Phi_prev(i, j) = dcmplx(exp(-0.5d0*(x*x+y*y)), 0d0)

                ! External potential
                select case (mode)
                case (0)
                    ! Harmonic Oscillator Trap
                    Pot(i, j) = dcmplx(0.5d0*(x*x+gamma*gamma*y*y), 0d0)
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
                    Pot(i, j) = 200d0*(1d0+tanh(2d0*(sqrt(x*x+y*y)-8.5d0)))
                    if (1.4d0 < x .and. x < 2d0 .and. 1.4d0 < y .and. y < 2d0 .or. .true.) then
                        Pot(i, j) = Pot(i, j) + 2d0*60d0*(1d0+tanh(sqrt((x-1.7d0)**2d0+(y-1.7d0)**2d0)))
                    end if
                case default
                    stop "Invalid mode of external potential"
                end select
            end do
        end do
    end subroutine initialize
end module