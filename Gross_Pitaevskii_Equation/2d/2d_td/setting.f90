module setting
    use constants
    implicit none
contains
    ! Initialize wave functions and potential
    subroutine initialize(Phi, Pot)
        complex(kind(0d0)),intent(out)  :: Phi(0:N, 0:N)
        double precision,intent(out)    :: Pot(0:N, 0:N)
        double precision                :: x, y, dummy, real_part, imag_part
        integer                         :: i, j
        ! Input wave function data file
        character(*),parameter          :: fn_input = "data_input.txt"
        ! sigma : Width of Gaussian's wave packet formed potential
        double precision,parameter      :: sigma = 0.5d0
        ! mode  : Specify type of potential forms
        integer,parameter               :: mode = 5
        ! R_0   : Radius of circle or box half width
        double precision,parameter      :: R_0 = 3d0

        if (.true.) then
            if (access(fn_input, "") > 0) then
                print *, "Input file '", fn_input, "' cannot be found"
                stop
            end if

            open(15, file=fn_input)
            do j = 0, N
                do i = 0, N
                    read (15, *) dummy, dummy, dummy, real_part, imag_part
                    Phi(i,j) = dcmplx(real_part, imag_part)
                end do
                read (15, *) 
            end do
            close(15)
        else 
            do j = 0, N
                y = -xmax + dh * j
                do i = 0, N
                    x = -xmax + dh * i

                    Phi(i,j) = dcmplx(exp(-0.5*(x*x+gamma*gamma*y*y)), 0d0)
                end do
            end do
        end if

        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i

                ! External potential
                select case (mode)
                case (0)
                    ! Harmonic Oscillator Trap
                    Pot(i, j) = 0.5d0*(x*x+gamma*gamma*y*y)
                case (1)
                    ! Harmonic Oscillator Trap and Very Narrow Gaussian-shaped Wall at the center
                    Pot(i, j) = 0.5d0*(x*x+gamma*gamma*y*y) + 100d0*exp(-0.5d0*(x*x+y*y)/(sigma*sigma))
                case (2)
                    ! Box Trap
                    if (abs(x) < R_0 .and. abs(y) < R_0) then
                        Pot(i, j) = -5d0
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
                    ! Axially-symmetry Harmonic Oscillator Potential
                    Pot(i, j) = 0.5d0*(x*x*2d0+gamma*gamma*y*y*0.06d0)*0.1d0
                case (5)
                    ! Pinning Grid with circulary symmetric trap
                    Pot(i, j) = pinning_potential(x, y, 200d0, 120d0, 6.5d0, -1.5d0, 1.5d0, 4d0)
                case default
                    stop "Invalid mode of external potential"
                end select
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
        V = 0.5d0*Vmax*(tanh(2d0*(r-R0))+1d0)
        ! Pinning trap located at (x0, y0)
        V = V + V0*(tanh(delta*rdiff)-1d0)
    end function

  ! Vary potential form depending on time
  subroutine vary_potential(Pot, Pot_TD, iter, iter_max)
    integer,intent(in)             :: iter, iter_max
    double precision,intent(in)    :: Pot(0:N, 0:N)
    double precision,intent(inout) :: Pot_TD(0:N, 0:N)
    integer                        :: i, j
    double precision               :: x_s, y_s, x, y, t
    double precision               :: V_0 = 100d0, V
    ! < Circular Stirring Options >
    ! R_0   : Radius of circular stirring
    double precision,parameter     :: R_0 = 1.8d0
    double precision,parameter     :: stirOMEGA = 0d0
    ! < Linear Stirring Option >
    ! v_x, v_y : Velocity of stirring (defined later)
    double precision               :: v_x, v_y
    ! gradually increase/descrease the intensity of the circularly stirring potential to avoid transient effects
    double precision               :: fade_in_start_clock  = 0d0
    double precision               :: fade_in_end_clock    = 3d0
    double precision               :: fade_out_start_clock = 6d0
    double precision               :: fade_out_end_clock   = 9d0
    ! < Common Options >
    ! sigma : Stirring potential width
    double precision               :: sigma = 0.35d0
    ! mode  : Specify type of stirring
    integer,parameter              :: mode = 0
    v_x = 2d0*xmax/iter_max
    v_y = v_x

    t = dt * iter
    if (t < fade_in_start_clock) then
        V = 0d0
    end if
    if (fade_in_start_clock < t .and. t < fade_in_end_clock) then
        V = V_0 * (t - fade_in_start_clock)/(fade_in_end_clock - fade_in_start_clock)
        if (V > V_0) then
            V = V_0
        end if
    end if
    if (fade_in_end_clock < t .and. t < fade_out_start_clock) then
        V = V_0
    end if
    if (fade_out_start_clock < t .and. t < fade_out_end_clock) then
        V = V_0 - V_0 * (t - fade_out_start_clock)/(fade_out_end_clock - fade_out_start_clock)
        if (V < 0d0) then
            V = 0d0
        end if
    end if
    if (fade_out_end_clock < t) then
        V = 0d0
    end if

    select case (mode)
    case (0)
        ! Non Stirring
        Pot_TD(:, :) = Pot(:, :)
    case (1)
        ! Linear Stirring
        do j = 0, N
            y_s = v_y * iter - xmax
            y = -xmax + dh*j
            do i = 0, N
                x_s = v_x * iter - xmax
                x = -xmax + dh*i
                Pot_TD(i, j) = Pot(i, j) + V*exp(-0.5d0*((x-x_s)**2d0+(y-y_s)**2d0)/sigma**2d0)
            end do
        end do
    case (2)
        ! Circular Stirring
        do j = 0, N
            y_s = R_0*sin(stirOMEGA*iter)
            y = -xmax + dh*j
            do i = 0, N
                x_s = R_0*cos(stirOMEGA*iter)
                x = -xmax + dh*i

                Pot_TD(i, j) = Pot(i, j) + V*exp(-0.5d0*((x-x_s)**2d0+(y-y_s)**2d0)/sigma**2d0)
            end do
        end do
    case default
        stop "Invalid mode of stirring"
    end select
  end subroutine vary_potential
end module