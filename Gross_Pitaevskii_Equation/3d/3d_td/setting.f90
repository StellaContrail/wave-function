module setting
    implicit none
contains
    ! Initialization
    subroutine initialize(Phi_next, Phi_prev, Pot, N, dh, xmax, gamma_y, gamma_z)
        integer,intent(in)              :: N
        complex(kind(0d0)),intent(out)  :: Phi_next(0:N, 0:N, 0:N), Phi_prev(0:N, 0:N, 0:N)
        double precision,intent(out)    :: Pot(0:N, 0:N, 0:N)
        double precision,intent(in)     :: dh, xmax, gamma_y, gamma_z
        double precision                :: x, y, z, dummy, real_part, imag_part
        integer                         :: i, j, k
        ! Input wave function data file
        character(*),parameter          :: fn_input = "data_input.txt"
        ! sigma : Width of Gaussian's wave packet formed potential
        double precision,parameter      :: sigma = 0.5d0
        ! mode  : Specify type of potential forms
        integer,parameter               :: mode = 0
        Phi_next(:, :, :) = dcmplx(0d0, 0d0)

        open(30, file=fn_input)
        do k = 0, N
            do j = 0, N
                do i = 0, N
                    read (30, *) dummy, dummy, dummy, dummy, real_part, imag_part
                    Phi_prev(i,j,k) = dcmplx(real_part, imag_part)
                end do
                read (30, *) 
            end do
        end do
        close(30)

        do k = 0, N
            z = -xmax + dh*k
            do j = 0, N
                y = -xmax + dh*j
                do i = 0, N
                    x = -xmax + dh*i

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
                        ! Axially-symmetry Harmonic Oscillator Potential
                        Pot(i,j,k) = 0.5d0*(x*x*2d0+gamma_y*gamma_y*y*y*0.06d0+gamma_z*gamma_z*z*z*2d0)
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

    ! Vary potential form depending on time
    subroutine vary_potential(Pot, Pot_TD, N, dh, pi, iter, total_iter, xmax)
        integer,intent(in)             :: N, iter, total_iter
        double precision,intent(in)    :: Pot(0:N, 0:N, 0:N), dh, xmax, pi
        double precision,intent(inout) :: Pot_TD(0:N, 0:N, 0:N)
        integer                        :: i, j, k
        double precision               :: x_s, y_s, z_s, x, y, z
        ! < Circular Stirring Options >
        ! R_0   : Radius of circular stirring
        double precision,parameter     :: R_0 = 2d0
        ! OMEGA : Angular velocity of circular stirring (defined later)
        double precision               :: OMEGA
        ! count : How many times does it stir circulaly
        integer                        :: count = 1
        ! < Linear Stirring Option >
        ! v_x, v_y : Velocity of stirring (defined later)
        double precision               :: v_x, v_y, v_z
        ! < Common Options >
        ! sigma : Stirring potential width
        double precision               :: sigma = 0.5d0
        ! mode  : Specify type of stirring
        integer,parameter              :: mode = 0
        OMEGA = 2d0*pi/(dble(total_iter)/count)
        v_x = 2d0*xmax/total_iter
        v_y = v_x
        v_z = v_x

        select case (mode)
        case (0)
            ! Non Stirring
            Pot_TD(:, :, :) = Pot(:, :, :)
        case (1)
            ! Linear Stirring
            do k = 0, N
                z_s = v_z * iter - xmax
                z = -xmax + dh*k
                do j = 0, N
                    y_s = v_y * iter - xmax
                    y = -xmax + dh*j
                    do i = 0, N
                        x_s = v_x * iter - xmax
                        x = -xmax + dh*i

                        Pot_TD(i, j, k) = Pot(i, j, k) + 10d0*exp(-0.5d0*((x-x_s)**2d0+(y-y_s)**2d0+(z-z_s)**2d0)/sigma**2d0)
                    end do
                end do
            end do
        case (2)
            ! Circular Stirring
            do k = 0, N
                z = -xmax + dh*k
                z_s = z
                do j = 0, N
                    y_s = R_0*sin(OMEGA*iter)
                    y = -xmax + dh*j
                    do i = 0, N
                        x_s = R_0*cos(OMEGA*iter)
                        x = -xmax + dh*i
                        
                        Pot_TD(i, j, k) = Pot(i, j, k) + 10d0*exp(-0.5d0*((x-x_s)**2d0+(y-y_s)**2d0+(z-z_s)**2d0)/sigma**2d0)
                    end do
                end do
            end do
        case default
            stop "Invalid mode of stirring"
        end select
    end subroutine vary_potential
end module