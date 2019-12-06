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
  subroutine initialize(Phi_next, Phi_prev, Pot, N, dh, xmax, gamma)
        integer,intent(in)              :: N
        complex(kind(0d0)),intent(out)  :: Phi_next(0:N, 0:N), Phi_prev(0:N, 0:N)
        double precision,intent(out)    :: Pot(0:N, 0:N)
        double precision,intent(in)     :: dh, xmax, gamma
        double precision                :: x, y, dummy, real_part, imag_part
        integer                         :: i, j
        Phi_next(:, :) = dcmplx(0d0, 0d0)

        if (access("data_input.txt", "") > 0) then
            stop "Input file 'data_input.txt' cannot be found"
        end if
        open(30, file="data_input.txt")
        do j = 0, N
            do i = 0, N
                read (30, *) dummy, dummy, dummy, real_part, imag_part
                Phi_prev(i,j) = dcmplx(real_part, imag_part)
            end do
            read (30, *) 
        end do
        close(30)

        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i

                ! External potential
                !Pot(i, j) = 0.5d0*(x*x+gamma*gamma*y*y)
                if (abs(x) < 10d0 .and. abs(y) < 10d0) then
                    Pot(i, j) = -5d0
                else
                    Pot(i, j) = 0d0
                end if
            end do
        end do
  end subroutine initialize

  subroutine vary_potential(Pot, Pot_TD, N, dh, iter, xmax)
    integer,intent(in)             :: N, iter
    double precision,intent(in)    :: Pot(0:N, 0:N), dh, xmax
    double precision,intent(inout) :: Pot_TD(0:N, 0:N)
    integer                        :: i, j
    double precision               :: R_0 = 2d0
    double precision               :: OMEGA = 2d0*acos(-1d0)/2000
    double precision               :: v_x, v_y
    double precision               :: x_s, y_s, x, y
    double precision               :: sigma = 0.5d0
    integer,parameter              :: mode = 1

    v_x = 2d0*xmax/10000d0
    v_y = v_x

    if (mode == 1) then
        ! Linear Stirring
        do j = 0, N
            y_s = v_y * iter - xmax
            y = -xmax + dh*j
            do i = 0, N
                x_s = v_x * iter - xmax
                x = -xmax + dh*i
        
                Pot_TD(i, j) = Pot(i, j) + 10d0*exp(-0.5d0*((x-x_s)**2d0+(y-y_s)**2d0)/sigma**2d0)
            end do
        end do
    else if (mode == 0) then
        ! Circular Stirring
        do j = 0, N
            y_s = R_0*sin(OMEGA*iter)
            y = -xmax + dh*j
            do i = 0, N
                x_s = R_0*cos(OMEGA*iter)
                x = -xmax + dh*i
                Pot_TD(i, j) = Pot(i, j) + 10d0*exp(-0.5d0*((x-x_s)**2d0+(y-y_s)**2d0)/sigma**2d0)
            end do
        end do
    else
        stop "Invalid mode was given"
    end if
  end subroutine vary_potential
end module