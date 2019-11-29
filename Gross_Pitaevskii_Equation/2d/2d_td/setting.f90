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
                Pot(i, j) = 0.5d0*(x*x+gamma*gamma*y*y)
            end do
        end do
  end subroutine initialize
end module