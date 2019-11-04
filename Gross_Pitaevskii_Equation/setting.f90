! Functions/Subroutines which set up configurations for a physical system
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
        complex(kind(0d0)),intent(out)  :: Phi_next(0:N), Phi_prev(0:N)
        double precision,intent(out)    :: Pot(0:N)
        double precision,intent(in)     :: dh, xmax
        integer i
        double precision x
        Phi_next(:) = dcmplx(0d0, 0d0)
        do i = 0, N
            x = -xmax + dh*i
            Pot(i) = 0.5d0*x*x
            if (i == 0 .or. i == N) then
                ! Boundary Condition (Fixed)
                Phi_prev(i) = dcmplx(0d0, 0d0)
            else
                ! Suppose initial wave function is set to be harmonic oscillator
                Phi_prev(i) = dcmplx(exp(-Pot(i)), 0d0)
            end if
        end do
    end subroutine
end module