module io
    implicit none
contains
    subroutine output_real(unit, f, N, dh, xmax)
        integer,intent(in)            :: unit, N
        double precision,intent(in)   :: dh, xmax
        double precision,intent(in) :: f(0:N)
        double precision x
        integer i
        do i = 0, N
            x = dh * i -0.5d0*dh*N
            write (unit, '(*(F10.5, X))') x, abs(f(i))**2d0
        end do
        write (unit, *)
    end subroutine output_real

    subroutine output_potential(unit, f, N, dh, xmax)
        integer,intent(in)            :: unit, N
        double precision,intent(in)   :: dh, xmax
        double precision,intent(in)   :: f(0:N)
        double precision x
        integer i
        do i = 0, N
            x = dh * i -0.5d0*dh*N
            write (unit, '(*(F10.5, X))') x, f(i)
        end do
        write (unit, *)
    end subroutine output_potential
end module io