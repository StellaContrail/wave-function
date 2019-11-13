module io
    implicit none
contains
    ! Output to a file
    subroutine output(unit, f, N, dh, xmax)
        integer,intent(in)            :: unit, N
        double precision,intent(in)   :: dh, xmax
        complex(kind(0d0)),intent(in) :: f(0:N)
        double precision x
        integer i
        do i = 0, N
            x = -xmax + dh * i
            write (unit, '(*(F10.5, X))') x, real(f(i)), aimag(f(i)), abs(f(i))**2d0
        end do
        write (unit, *)
    end subroutine output
end module io