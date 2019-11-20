module io
    implicit none
contains
    ! Output to a file
    subroutine output(unit, f, N, dh, xmax)
        integer,intent(in)            :: unit, N
        double precision,intent(in)   :: dh, xmax
        complex(kind(0d0)),intent(in) :: f(0:N)
        double precision x, real_part, imag_part, prob
        integer i
        do i = 0, N
            x = -xmax + dh * i
            real_part = real(f(i))
            imag_part = aimag(f(i))
            prob      = abs(f(i))**2d0
            if (prob < 1d-15) then
                write (unit, '(*(F10.5, X))') x, real_part, imag_part, prob, 0d0
            else
                write (unit, '(*(F10.5, X))') x, real_part, imag_part, prob, datan2(aimag(f(i)), real(f(i)))
            end if
        end do
        write (unit, *)
    end subroutine output
end module io