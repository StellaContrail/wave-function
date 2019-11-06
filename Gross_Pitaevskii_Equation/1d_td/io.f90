module io
    implicit none
contains
    ! Output to a file
  subroutine output(unit, f, N, dh, xmax)
    integer,intent(in)            :: unit, N
    double precision,intent(in)   :: dh, xmax
    complex(kind(0d0)),intent(in) :: f(1:N)
    double precision x
    integer i
    do i = 1, N
       x = -xmax + dh * i
<<<<<<< HEAD
       write (unit, '(*(F15.5, X))') x, real(f(i)), aimag(f(i)), abs(f(i))**2d0
=======
       write (unit, '(F15.5, X, F15.5)') x, abs(f(i))**2d0
>>>>>>> 3096838a4ae31bf9300dfcf01b45c983364a4c68
    end do
    write (unit, *)
  end subroutine output
end module io