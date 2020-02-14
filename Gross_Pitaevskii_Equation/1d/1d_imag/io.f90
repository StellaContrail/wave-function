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

    ! Output to a file
    subroutine output_real(unit, f, N, dh, xmax)
        integer,intent(in)            :: unit, N
        double precision,intent(in)   :: dh, xmax
        double precision,intent(in)   :: f(0:N)
        double precision x
        integer i
        do i = 0, N
            x = -xmax + dh * i
            write (unit, '(*(F10.5, X))') x, f(i)
        end do
        write (unit, *)
    end subroutine output_real

    ! Output probability current to a file
    subroutine output_current(unit, f, N, dh, xmax)
        integer,intent(in)            :: unit, N
        double precision,intent(in)   :: dh, xmax, f(0:N)
        double precision x
        integer i
        do i = 0, N
            x = -xmax + dh * i
            write (unit, '(*(F10.5, X))') x, f(i)
        end do
        write (unit, *)
    end subroutine

    ! Print string to display
    ! string : content to print
    ! enable : whether to enable print feature
    ! mode   : print mode 'A'lways or 'E'very some iterations
    ! every  : (necessary when mode is 'E')
    ! iter   : (necessary when mode is 'E')
    ! When mode is set to be 'E', print every "every" with present "iter"ation
    subroutine print_ex(string, enable, mode, every, iter)
        character(:),allocatable,intent(in)     :: string
        character(len=1),intent(in)             :: mode
        logical,intent(in)                      :: enable
        integer,optional,intent(in)             :: every, iter
        if (enable) then
            select case(mode)
                case ('a', 'A')
                    print '(A)', string
                case ('e', 'E')
                    if (mod(iter, every) == 0) then
                        print '(A)', string
                    end if
            end select
        end if
    end subroutine print_ex
end module io