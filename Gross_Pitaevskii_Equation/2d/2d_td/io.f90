module io
    implicit none
contains
    ! Output to a file
    subroutine output(unit, f, N, dh, xmax)
        integer,intent(in)            :: unit, N
        double precision,intent(in)   :: dh, xmax
        complex(kind(0d0)),intent(in) :: f(0:N, 0:N)
        double precision x, y
        integer i, j

        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                write (unit, '(*(F10.5, X))') x, y, abs(f(i,j))**2d0, dble(f(i,j)), aimag(f(i,j))
            end do
            write (unit, *)
        end do
    end subroutine output

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