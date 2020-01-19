module io
    use constants
    implicit none

    interface output
        module procedure output_complex, output_real, output_flux, output_widths
    end interface
contains
    ! Save double precision complex wave function
    subroutine output_complex(unit, f)
        integer,intent(in)            :: unit
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
    end subroutine output_complex

    ! Save potential forms
    subroutine output_real(unit, f)
        integer,intent(in)            :: unit
        double precision,intent(in)   :: f(0:N, 0:N)
        double precision              :: x, y
        integer                       :: i, j

        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                write (unit, '(*(F10.5, X))') x, y, f(i,j)**2d0, f(i,j), 0d0
            end do
            write (unit, *)
        end do
    end subroutine output_real

    ! Save probability current
    subroutine output_flux(unit, Flux)
        integer,intent(in)          :: unit
        double precision,intent(in) :: Flux(0:N,0:N,1:2)
        double precision            :: x, y
        integer                     :: i, j
        double precision,parameter  :: SCALE = 1d0!1000d0
        
        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                write (unit, *) x, y, SCALE*Flux(i,j,1), SCALE*Flux(i,j,2)
            end do
            write (unit, *)
        end do
    end subroutine

    ! Save rotation of flux vectors
    subroutine output_rotation(unit, Rot)
        integer,intent(in)          :: unit
        double precision,intent(in) :: Rot(0:N,0:N)
        double precision            :: x, y, z
        integer                     :: i, j, k
        double precision,parameter  :: SCALE = 1d0!1000d0
        
        ! Be careful of putting two blank lines. Gnuplot behaves wrongly when specifying every keyword.
        do k = 0, 1
            z = -xmax + dh * k
            do j = 0, N
                y = -xmax + dh * j
                do i = 0, N
                    x = -xmax + dh * i

                    ! X Y Z ROT_X ROT_Y ROT_Z
                    if (k == 0) then
                        write (unit, *) x, y, 0d0, 0d0, 0d0, SCALE*Rot(i,j)
                    else
                        write (unit, *) x, y, z, 0d0, 0d0, 0d0
                    end if
                end do
                write (unit, *)
            end do
        end do
    end subroutine

    subroutine output_widths(unit, iter, width_x, width_y)
        integer,intent(in)          :: unit, iter
        double precision,intent(in) :: width_x, width_y
        write (unit, *) dt*iter, width_x, width_y
    end subroutine
end module io