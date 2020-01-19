! I/O Procedures
module io
    use constants
    implicit none

    interface output
        module procedure output_complex, output_real, output_flux
    end interface
contains
    ! Save double precision complex wave function
    subroutine output_complex(filename, f)
        complex(kind(0d0)),intent(in) :: f(0:N, 0:N)
        character(*),intent(in)       :: filename
        double precision x, y
        integer i, j
        open(10, file=filename)

        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                write (10, '(*(F10.5, X))') x, y, abs(f(i,j))**2d0, dble(f(i,j)), aimag(f(i,j))
            end do
            write (10, *)
        end do

        close(10)
    end subroutine output_complex

    ! Save double precision real wave function
    subroutine output_real(filename, f)
        double precision,intent(in)   :: f(0:N, 0:N)
        character(*),intent(in)       :: filename
        double precision x, y
        integer i, j
        open(11, file=filename)

        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                write (11, '(*(F10.5, X))') x, y, f(i,j)**2d0, f(i,j), 0d0
            end do
            write (11, *)
        end do

        close(11)
    end subroutine output_real
    
    ! Save probability current
    subroutine output_flux(filename, Flux)
        double precision,intent(in) :: Flux(0:N,0:N,1:2)
        character(*),intent(in)     :: filename
        double precision            :: x, y
        integer                     :: i, j
        double precision,parameter  :: SCALE = 1d0
        open(12, file=filename)
        
        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                write (12, *) x, y, SCALE*Flux(i,j,1), SCALE*Flux(i,j,2)
            end do
            write (12, *)
        end do

        close(12)
    end subroutine
end module io