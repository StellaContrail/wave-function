! I/O Procedures
module io
    use constants
    implicit none

    interface output
        module procedure output_complex, output_real, output_complex_unit, output_real_unit
    end interface
contains
    ! Save double precision complex wave function
    subroutine output_complex(filename, f)
        complex(kind(0d0)),intent(in) :: f(0:N)
        character(*),intent(in)       :: filename
        double precision x
        integer i
        open(10, file=filename)

        do i = 0, N
            x = -xmax + dh * i
            write (10, '(*(F10.5, X))') x, abs(f(i))**2d0, dble(f(i)), aimag(f(i))
        end do
        write (10, *)

        close(10)
    end subroutine output_complex
    subroutine output_complex_unit(unit, f)
        complex(kind(0d0)),intent(in) :: f(0:N)
        integer,intent(in)            :: unit
        double precision x
        integer i
        do i = 0, N
            x = -xmax + dh * i
            write (unit, '(*(F10.5, X))') x, abs(f(i))**2d0, dble(f(i)), aimag(f(i))
        end do
        write (unit, *)
    end subroutine output_complex_unit

    ! Save double precision real wave function
    subroutine output_real(filename, f)
        double precision,intent(in)   :: f(0:N)
        character(*),intent(in)       :: filename
        double precision x
        integer i
        open(11, file=filename)

        do i = 0, N
            x = -xmax + dh * i
            write (11, *) x, f(i)**2d0, f(i), 0d0
        end do
        write (11, *)

        close(11)
    end subroutine output_real
    subroutine output_real_unit(unit, f)
        double precision,intent(in)   :: f(0:N)
        integer,intent(in)            :: unit
        double precision x
        integer i

        do i = 0, N
            x = -xmax + dh * i
            write (unit, *) x, f(i)**2d0, f(i), 0d0
        end do
        write (unit, *)
    end subroutine output_real_unit

    ! Save potential
    subroutine output_potential(filename, Pot)
        double precision,intent(in)   :: Pot(0:N)
        character(*),intent(in)       :: filename
        double precision x
        integer i
        open(11, file=filename)

        do i = 0, N
            x = -xmax + dh * i
            write (11, *) x, Pot(i)
        end do
        write (11, *)

        close(11)
    end subroutine
    subroutine output_potential_unit(unit, Pot)
        double precision,intent(in)   :: Pot(0:N)
        integer,intent(in)            :: unit
        double precision x
        integer i

        do i = 0, N
            x = -xmax + dh * i
            write (unit, *) x, Pot(i)
        end do
        write (unit, *)
    end subroutine
    
    ! Save probability current
    subroutine output_flux(filename, Flux)
        double precision,intent(in) :: Flux(0:N)
        character(*),intent(in)     :: filename
        double precision            :: x
        integer                     :: i
        double precision,parameter  :: SCALE = 1d0
        open(12, file=filename)
        
        do i = 0, N
            x = -xmax + dh * i
            write (12, *) x, SCALE*Flux(i)
        end do
        write (12, *)

        close(12)
    end subroutine
    subroutine output_flux_unit(unit, Flux)
        double precision,intent(in) :: Flux(0:N)
        integer,intent(in)          :: unit
        double precision            :: x
        integer                     :: i
        double precision,parameter  :: SCALE = 1d0

        do i = 0, N
            x = -xmax + dh * i
            write (unit, *) x, SCALE*Flux(i)
        end do
        write (unit, *)
    end subroutine
end module io