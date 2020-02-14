! I/O Procedures
module io
    use constants
    implicit none

    interface output
        module procedure output_complex, output_real, output_flux, output_complex_unit, output_real_unit, output_flux_unit
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
    subroutine output_complex_unit(unit, f)
        complex(kind(0d0)),intent(in) :: f(0:N, 0:N)
        integer,intent(in)            :: unit
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
    end subroutine output_complex_unit

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
                write (11, *) x, y, f(i,j)**2d0, f(i,j), 0d0
            end do
            write (11, *)
        end do

        close(11)
    end subroutine output_real
    subroutine output_real_unit(unit, f)
        double precision,intent(in)   :: f(0:N, 0:N)
        integer,intent(in)            :: unit
        double precision x, y
        integer i, j
        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                write (unit, *) x, y, f(i,j)**2d0, f(i,j), 0d0
            end do
            write (unit, *)
        end do
    end subroutine output_real_unit

    ! Save potential
    subroutine output_potential(filename, Pot)
        double precision,intent(in)   :: Pot(0:N, 0:N)
        character(*),intent(in)       :: filename
        double precision x, y
        integer i, j
        open(11, file=filename)

        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                write (11, *) x, y, Pot(i,j)*ENERGY_UNIT_IN_DIMENSIONLESS_GPE/(1.603d-19) ! eV
            end do
            write (11, *)
        end do

        close(11)
    end subroutine
    subroutine output_potential_unit(unit, Pot)
        double precision,intent(in)   :: Pot(0:N, 0:N)
        integer,intent(in)            :: unit
        double precision x, y
        integer i, j
        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                write (unit, *) x, y, Pot(i,j)*ENERGY_UNIT_IN_DIMENSIONLESS_GPE/(1.603d-19)
            end do
            write (unit, *)
        end do
    end subroutine
    
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
    subroutine output_flux_unit(unit, Flux)
        double precision,intent(in) :: Flux(0:N,0:N,1:2)
        integer,intent(in)          :: unit
        double precision            :: x, y
        integer                     :: i, j
        double precision,parameter  :: SCALE = 1d0
        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                write (unit, *) x, y, SCALE*Flux(i,j,1), SCALE*Flux(i,j,2)
            end do
            write (unit, *)
        end do
    end subroutine
end module io