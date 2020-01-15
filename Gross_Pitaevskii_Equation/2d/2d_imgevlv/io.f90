! I/O Procedures
module io
    implicit none
contains
    ! Save double precision complex wave function
    subroutine output_complex(unit, f, N, dh, xmax)
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
        write (unit, *)
    end subroutine output_complex

    ! Save double precision real wave function
    subroutine output_real(unit, f, N, dh, xmax)
        integer,intent(in)            :: unit, N
        double precision,intent(in)   :: dh, xmax
        double precision,intent(in)   :: f(0:N, 0:N)
        double precision x, y
        integer i, j

        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                write (unit, '(*(F10.5, X))') x, y, f(i,j)**2d0, f(i,j), 0d0
            end do
            write (unit, *)
        end do
        write (unit, *)
    end subroutine output_real

    ! Save double precision real data
    subroutine output_potential(unit, f, N, dh, xmax)
        integer,intent(in)            :: unit, N
        double precision,intent(in)   :: dh, xmax
        complex(kind(0d0)),intent(in) :: f(0:N, 0:N)
        double precision x, y
        integer i, j

        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                write (unit, '(*(F10.5, X))') x, y, dble(f(i,j)), aimag(f(i,j))
            end do
            write (unit, *)
        end do
        write (unit, *)
    end subroutine output_potential
    
    ! Save probability current
    subroutine output_flux(unit, Flux, N, dh, xmax)
        integer,intent(in)          :: unit, N
        double precision,intent(in) :: dh, xmax
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
end module io