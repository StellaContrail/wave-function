module io
    implicit none
contains
    ! Save double precision complex wave function
    subroutine output_complex(unit, f, N, dh, xmax)
        integer,intent(in)            :: unit, N
        double precision,intent(in)   :: dh, xmax
        complex(kind(0d0)),intent(in) :: f(0:N,0:N,0:N)
        double precision              :: x, y, z
        integer                       :: i, j, k

        do k = 0, N
            z = -xmax + dh * k
            do j = 0, N
                y = -xmax + dh * j
                do i = 0, N
                    x = -xmax + dh * i

                    write (unit, '(*(F10.5, X))') x, y, z, abs(f(i,j,k))**2d0, dble(f(i,j,k)), aimag(f(i,j,k))
                end do
                write (unit, *)
            end do
        end do
    end subroutine output_complex

    subroutine output_potential_cutout(unit, Pot, N, dh, xmax, x_cutout, y_cutout, z_cutout) 
        integer,intent(in)          :: unit, N, x_cutout, y_cutout, z_cutout
        double precision,intent(in) :: dh, Pot(0:N,0:N,0:N), xmax
        integer                     :: i, j, k
        double precision            :: x, y, z

        x = -xmax + dh*x_cutout
        do k = 0, N
            z = -xmax + dh*k
            do j = 0, N
                y = -xmax + dh*j

                write (unit, '(*(F10.5, X))') x, y, z, Pot(x_cutout,j,k)
            end do
            write (unit, *)
        end do

        y = -xmax + dh*y_cutout
        do k = 0, N
            z = -xmax + dh*k
            do i = 0, N
                x = -xmax + dh*i

                write (unit, '(*(F10.5, X))') x, y, z, Pot(i,y_cutout,k)
            end do
            write (unit, *)
        end do

        z = -xmax + dh*z_cutout
        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i

                write (unit, '(*(F10.5, X))') x, y, z, Pot(i,j,z_cutout)
            end do
            write (unit, *)
        end do
    end subroutine

    subroutine output_potential_raw(unit, Pot, N, dh, xmax) 
        integer,intent(in)          :: unit, N
        double precision,intent(in) :: dh, Pot(0:N,0:N,0:N), xmax
        integer                     :: i, j, k
        double precision            :: x, y, z

        do k = 0, N
            z = -xmax + dh*k
            do j = 0, N
                y = -xmax + dh*j
                do i = 0, N
                    x = -xmax + dh*i

                    write (unit, '(*(F10.5, X))') x, y, z, Pot(i,j,k)
                end do
                write (unit, *)
            end do
        end do
    end subroutine

    ! Save double precision complex wave function
    subroutine output_real(unit, f, N, dh, xmax)
        integer,intent(in)          :: unit, N
        double precision,intent(in) :: dh, xmax
        double precision,intent(in) :: f(0:N,0:N,0:N)
        double precision            :: x, y, z
        integer                     :: i, j, k

        do k = 0, N
            z = -xmax + dh * k
            do j = 0, N
                y = -xmax + dh * j
                do i = 0, N
                    x = -xmax + dh * i

                    write (unit, '(*(F10.5, X))') x, y, z, f(i,j,k)**2d0, f(i,j,k), 0d0
                end do
                write (unit, *)
            end do
        end do
    end subroutine output_real

    ! Save projection of wave function
    subroutine output_projection(unit, f, N, dh, xmax)
        integer,intent(in)          :: unit, N
        double precision,intent(in) :: dh, xmax
        double precision,intent(in) :: f(0:N,0:N,0:N)
        double precision            :: prob(0:N,0:N,0:N)
        double precision            :: prob_proj(0:N,0:N)
        double precision            :: x, y
        integer                     :: i, j, k
        prob(:,:,:) = abs(f(:,:,:))**2d0

        prob_proj(:,:) = 0d0
        do k = 0, N
            do j = 0, N
                do i = 0, N
                    if (k == 0 .or. k == N) then
                        prob_proj(i,j) = prob_proj(i,j) + 0.5d0*prob(i,j,k)*dh
                    else 
                        prob_proj(i,j) = prob_proj(i,j) + prob(i,j,k)*dh
                    end if
                end do
            end do
        end do

        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                
                write (unit, '(*(F10.5, X))') x, y, prob_proj(i,j)
            end do
            write (unit, *)
        end do
    end subroutine output_projection

    ! Save cutout plane of wave function
    subroutine output_cutout(unit, f, N, dh, xmax, index)
        integer,intent(in)          :: unit, N, index
        double precision,intent(in) :: dh, xmax
        double precision,intent(in) :: f(0:N, 0:N, 0:N)
        double precision            :: x, y
        integer                     :: i, j

        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i

                write (unit, '(*(F10.5, X))') x, y, f(i,j,index)
            end do
            write (unit, *)
        end do
    end subroutine output_cutout
end module io