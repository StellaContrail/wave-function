module io
    implicit none
contains
    ! Output to a file
    subroutine output(unit, f, N, dh, xmax)
        integer,intent(in)            :: unit, N
        double precision,intent(in)   :: dh, xmax
        complex(kind(0d0)),intent(in) :: f(0:N,0:N,0:N)
        double precision x, y, z
        integer i, j, k

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
            write (unit, *)
        end do
        write (unit, *)
    end subroutine output

    ! Output to a file
    subroutine output_projection(unit, f, N, dh, xmax)
        integer,intent(in)            :: unit, N
        double precision,intent(in)   :: dh, xmax
        complex(kind(0d0)),intent(in) :: f(0:N,0:N,0:N)
        double precision              :: prob(0:N,0:N,0:N)
        double precision              :: prob_proj(0:N,0:N)
        double precision x, y, z
        integer i, j, k
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

    ! Output to a file
    subroutine output_projection_real(unit, f, N, dh, xmax)
        integer,intent(in)            :: unit, N
        double precision,intent(in)   :: dh, xmax
        double precision,intent(in)   :: f(0:N,0:N,0:N)
        double precision              :: prob(0:N,0:N,0:N)
        double precision              :: prob_proj(0:N,0:N)
        double precision x, y, z
        integer i, j, k
        prob(:,:,:) = abs(f(:,:,:))**2d0

        prob_proj(:,:) = 0d0
        do k = 0, N
            do j = 0, N
                do i = 0, N
                    if (k == 0 .or. k ==N) then
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
    end subroutine output_projection_real

    ! Output to a file
    subroutine output_cutout(unit, f, N, dh, xmax, index)
        integer,intent(in)            :: unit, N, index
        double precision,intent(in)   :: dh, xmax
        double precision,intent(in) :: f(0:N, 0:N, 0:N)
        double precision x, y, z
        integer i, j, k

        do j = 0, N
            y = -xmax + dh * j
            do i = 0, N
                x = -xmax + dh * i
                write (unit, '(*(F10.5, X))') x, y, f(i,j,index)
            end do
            write (unit, *)
        end do
    end subroutine output_cutout

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