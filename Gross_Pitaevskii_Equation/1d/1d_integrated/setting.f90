! Set up wave functions and potential
module setting
    use constants
    use mathf
    implicit none
contains
    ! Initialize wave functions and potential
    subroutine initialize(Pot, mode, Phi)
        complex(kind(0d0)),intent(out),optional  :: Phi(0:N)
        double precision,intent(out)             :: Pot(0:N)
        integer,intent(in)                       :: mode
        double precision                         :: x    ! x~, y~
        integer                                  :: i

        do i = 0, N
            x = -xmax + dh*i
            if (present(Phi)) then
                ! Initial trial wave function
                Phi(i) = exp(-0.5d0*x*x)
            end if

            ! External potential
            select case (mode)
            case (0)
                ! Harmonic Oscillator Trap
                Pot(i) = 0.5d0*x**2d0
            case (1)
                ! Harmonic Oscillator Trap and Very Narrow Gaussian-shaped Wall at the center
                Pot(i) = 0.5d0*x**2d0
                Pot(i) = Pot(i) + 100d0*exp(-0.5d0*x*x/0.5d0**2d0)
            case (2)
                ! Box Trap
                if (abs(x) < 3d0) then
                    Pot(i) = -50d0
                else
                    Pot(i) = 0d0
                end if
            case (3)
                ! Free particle
                Pot(i) = 0d0
            case default
                stop "Invalid mode of external potential"
            end select
        end do
    end subroutine initialize
end module