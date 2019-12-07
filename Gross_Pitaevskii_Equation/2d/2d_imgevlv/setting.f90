module setting
    implicit none
contains
    ! Initialize the wave functions
    ! Phi_next : Wave function of next space
    ! Phi_prev : Wave function of previous space
    ! Pot      : Potential function
    ! N        : Dimension of space exclusing the first element
    ! dh       : Step distance of space
    ! xmax     : Max x position
  subroutine initialize(Phi_next, Phi_prev, Pot, N, dh, xmax, gamma)
      integer,intent(in)              :: N
      complex(kind(0d0)),intent(out)  :: Phi_next(0:N, 0:N), Phi_prev(0:N, 0:N)
      double precision,intent(out)    :: Pot(0:N, 0:N)
      double precision,intent(in)     :: dh, xmax, gamma
      double precision                :: x, y
      double precision,parameter      :: sigma = 0.5d0
      integer                         :: i,j
      double precision                :: SCALE = 0.1d0
      integer,parameter               :: mode = 4
      Phi_next(:, :) = dcmplx(0d0, 0d0)

      do j = 0, N
         y = -xmax + dh*j
         do i = 0, N
            x = -xmax + dh*i

            ! External potential
            if (mode == 1) then
                Pot(i, j) = 0.5d0*(x*x+gamma*gamma*y*y)*SCALE! + 100d0*exp(-0.5d0*(x*x+y*y)/(sigma*sigma))
            else if (mode == 2) then
                if (abs(x) < 10d0 .and. abs(y) < 10d0) then
                    Pot(i, j) = -5d0
                else
                    Pot(i, j) = 0d0
                end if
            else if (mode == 3) then
                if (sqrt(x**2d0+y**2d0) < 10d0) then
                    Pot(i, j) = -5d0
                else
                    Pot(i, j) = 0d0
                end if
            else if (mode == 4) then
                Pot(i, j) = 0.5d0*(x*x*2d0+gamma*gamma*y*y*0.06d0)*SCALE
                if (abs(y) < 0.4d0 .and. abs(x) < 5d0) then
                    Pot(i, j) = 100d0
                end if
            end if

            ! Initial trial wave function
            Phi_prev(i, j) = dcmplx(exp(-0.5d0*(x*x+y*y)), 0d0)
         end do
      end do
  end subroutine initialize
end module