! Mathematical Procedures
module mathf
    implicit none
contains
    ! Integration
    subroutine integrate(f, N, dh, sum)
        integer,intent(in)           :: N
        double precision,intent(in)  :: f(0:N, 0:N), dh
        double precision,intent(out) :: sum
        double precision             :: sum_temp(0:N)
        integer                      :: i
        sum = 0d0
        sum_temp(:) = 0d0
        do i = 0, N
            if (i == 0 .or. i == N) then
                sum_temp(:) = sum_temp(:) + 0.5d0*f(i,:)*dh
            else
                sum_temp(:) = sum_temp(:) + f(i,:)*dh
            end if
        end do
        do i = 0, N
            if (i == 0 .or. i == N) then
                sum = sum + 0.5d0*sum_temp(i)*dh
            else
                sum = sum + sum_temp(i)*dh
            end if
        end do
    end subroutine

    ! Normalize
    subroutine normalize(f, N, dh)
        integer,intent(in)             :: N
        double precision,intent(in)    :: dh
        double precision,intent(inout) :: f(0:N, 0:N)
        double precision               :: sum
        call integrate(f(:, :)**2d0, N, dh, sum)
        f(:, :) = f(:, :) / sqrt(sum)
    end subroutine normalize

    ! Imaginary-time propagation
    subroutine evolve(Phi_old, N, dt, dh, epsilon, kappa, density, Pot, Phi_next)
        integer,intent(in)           :: N
        double precision,intent(in)  :: dt, dh, epsilon, kappa, density(0:N,0:N), Pot(0:N,0:N)
        double precision,intent(in)  :: Phi_old(0:N,0:N)
        double precision,intent(out) :: Phi_next(0:N,0:N)
        double precision             :: temp(0:N,0:N), Atemp(0:N,0:N)
        integer                      :: i
        ! First term of Taylor expansion
        temp(:,:)     = Phi_old(:,:)
        Phi_next(:,:) = temp(:,:)
        ! Other terms of Taylor expansion
        do i = 1, 10
            call apply_hamiltonian(temp, N, dh, epsilon, kappa, density, Pot, Atemp)
            temp(:,:)     = -Atemp(:,:)*dt/(epsilon*i)
            Phi_next(:,:) = Phi_next(:,:) + temp(:,:)
        end do
    end subroutine evolve

    ! Calculate HPhi (H:Hamiltonian, Phi:Wave function)
    subroutine apply_hamiltonian(Phi, N, dh, epsilon, kappa, density, Pot, HPhi)
        integer,intent(in)           :: N
        double precision,intent(in)  :: Phi(0:N,0:N)
        double precision,intent(out) :: HPhi(0:N,0:N)
        double precision,intent(in)  :: dh, epsilon, kappa
        double precision,intent(in)  :: Pot(0:N,0:N), density(0:N,0:N)
        integer                      :: i, j
        HPhi(:,:) = 0d0
        ! Laplacian part (Five Point Stencil)
        do j = 0, N
            do i = 0, N
                HPhi(i,j) = -60d0*Phi(i,j)
                if (0 < i) then
                    HPhi(i,j) = HPhi(i,j) + 16d0*Phi(i-1,j)
                end if
                if (1 < i) then
                    HPhi(i,j) = HPhi(i,j) - Phi(i-2,j)
                end if
                if (i < N-1) then
                    HPhi(i,j) = HPhi(i,j) - Phi(i+2,j)
                end if
                if (i < N) then
                    HPhi(i,j) = HPhi(i,j) + 16d0*Phi(i+1,j)
                end if

                if (0 < j) then
                    HPhi(i,j) = HPhi(i,j) + 16d0*Phi(i,j-1)
                end if
                if (1 < j) then
                    HPhi(i,j) = HPhi(i,j) - Phi(i,j-2)
                end if
                if (j < N-1) then
                    HPhi(i,j) = HPhi(i,j) - Phi(i,j+2)
                end if
                if (j < N) then
                    HPhi(i,j) = HPhi(i,j) + 16d0*Phi(i,j+1)
                end if
            end do
        end do
        HPhi(:,:) = HPhi(:,:) / (12d0*dh*dh)
        
        do j = 0, N
            do i = 0, N
                ! Kinetic energy
                HPhi(i,j) = -0.5d0*epsilon*epsilon*HPhi(i,j)

                ! External potential
                HPhi(i,j) = HPhi(i,j) + Pot(i,j)*Phi(i,j)

                ! Nonlinear part
                HPhi(i,j) = HPhi(i,j) + kappa*density(i,j)*Phi(i,j)
            end do
        end do
    end subroutine

    ! Solve Energy Expected Value
    subroutine solve_energy(Phi, Pot, N, epsilon, kappa, mu, dh)
        integer,intent(in)           :: N
        double precision,intent(in)  :: dh, epsilon, kappa
        double precision,intent(in)  :: Pot(0:N,0:N)
        double precision,intent(in)  :: Phi(0:N,0:N)
        double precision,intent(out) :: mu
        double precision             :: HPhi(0:N,0:N), sum_temp(0:N), sum
        integer                      :: i
        mu = 0d0
        call apply_hamiltonian(Phi, N, dh, epsilon, kappa, Phi**2d0, Pot, HPhi)
        sum_temp(:) = 0d0
        do i = 0, N
            if (i == 0 .or. i == N) then
                sum_temp(:) = sum_temp(:) + 0.5d0*Phi(i,:)*HPhi(i,:)*dh
            else 
                sum_temp(:) = sum_temp(:) + Phi(i,:)*HPhi(i,:)*dh
            end if
        end do
        sum = 0d0
        do i = 0, N
            if (i == 0 .or. i == N) then
                sum = sum + 0.5d0*sum_temp(i)*dh
            else
                sum = sum + sum_temp(i)*dh
            end if
        end do
        mu = sum
    end subroutine

    ! Shift wave fuction's phase partially
    subroutine shift_phase(Phi_IN, N, x_start, x_end, y_start, y_end, Phi_OUT, iu, angle)
        integer,intent(in)             :: x_start, x_end, y_start, y_end, N
        double precision,intent(in)    :: angle
        double precision,intent(in)    :: Phi_IN(0:N,0:N)
        complex(kind(0d0)),intent(in)  :: iu
        complex(kind(0d0)),intent(out) :: Phi_OUT(0:N,0:N)
        integer                        :: i
        Phi_OUT(:,:) = Phi_IN(:,:)
        do i = x_start, x_end
            Phi_OUT(i,y_start:y_end) = exp(iu*angle) * Phi_IN(i,y_start:y_end)
        end do
    end subroutine
end module mathf