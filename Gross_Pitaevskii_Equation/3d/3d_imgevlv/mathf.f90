! Mathematical Procedures
module mathf
    implicit none
contains
    ! Integration
    subroutine integrate(f, N, dh, sum)
        integer,intent(in)           :: N
        double precision,intent(in)  :: f(0:N, 0:N, 0:N), dh
        double precision,intent(out) :: sum
        double precision             :: sum_temp2(0:N, 0:N), sum_temp(0:N)
        integer                      :: i, j, k
        sum = 0d0
        sum_temp(:) = 0d0
        sum_temp2(:,:) = 0d0

        do j = 0, N
            do i = 0, N
                if (i == 0 .or. i == N) then
                    sum_temp2(j,:) = sum_temp2(j,:) + 0.5d0*f(i,j,:)*dh
                else
                    sum_temp2(j,:) = sum_temp2(j,:) + f(i,j,:)*dh
                end if
            end do
        end do

        do j = 0, N
            if (j == 0 .or. j == N) then
                sum_temp(:) = sum_temp(:) + 0.5d0*sum_temp2(j,:)*dh
            else
                sum_temp(:) = sum_temp(:) + sum_temp2(j,:)*dh
            end if
        end do

        do k = 0, N
            if (k == 0 .or. k == N) then
                sum = sum + 0.5d0*sum_temp(k)*dh
            else
                sum = sum + sum_temp(k)*dh
            end if
        end do
    end subroutine

    ! Normalization
    subroutine normalize(f, N, dh)
        integer,intent(in)              :: N
        double precision,intent(in)     :: dh
        double precision,intent(inout)  :: f(0:N, 0:N, 0:N)
        double precision                :: sum
        call integrate(f(:, :, :)**2d0, N, dh, sum)
        f(:, :, :) = f(:, :, :) / sqrt(sum)
    end subroutine normalize

    ! Imaginary-time propagation
    subroutine evolve(Phi_old, N, dt, dh, epsilon, kappa, density, Pot, Phi_next)
        integer,intent(in)           :: N
        double precision,intent(in)  :: dt, dh, epsilon, kappa, density(0:N,0:N,0:N), Pot(0:N,0:N,0:N)
        double precision,intent(in)  :: Phi_old(0:N,0:N,0:N)
        double precision,intent(out) :: Phi_next(0:N,0:N,0:N)
        integer                      :: i
        double precision             :: temp(0:N,0:N,0:N), Atemp(0:N,0:N,0:N)
        ! First term of Taylor expansion
        temp(:,:,:)     = Phi_old(:,:,:)
        Phi_next(:,:,:) = temp(:,:,:)
        ! Other terms of Taylor expansion
        do i = 1, 10
            call apply_hamiltonian(temp, N, dh, epsilon, kappa, density, Pot, Atemp)
            ! Atemp = H*LastTerm
            temp(:,:,:)   = -Atemp(:,:,:)*dt/(epsilon*i)
            Phi_next(:,:,:) = Phi_next(:,:,:) + temp(:,:,:)
        end do
    end subroutine evolve

    ! Calculate HPhi (H:Hamiltonian, Phi:Wave function)
    subroutine apply_hamiltonian(Phi, N, dh, epsilon, kappa, density, Pot, HPhi)
        integer,intent(in)           :: N
        double precision,intent(in)  :: Phi(0:N,0:N,0:N)
        double precision,intent(out) :: HPhi(0:N,0:N,0:N)
        integer                      :: i, j, k
        double precision,intent(in)  :: dh, epsilon, kappa
        double precision,intent(in)  :: Pot(0:N,0:N,0:N), density(0:N,0:N,0:N)
        HPhi(:,:,:) = 0d0
        ! Laplacian part (Five Point Stencil)
        do k = 0, N
            do j = 0, N
                do i = 0, N
                    HPhi(i,j,k) = -90d0*Phi(i,j,k)
                    if (0 < i) then
                        HPhi(i,j,k) = HPhi(i,j,k) + 16d0*Phi(i-1,j,k)
                    end if
                    if (1 < i) then
                        HPhi(i,j,k) = HPhi(i,j,k) - Phi(i-2,j,k)
                    end if
                    if (i < N-1) then
                        HPhi(i,j,k) = HPhi(i,j,k) - Phi(i+2,j,k)
                    end if
                    if (i < N) then
                        HPhi(i,j,k) = HPhi(i,j,k) + 16d0*Phi(i+1,j,k)
                    end if

                    if (0 < j) then
                        HPhi(i,j,k) = HPhi(i,j,k) + 16d0*Phi(i,j-1,k)
                    end if
                    if (1 < j) then
                        HPhi(i,j,k) = HPhi(i,j,k) - Phi(i,j-2,k)
                    end if
                    if (j < N-1) then
                        HPhi(i,j,k) = HPhi(i,j,k) - Phi(i,j+2,k)
                    end if
                    if (j < N) then
                        HPhi(i,j,k) = HPhi(i,j,k) + 16d0*Phi(i,j+1,k)
                    end if

                    if (0 < k) then
                        HPhi(i,j,k) = HPhi(i,j,k) + 16d0*Phi(i,j,k-1)
                    end if
                    if (1 < k) then
                        HPhi(i,j,k) = HPhi(i,j,k) - Phi(i,j,k-2)
                    end if
                    if (k < N-1) then
                        HPhi(i,j,k) = HPhi(i,j,k) - Phi(i,j,k+2)
                    end if
                    if (k < N) then
                        HPhi(i,j,k) = HPhi(i,j,k) + 16d0*Phi(i,j,k+1)
                    end if
                end do
            end do
        end do
        HPhi(:,:,:) = HPhi(:,:,:) / (12d0*dh*dh)
        
        do k = 0, N
            do j = 0, N
                do i = 0, N
                    ! Kinetic energy
                    HPhi(i,j,k) = -0.5d0*epsilon*epsilon*HPhi(i,j,k)

                    ! External potential
                    HPhi(i,j,k) = HPhi(i,j,k) + Pot(i,j,k)*Phi(i,j,k)

                    ! Nonlinear part
                    HPhi(i,j,k) = HPhi(i,j,k) + kappa*density(i,j,k)*Phi(i,j,k)
                end do
            end do
        end do
    end subroutine

    ! Solve Energy Expected Value
    subroutine solve_energy(Phi, Pot, N, epsilon, kappa, mu, dh)
        integer,intent(in)           :: N
        double precision,intent(in)  :: dh, epsilon, kappa
        double precision,intent(in)  :: Pot(0:N,0:N,0:N)
        double precision,intent(in)  :: Phi(0:N,0:N,0:N)
        double precision,intent(out) :: mu
        double precision             :: HPhi(0:N,0:N,0:N), sum_temp2(0:N,0:N), sum_temp(0:N), sum
        integer                      :: i, j, k
        call apply_hamiltonian(Phi, N, dh, epsilon, kappa, abs(Phi)**2d0, Pot, HPhi)

        sum_temp2(:,:) = 0d0
        do j = 0, N
            do i = 0, N
                if (i == 0 .or. i == N) then
                    sum_temp2(j,:) = sum_temp2(j,:) + 0.5d0*Phi(i,j,:)*HPhi(i,j,:)*dh
                else 
                    sum_temp2(j,:) = sum_temp2(j,:) + Phi(i,j,:)*HPhi(i,j,:)*dh
                end if
            end do
        end do
        sum_temp(:) = 0d0
        do j = 0, N
            if (j == 0 .or. j == N) then
                sum_temp(:) = sum_temp(:) + 0.5d0*sum_temp2(j,:)*dh
            else
                sum_temp(:) = sum_temp(:) + sum_temp2(j,:)*dh
            end if
        end do
        sum = 0d0
        do k = 0, N
            if (k == 0 .or. k == N) then
                sum = sum + 0.5d0*sum_temp(k)*dh
            else
                sum = sum + sum_temp(k)*dh
            end if
        end do

        mu = sum
    end subroutine
    
    ! Shift wave fuction's phase partially
    subroutine shift_phase(Phi_IN, N, x_start, x_end, y_start, y_end, z_start, z_end, Phi_OUT, iu, angle)
        integer,intent(in)             :: x_start, x_end, y_start, y_end, z_start, z_end, N
        double precision,intent(in)    :: angle
        double precision,intent(in)    :: Phi_IN(0:N,0:N,0:N)
        complex(kind(0d0)),intent(in)  :: iu
        complex(kind(0d0)),intent(out) :: Phi_OUT(0:N,0:N,0:N)
        Phi_OUT(:,:,:) = Phi_IN(:,:,:)

        Phi_OUT(x_start:x_end,y_start:y_end,z_start:z_end) = exp(iu*angle) * Phi_IN(x_start:x_end,y_start:y_end,z_start:z_end)
    end subroutine
end module mathf