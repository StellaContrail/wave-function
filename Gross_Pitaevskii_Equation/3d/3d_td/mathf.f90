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
        integer,intent(in)                :: N
        double precision,intent(in)       :: dh
        complex(kind(0d0)),intent(inout)  :: f(0:N, 0:N, 0:N)
        double precision                  :: sum
        call integrate(abs(f(:, :, :))**2d0, N, dh, sum)
        f(:, :, :) = f(:, :, :) / sqrt(sum)
    end subroutine normalize

    ! Real-time propagation
    subroutine evolve(Phi_old, N, dt, dh, epsilon, kappa, iu, density, Pot, Phi_next)
        integer,intent(in)             :: N
        double precision,intent(in)    :: dt, dh, epsilon, kappa, density(0:N,0:N,0:N), Pot(0:N,0:N,0:N)
        complex(kind(0d0)),intent(in)  :: Phi_old(0:N,0:N,0:N), iu
        complex(kind(0d0)),intent(out) :: Phi_next(0:N,0:N,0:N)
        integer                        :: i
        complex(kind(0d0))             :: temp(0:N,0:N,0:N), Atemp(0:N,0:N,0:N)
        ! First term of Taylor expansion
        temp(:,:,:)     = Phi_old(:,:,:)
        Phi_next(:,:,:) = temp(:,:,:)
        ! Other terms of Taylor expansion
        do i = 1, 10
            call apply_hamiltonian(temp, N, dh, epsilon, kappa, density, Pot, Atemp)
            temp(:,:,:)   = -iu*Atemp(:,:,:)*dt/(epsilon*i)
            Phi_next(:,:,:) = Phi_next(:,:,:) + temp(:,:,:)
        end do
    end subroutine evolve

    ! Calculate HPhi (H:Hamiltonian, Phi:Wave function)
    subroutine apply_hamiltonian(Phi, N, dh, epsilon, kappa, density, Pot, HPhi)
        integer,intent(in)             :: N
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N,0:N)
        complex(kind(0d0)),intent(out) :: HPhi(0:N,0:N,0:N)
        integer                        :: i, j, k
        double precision,intent(in)    :: dh, epsilon, kappa
        double precision,intent(in)    :: Pot(0:N,0:N,0:N), density(0:N,0:N,0:N)
        HPhi(:,:,:) = dcmplx(0d0, 0d0)
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
        integer,intent(in)             :: N
        double precision,intent(in)    :: dh, epsilon, kappa
        double precision,intent(in)    :: Pot(0:N,0:N,0:N)
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N,0:N)
        double precision,intent(out)   :: mu
        complex(kind(0d0))             :: HPhi(0:N,0:N,0:N), sum_temp2(0:N,0:N), sum_temp(0:N), sum
        integer                        :: i, j, k
        call apply_hamiltonian(Phi, N, dh, epsilon, kappa, abs(Phi)**2d0, Pot, HPhi)

        sum_temp2(:,:) = dcmplx(0d0, 0d0)
        do j = 0, N
            do i = 0, N
                if (i == 0 .or. i == N) then
                    sum_temp2(j,:) = sum_temp2(j,:) + 0.5d0*conjg(Phi(i,j,:))*HPhi(i,j,:)*dh
                else 
                    sum_temp2(j,:) = sum_temp2(j,:) + conjg(Phi(i,j,:))*HPhi(i,j,:)*dh
                end if
            end do
        end do
        sum_temp(:) = dcmplx(0d0, 0d0)
        do j = 0, N
            if (j == 0 .or. j == N) then
                sum_temp(:) = sum_temp(:) + 0.5d0*sum_temp2(j,:)*dh
            else
                sum_temp(:) = sum_temp(:) + sum_temp2(j,:)*dh
            end if
        end do
        sum = dcmplx(0d0, 0d0)
        do k = 0, N
            if (k == 0 .or. k == N) then
                sum = sum + 0.5d0*sum_temp(k)*dh
            else
                sum = sum + sum_temp(k)*dh
            end if
        end do

        mu = dble(sum)
    end subroutine

    ! Calculate Probability Current
    subroutine calc_flux(Phi, N, mass, dh, hbar, Flux)
        integer,intent(in)            :: N
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N,0:N)
        double precision,intent(in)   :: dh, hbar, mass
        double precision,intent(out)  :: Flux(0:N,0:N,0:N,1:3)
        integer                       :: i, j, k

        ! For point (0,0,0)
        Flux(0,0,0,1) = aimag(conjg(Phi(0,0,0))*(Phi(1,0,0)-Phi(0,0,0))/dh)
        Flux(0,0,0,2) = aimag(conjg(Phi(0,0,0))*(Phi(0,1,0)-Phi(0,0,0))/dh)
        Flux(0,0,0,3) = aimag(conjg(Phi(0,0,0))*(Phi(0,0,1)-Phi(0,0,0))/dh)
        ! For point (N,0,0)
        Flux(N,0,0,1) = aimag(conjg(Phi(N,0,0))*(Phi(N,0,0)-Phi(N-1,0,0))/dh)
        Flux(N,0,0,2) = aimag(conjg(Phi(N,0,0))*(Phi(N,1,0)-Phi(N,0,0))/dh)
        Flux(N,0,0,3) = aimag(conjg(Phi(N,0,0))*(Phi(N,0,1)-Phi(N,0,0))/dh)
        ! For point (0,N,0)
        Flux(0,N,0,1) = aimag(conjg(Phi(0,N,0))*(Phi(1,N,0)-Phi(0,N,0))/dh)
        Flux(0,N,0,2) = aimag(conjg(Phi(0,N,0))*(Phi(0,N,0)-Phi(0,N-1,0))/dh)
        Flux(0,N,0,3) = aimag(conjg(Phi(0,N,0))*(Phi(0,N,1)-Phi(0,N,0))/dh)
        ! For point (N,N,0)
        Flux(N,N,0,1) = aimag(conjg(Phi(N,N,0))*(Phi(N,N,0)-Phi(N-1,N,0))/dh)
        Flux(N,N,0,2) = aimag(conjg(Phi(N,N,0))*(Phi(N,N,0)-Phi(N,N-1,0))/dh)
        Flux(N,N,0,3) = aimag(conjg(Phi(N,N,0))*(Phi(N,N,1)-Phi(N,N,0))/dh)
        ! For point (0,0,N)
        Flux(0,0,0,1) = aimag(conjg(Phi(0,0,N))*(Phi(1,0,N)-Phi(0,0,N))/dh)
        Flux(0,0,0,2) = aimag(conjg(Phi(0,0,N))*(Phi(0,1,N)-Phi(0,0,N))/dh)
        Flux(0,0,0,3) = aimag(conjg(Phi(0,0,N))*(Phi(0,0,N)-Phi(0,0,N-1))/dh)
        ! For point (N,0,N)
        Flux(N,0,N,1) = aimag(conjg(Phi(N,0,N))*(Phi(N,0,N)-Phi(N-1,0,N))/dh)
        Flux(N,0,N,2) = aimag(conjg(Phi(N,0,N))*(Phi(N,1,N)-Phi(N,0,N))/dh)
        Flux(N,0,N,3) = aimag(conjg(Phi(N,0,N))*(Phi(N,0,N)-Phi(N,0,N-1))/dh)
        ! For point (0,N,N)
        Flux(0,N,N,1) = aimag(conjg(Phi(0,N,N))*(Phi(1,N,N)-Phi(0,N,N))/dh)
        Flux(0,N,N,2) = aimag(conjg(Phi(0,N,N))*(Phi(0,N,N)-Phi(0,N-1,N))/dh)
        Flux(0,N,N,3) = aimag(conjg(Phi(0,N,N))*(Phi(0,N,N)-Phi(0,N,N-1))/dh)
        ! For point (N,N,N)
        Flux(N,N,N,1) = aimag(conjg(Phi(N,N,N))*(Phi(N,N,N)-Phi(N-1,N,N))/dh)
        Flux(N,N,N,2) = aimag(conjg(Phi(N,N,N))*(Phi(N,N,N)-Phi(N,N-1,N))/dh)
        Flux(N,N,N,3) = aimag(conjg(Phi(N,N,N))*(Phi(N,N,N)-Phi(N,N,N-1))/dh)

        ! For lines on XY Plane of Z=0
        do j = 1, N-1
            ! For range (0,1~N-1,0)
            Flux(0,j,0,1) = aimag(conjg(Phi(0,j,0))*(Phi(1,j,0)-Phi(0,j,0))/dh)
            Flux(0,j,0,2) = aimag(conjg(Phi(0,j,0))*(Phi(0,j+1,0)-Phi(0,j-1,0))/(2d0*dh))
            Flux(0,j,0,3) = aimag(conjg(Phi(0,j,0))*(Phi(0,j,1)-Phi(0,j,0))/dh)

            ! For range (N,1~N-1,0)
            Flux(N,j,0,1) = aimag(conjg(Phi(N,j,0))*(Phi(N,j,0)-Phi(N-1,j,0))/dh)
            Flux(N,j,0,2) = aimag(conjg(Phi(N,j,0))*(Phi(N,j+1,0)-Phi(N,j-1,0))/(2d0*dh))
            Flux(N,j,0,3) = aimag(conjg(Phi(N,j,0))*(Phi(N,j,1)-Phi(N,j,0))/dh)
        end do
        do i = 1, N-1
            ! For range (1~N-1,0,0)
            Flux(i,0,0,1) = aimag(conjg(Phi(i,0,0))*(Phi(i+1,0,0)-Phi(i-1,0,0))/(2d0*dh))
            Flux(i,0,0,2) = aimag(conjg(Phi(i,0,0))*(Phi(i,1,0)-Phi(i,0,0))/dh)
            Flux(i,0,0,3) = aimag(conjg(Phi(i,0,0))*(Phi(i,0,1)-Phi(i,0,0))/dh)

            ! For range (1~N-1,N,0)
            Flux(i,N,0,1) = aimag(conjg(Phi(i,N,0))*(Phi(i+1,N,0)-Phi(i-1,N,0))/(2d0*dh))
            Flux(i,N,0,2) = aimag(conjg(Phi(i,N,0))*(Phi(i,N,0)-Phi(i,N-1,0))/dh)
            Flux(i,N,0,3) = aimag(conjg(Phi(i,N,0))*(Phi(i,N,1)-Phi(i,N,0))/dh)
        end do

        ! For lines on XY Plane of Z=N
        do j = 1, N-1
            ! For range (0,1~N-1,N)
            Flux(0,j,N,1) = aimag(conjg(Phi(0,j,N))*(Phi(1,j,N)-Phi(0,j,N))/dh)
            Flux(0,j,N,2) = aimag(conjg(Phi(0,j,N))*(Phi(0,j+1,N)-Phi(0,j-1,N))/(2d0*dh))
            Flux(0,j,N,3) = aimag(conjg(Phi(0,j,N))*(Phi(0,j,N)-Phi(0,j,N-1))/dh)

            ! For range (N,1~N-1,N)
            Flux(N,j,N,1) = aimag(conjg(Phi(N,j,0))*(Phi(N,j,N)-Phi(N-1,j,N))/dh)
            Flux(N,j,N,2) = aimag(conjg(Phi(N,j,0))*(Phi(N,j+1,N)-Phi(N,j-1,N))/(2d0*dh))
            Flux(N,j,N,3) = aimag(conjg(Phi(N,j,0))*(Phi(N,j,N)-Phi(N,j,N-1))/dh)
        end do
        do i = 1, N-1
            ! For range (1~N-1,0,N)
            Flux(i,0,N,1) = aimag(conjg(Phi(i,0,N))*(Phi(i+1,0,N)-Phi(i-1,0,N))/(2d0*dh))
            Flux(i,0,N,2) = aimag(conjg(Phi(i,0,N))*(Phi(i,1,N)-Phi(i,0,N))/dh)
            Flux(i,0,N,3) = aimag(conjg(Phi(i,0,N))*(Phi(i,0,N)-Phi(i,0,N-1))/dh)

            ! For range (1~N-1,N,N)
            Flux(i,N,N,1) = aimag(conjg(Phi(i,N,N))*(Phi(i+1,N,N)-Phi(i-1,N,N))/(2d0*dh))
            Flux(i,N,N,2) = aimag(conjg(Phi(i,N,N))*(Phi(i,N,N)-Phi(i,N-1,N))/dh)
            Flux(i,N,N,3) = aimag(conjg(Phi(i,N,N))*(Phi(i,N,N)-Phi(i,N,N-1))/dh)
        end do

        ! For lines on YZ Plane of X=0
        do j = 1, N-1
            ! For range (0,0,1~N-1)
            Flux(0,0,j,1) = aimag(conjg(Phi(0,0,j))*(Phi(1,0,j)-Phi(0,0,j))/dh)
            Flux(0,0,j,2) = aimag(conjg(Phi(0,0,j))*(Phi(0,1,j)-Phi(0,0,j))/dh)
            Flux(0,0,j,3) = aimag(conjg(Phi(0,0,j))*(Phi(0,0,j+1)-Phi(0,0,j-1))/(2d0*dh))

            ! For range (0,N,1~N-1)
            Flux(0,N,j,1) = aimag(conjg(Phi(0,N,j))*(Phi(1,N,j)-Phi(0,N,j))/dh)
            Flux(0,N,j,2) = aimag(conjg(Phi(0,N,j))*(Phi(0,N,j)-Phi(0,N-1,j))/dh)
            Flux(0,N,j,3) = aimag(conjg(Phi(0,N,j))*(Phi(0,N,j+1)-Phi(0,N,j-1))/(2d0*dh))
        end do

        ! For lines on YZ Plane of X=N
        do j = 1, N-1
            ! For range (N,0,1~N-1)
            Flux(N,0,j,1) = aimag(conjg(Phi(N,0,j))*(Phi(N,0,j)-Phi(N-1,0,j))/dh)
            Flux(N,0,j,2) = aimag(conjg(Phi(N,0,j))*(Phi(N,1,j)-Phi(N,0,j))/dh)
            Flux(N,0,j,3) = aimag(conjg(Phi(N,0,j))*(Phi(N,0,j+1)-Phi(N,0,j-1))/(2d0*dh))

            ! For range (N,N,1~N-1)
            Flux(N,N,j,1) = aimag(conjg(Phi(N,N,j))*(Phi(N,N,j)-Phi(N-1,N,j))/dh)
            Flux(N,N,j,2) = aimag(conjg(Phi(N,N,j))*(Phi(N,N,j)-Phi(N,N-1,j))/dh)
            Flux(N,N,j,3) = aimag(conjg(Phi(N,N,j))*(Phi(N,N,j+1)-Phi(N,N,j-1))/(2d0*dh))
        end do

        ! For XY Plane without lines
        do j = 1, N-1
            do i = 1, N-1
                ! For Z=0
                Flux(i,j,0,1) = aimag(conjg(Phi(i,j,0))*(Phi(i+1,j,0)-Phi(i-1,j,0))/(2d0*dh))
                Flux(i,j,0,2) = aimag(conjg(Phi(i,j,0))*(Phi(i,j+1,0)-Phi(i,j-1,0))/(2d0*dh))
                Flux(i,j,0,3) = aimag(conjg(Phi(i,j,0))*(Phi(i,j,1)-Phi(i,j,0))/dh)
                ! For Z=N
                Flux(i,j,N,1) = aimag(conjg(Phi(i,j,N))*(Phi(i+1,j,N)-Phi(i-1,j,N))/(2d0*dh))
                Flux(i,j,N,2) = aimag(conjg(Phi(i,j,N))*(Phi(i,j+1,N)-Phi(i,j-1,N))/(2d0*dh))
                Flux(i,j,N,3) = aimag(conjg(Phi(i,j,N))*(Phi(i,j,N)-Phi(i,j,N-1))/dh)
            end do
        end do

        ! For YZ Plane without lines
        do k = 1, N-1
            do j = 1, N-1
                ! For X=0
                Flux(0,j,k,1) = aimag(conjg(Phi(0,j,k))*(Phi(1,j,k)-Phi(0,j,k))/dh)
                Flux(0,j,k,2) = aimag(conjg(Phi(0,j,k))*(Phi(0,j+1,k)-Phi(0,j-1,k))/(2d0*dh))
                Flux(0,j,k,3) = aimag(conjg(Phi(0,j,k))*(Phi(0,j,k+1)-Phi(0,j,k-1))/(2d0*dh))
                ! For X=N
                Flux(N,j,k,1) = aimag(conjg(Phi(N,j,k))*(Phi(N,j,k)-Phi(N-1,j,k))/dh)
                Flux(N,j,k,2) = aimag(conjg(Phi(N,j,k))*(Phi(N,j+1,k)-Phi(N,j-1,k))/(2d0*dh))
                Flux(N,j,k,3) = aimag(conjg(Phi(N,j,k))*(Phi(N,j,k+1)-Phi(N,j,k-1))/(2d0*dh))
            end do
        end do

        ! For ZX Plane without lines
        do k = 1, N-1
            do i = 1, N-1
                ! For Y=0
                Flux(i,0,k,1) = aimag(conjg(Phi(i,0,k))*(Phi(i+1,0,k)-Phi(i-1,0,k))/(2d0*dh))
                Flux(i,0,k,2) = aimag(conjg(Phi(i,0,k))*(Phi(i,1,k)-Phi(i,0,k))/dh)
                Flux(i,0,k,3) = aimag(conjg(Phi(i,0,k))*(Phi(i,0,k+1)-Phi(i,0,k-1))/(2d0*dh))
                ! For Y=N
                Flux(i,0,k,1) = aimag(conjg(Phi(i,N,k))*(Phi(i+1,N,k)-Phi(i-1,N,k))/(2d0*dh))
                Flux(i,0,k,2) = aimag(conjg(Phi(i,N,k))*(Phi(i,N,k)-Phi(i,N-1,k))/dh)
                Flux(i,0,k,3) = aimag(conjg(Phi(i,N,k))*(Phi(i,N,k+1)-Phi(i,N,k-1))/(2d0*dh))
            end do
        end do

        ! For kernel
        do k = 1, N-1
            do j = 1, N-1
                do i = 1, N-1
                    Flux(i,j,k,1) = aimag(conjg(Phi(i,j,k))*(Phi(i+1,j,k)-Phi(i-1,j,k))/(2d0*dh))
                    Flux(i,j,k,2) = aimag(conjg(Phi(i,j,k))*(Phi(i,j+1,k)-Phi(i,j-1,k))/(2d0*dh))
                    Flux(i,j,k,3) = aimag(conjg(Phi(i,j,k))*(Phi(i,j,k+1)-Phi(i,j,k-1))/(2d0*dh))
                end do
            end do
        end do

        !Flux(:,:,:) = Flux(:,:,:)*(hbar/mass)
    end subroutine

    ! Calculate rotation of flux vectors
    subroutine calc_rotation(Flux, N, dh, xmax, Rot)
        integer,intent(in)           :: N
        double precision,intent(in)  :: Flux(0:N,0:N,0:N,1:3), dh, xmax
        double precision,intent(out) :: Rot(0:N,0:N,0:N,1:3)
        integer                      :: i, j, k
        double precision             :: x, y, z

        ! XY Plane
        Rot(:,:,0,:) = 0d0
        Rot(:,:,N,:) = 0d0
        ! YZ Plane
        Rot(0,:,:,:) = 0d0
        Rot(N,:,:,:) = 0d0
        ! ZX Plane
        Rot(:,0,:,:) = 0d0
        Rot(:,N,:,:) = 0d0

        ! For kernel
        do k = 1, N-1
            z = -xmax + dh*k
            do j = 1, N-1
                y = -xmax + dh*j
                do i = 1, N-1
                    x = -xmax + dh*i
                    
                    Rot(i,j,k,1) = (Flux(i,j+1,k,3)-Flux(i,j-1,k,3))/(2d0*dh) - (Flux(i,j,k+1,2)-Flux(i,j,k-1,2))/(2d0*dh)
                    Rot(i,j,k,2) = (Flux(i,j,k+1,1)-Flux(i,j,k-1,1))/(2d0*dh) - (Flux(i+1,j,k,3)-Flux(i-1,j,k,3))/(2d0*dh)
                    Rot(i,j,k,3) = (Flux(i+1,j,k,2)-Flux(i-1,j,k,2))/(2d0*dh) - (Flux(i,j+1,k,1)-Flux(i,j-1,k,1))/(2d0*dh)
                end do
            end do
        end do
    end subroutine


end module mathf