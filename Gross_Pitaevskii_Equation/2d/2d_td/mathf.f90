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

    ! Normalization
    subroutine normalize(f, N, dh)
        integer,intent(in)                :: N
        double precision,intent(in)       :: dh
        complex(kind(0d0)),intent(inout)  :: f(0:N, 0:N)
        double precision                  :: sum
        call integrate(abs(f(:, :))**2d0, N, dh, sum)
        f(:, :) = f(:, :) / sqrt(sum)
    end subroutine normalize

    ! Real-time propagation
    subroutine evolve(Phi_old, N, dt, dh, epsilon, kappa, iu, density, Pot, Phi_next)
        integer,intent(in)             :: N
        double precision,intent(in)    :: dt, dh, epsilon, kappa, density(0:N,0:N), Pot(0:N,0:N)
        complex(kind(0d0)),intent(in)  :: Phi_old(0:N,0:N), iu
        complex(kind(0d0)),intent(out) :: Phi_next(0:N,0:N)
        integer                        :: i
        complex(kind(0d0))             :: temp(0:N,0:N), Atemp(0:N,0:N)
        ! First term of Taylor expansion
        temp(:,:)     = Phi_old(:,:)
        Phi_next(:,:) = temp(:,:)
        ! Other terms of Taylor expansion
        do i = 1, 20
            call apply_hamiltonian(temp, N, dh, epsilon, kappa, density, Pot, Atemp)
            temp(:,:)   = -iu*Atemp(:,:)*dt/(epsilon*i)
            Phi_next(:,:) = Phi_next(:,:) + temp(:,:)
        end do
    end subroutine evolve

    ! Calculate HPhi (H:Hamiltonian, Phi:Wave function)
    subroutine apply_hamiltonian(Phi, N, dh, epsilon, kappa, density, Pot, HPhi)
        integer,intent(in)             :: N
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N)
        complex(kind(0d0)),intent(out) :: HPhi(0:N,0:N)
        integer                        :: i, j
        double precision,intent(in)    :: dh, epsilon, kappa
        double precision,intent(in)    :: Pot(0:N,0:N), density(0:N,0:N)
        HPhi(:,:) = dcmplx(0d0, 0d0)
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
        integer,intent(in)             :: N
        double precision,intent(in)    :: dh, epsilon, kappa
        double precision,intent(in)    :: Pot(0:N,0:N)
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N)
        double precision,intent(out)   :: mu
        complex(kind(0d0))             :: HPhi(0:N,0:N), sum_temp(0:N), sum
        integer                        :: i
        call apply_hamiltonian(Phi, N, dh, epsilon, kappa, abs(Phi)**2d0, Pot, HPhi)

        sum_temp(:) = dcmplx(0d0, 0d0)
        do i = 0, N
            if (i == 0 .or. i == N) then
                sum_temp(:) = sum_temp(:) + 0.5d0*conjg(Phi(i,:))*HPhi(i,:)*dh
            else 
                sum_temp(:) = sum_temp(:) + conjg(Phi(i,:))*HPhi(i,:)*dh
            end if
        end do
        sum = dcmplx(0d0, 0d0)
        do i = 0, N
            if (i == 0 .or. i == N) then
                sum = sum + 0.5d0*sum_temp(i)*dh
            else
                sum = sum + sum_temp(i)*dh
            end if
        end do

        mu = dble(sum)
    end subroutine

    ! Calculate Probability Current
    subroutine calc_flux(Phi, N, mass, dh, hbar, Flux)
        integer,intent(in)            :: N
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N)
        double precision,intent(in)   :: dh, hbar, mass
        double precision,intent(out)  :: Flux(0:N,0:N,1:2)
        integer                       :: i, j

        ! For point (0,0)
        Flux(0,0,1) = aimag(conjg(Phi(0,0))*(Phi(1,0)-Phi(0,0))/dh)
        Flux(0,0,2) = aimag(conjg(Phi(0,0))*(Phi(0,1)-Phi(0,0))/dh)
        ! For point (0,N)
        Flux(0,N,1) = aimag(conjg(Phi(0,N))*(Phi(1,N)-Phi(0,N))/dh)
        Flux(0,N,2) = aimag(conjg(Phi(0,N))*(Phi(0,N)-Phi(0,N-1))/dh)
        ! For point (N,0)
        Flux(N,0,1) = aimag(conjg(Phi(N,0))*(Phi(N,0)-Phi(N-1,0))/dh)
        Flux(N,0,2) = aimag(conjg(Phi(N,0))*(Phi(N,1)-Phi(N,0))/dh)
        ! For point (N,N)
        Flux(N,N,1) = aimag(conjg(Phi(N,N))*(Phi(N,N)-Phi(N-1,N))/dh)
        Flux(N,N,2) = aimag(conjg(Phi(N,N))*(Phi(N,N)-Phi(N,N-1))/dh)

        do j = 1, N-1
            Flux(0,j,1) = aimag(conjg(Phi(0,j))*(Phi(1,j)-Phi(0,j))/dh)
            Flux(0,j,2) = aimag(conjg(Phi(0,j))*(Phi(0,j+1)-Phi(0,j-1))/(2d0*dh))
            Flux(N,j,1) = aimag(conjg(Phi(N,j))*(Phi(N,j)-Phi(N-1,j))/dh)
            Flux(N,j,2) = aimag(conjg(Phi(N,j))*(Phi(N,j+1)-Phi(N,j-1))/(2d0*dh))
        end do

        do i = 1, N-1
            Flux(i,0,1) = aimag(conjg(Phi(i,0))*(Phi(i+1,0)-Phi(i-1,0))/(2d0*dh))
            Flux(i,0,2) = aimag(conjg(Phi(i,0))*(Phi(i,1)-Phi(i,0))/dh)
            Flux(i,N,1) = aimag(conjg(Phi(i,N))*(Phi(i+1,N)-Phi(i-1,N))/(2d0*dh))
            Flux(i,N,2) = aimag(conjg(Phi(i,N))*(Phi(i,N)-Phi(i,N-1))/dh)
        end do

        do j = 1, N-1
            do i = 1, N-1
                Flux(i,j,1) = aimag(conjg(Phi(i,j))*(Phi(i+1,j)-Phi(i-1,j))/(2d0*dh))
                Flux(i,j,2) = aimag(conjg(Phi(i,j))*(Phi(i,j+1)-Phi(i,j-1))/(2d0*dh))
            end do
        end do

        !Flux(:,:,:) = Flux(:,:,:)*(hbar/mass)
    end subroutine

    ! Calculate rotation of flux vectors
    subroutine calc_rotation(Flux, N, dh, xmax, Rot)
        integer,intent(in)           :: N
        double precision,intent(in)  :: Flux(0:N,0:N,1:2), dh, xmax
        double precision,intent(out) :: Rot(0:N,0:N)
        integer                      :: i, j
        double precision             :: x, y

        Rot(0,0) = (Flux(1,0,2)-Flux(0,0,2))/dh - (Flux(0,1,1)-Flux(0,0,1))/dh
        Rot(0,N) = (Flux(1,N,2)-Flux(0,N,2))/dh - (Flux(0,N,1)-Flux(0,N-1,1))/dh
        Rot(N,0) = (Flux(N,0,2)-Flux(N-1,0,2))/dh - (Flux(N,1,1)-Flux(N,0,1))/dh
        Rot(N,N) = (Flux(N,N,2)-Flux(N-1,N,2))/dh - (Flux(N,N,1)-Flux(N,N-1,1))/dh

        do j = 1, N-1
            Rot(0,j) = (Flux(1,j,2)-Flux(0,j,2))/dh - (Flux(0,j+1,1)-Flux(0,j-1,1))/(2d0*dh)
            Rot(N,j) = (Flux(N,j,2)-Flux(N-1,j,2))/dh - (Flux(N,j+1,1)-Flux(N,j-1,1))/(2d0*dh)
        end do

        do i = 1, N-1
            Rot(i,0) = (Flux(i+1,0,2)-Flux(i-1,0,2))/(2d0*dh) - (Flux(i,1,1)-Flux(i,0,1))/dh
            Rot(i,N) = (Flux(i+1,N,2)-Flux(i-1,N,2))/(2d0*dh) - (Flux(i,N,1)-Flux(i,N-1,1))/dh
        end do

        do j = 1, N-1
            y = -xmax + dh*j
            do i = 1, N-1
                x = -xmax + dh*i
                Rot(i,j) = (Flux(i+1,j,2)-Flux(i-1,j,2))/(2d0*dh) - (Flux(i,j+1,1)-Flux(i,j-1,1))/(2d0*dh)
            end do
        end do
    end subroutine

    ! Calculate LzPhi (Excluding Plank constant)
    subroutine calc_angular_momentum(Phi, N, xmax, dh, iu, LzPhi)
        integer,intent(in)             :: N
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N), iu
        double precision,intent(in)    :: dh, xmax
        complex(kind(0d0)),intent(out) :: LzPhi(0:N,0:N)
        integer                        :: i, j
        double precision               :: x, y
        ! Assume the existence probability at the edge equals zero
        LzPhi(0,0:N) = dcmplx(0d0, 0d0)
        LzPhi(N,0:N) = dcmplx(0d0, 0d0)
        LzPhi(0:N,0) = dcmplx(0d0, 0d0)
        LzPhi(0:N,N) = dcmplx(0d0, 0d0)

        do j = 1, N-1
            y = -xmax + dh*j
            do i = 1, N-1
                x = -xmax + dh*i

                LzPhi(i,j) = -iu*(x*(Phi(i,j+1)-Phi(i,j-1)) - y*(Phi(i+1,j)-Phi(i-1,j)))/(2d0*dh)
            end do
        end do
        LzPhi(:,:) = LzPhi(:,:)
    end subroutine

    ! Calculate Expected Angular Momentum Value (Excluding Plank constant)
    subroutine calc_angular_momentum_expected_value(Phi, N, dh, LzPhi, Lz)
        integer,intent(in)             :: N
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N)
        double precision,intent(in)    :: dh
        complex(kind(0d0)),intent(out) :: LzPhi(0:N,0:N)
        integer                        :: i
        double precision               :: Lz_temp
        double precision,intent(out)   :: Lz
        complex(kind(0d0))             :: sum_temp(0:N), sum_cmplx

        sum_temp(:) = dcmplx(0d0, 0d0)
        do i = 0, N
            if (i == 0 .or. i == N) then
                sum_temp(:) = sum_temp(:) + 0.5d0*conjg(Phi(i,:))*LzPhi(i,:)*dh
            else 
                sum_temp(:) = sum_temp(:) + conjg(Phi(i,:))*LzPhi(i,:)*dh
            end if
        end do
        sum_cmplx = dcmplx(0d0, 0d0)
        do i = 0, N
            if (i == 0 .or. i == N) then
                sum_cmplx = sum_cmplx + 0.5d0*sum_temp(i)*dh
            else
                sum_cmplx = sum_cmplx + sum_temp(i)*dh
            end if
        end do
        ! Check wether sum is almost real here (Not implemented yet)
        Lz = dble(sum_cmplx)
        if (Lz > 1d-6) then
            write (*, *) "Angular momentum is not properly calculated"
        end if
    end subroutine  
    
    subroutine calc_widths(Phi, N, dh, xmax, sigma_x, sigma_y)
        integer,intent(in)            :: N
        double precision,intent(in)   :: dh, xmax
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N)
        double precision,intent(out)  :: sigma_x, sigma_y
        double precision              :: x, y, density, EX, EY, EX2, EY2
        integer                       :: i, j
        EX  = 0d0
        EY  = 0d0
        EX2 = 0d0
        EY2 = 0d0
        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i

                density = abs(Phi(i,j))**2d0
                if (i == 0 .or. i == N) then
                    EX  = EX + 0.5d0*x*density*dh
                    EX2 = EX2 + 0.5d0*x*x*density*dh 
                else
                    EX  = EX + x*density*dh
                    EX2 = EX2 + x*x*density*dh
                end if
                if (j == 0 .or. j == N) then
                    EY  = EY + 0.5d0*y*density*dh
                    EY2 = EY2 + 0.5d0*y*y*density*dh
                else
                    EY  = EY + y*density*dh
                    EY2 = EY2 + y*y*density*dh
                end if
            end do
        end do

        sigma_x = sqrt(EX2 - EX**2d0)
        sigma_y = sqrt(EY2 - EY**2d0)
    end subroutine
end module mathf