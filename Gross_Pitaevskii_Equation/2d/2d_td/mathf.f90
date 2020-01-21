! Mathematical Procedures
module mathf
    use constants
    implicit none

    interface integrate
        module procedure integrate_real, integrate_complex
    end interface
contains
    ! Integration
    double precision function integrate_real(f)
        double precision,intent(in)  :: f(0:N, 0:N)
        double precision             :: sum
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
        integrate_real = sum
    end function
    complex(kind(0d0)) function integrate_complex(f)
        complex(kind(0d0)),intent(in)  :: f(0:N, 0:N)
        complex(kind(0d0))             :: sum
        complex(kind(0d0))             :: sum_temp(0:N)
        integer                        :: i
        sum = dcmplx(0d0, 0d0)
        sum_temp(:) = dcmplx(0d0, 0d0)
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
        integrate_complex = sum
    end function

    ! Normalization
    subroutine normalize(f)
        complex(kind(0d0)),intent(inout)  :: f(0:N, 0:N)
        f(:, :) = f(:, :) / sqrt(integrate(abs(f(:, :))**2d0))
    end subroutine normalize

    ! Real-time propagation
    subroutine evolve(Phi, Pot, LzPhi)
        double precision,intent(in)      :: Pot(0:N,0:N)
        complex(kind(0d0)),intent(in)    :: LzPhi(0:N,0:N) 
        complex(kind(0d0)),intent(inout) :: Phi(0:N,0:N)
        complex(kind(0d0))               :: Phi_old(0:N,0:N), Phi_new(0:N,0:N)
        integer                          :: i
        double precision                 :: density(0:N,0:N)
        complex(kind(0d0))               :: temp(0:N,0:N), Atemp(0:N,0:N)

        Phi_old(:,:) = Phi(:,:)
        ! First term of Taylor expansion
        temp(:,:)    = Phi_old(:,:)
        Phi_new(:,:) = temp(:,:)
        ! Other terms of Taylor expansion
        do i = 1, 10
            call apply_hamiltonian(temp, abs(Phi_old)**2d0, Pot, LzPhi, Atemp)
            temp(:,:)   = -iu*Atemp(:,:)*dt/(epsilon*i)
            Phi_new(:,:) = Phi_new(:,:) + temp(:,:)
        end do
        
        density(:,:) = 0.5d0*abs(Phi(:,:))**2d0 + 0.5d0*abs(Phi_new(:,:))**2d0

        ! First term of Taylor expansion
        temp(:,:)    = Phi_old(:,:)
        Phi_new(:,:) = temp(:,:)
        ! Other terms of Taylor expansion
        do i = 1, 10
            call apply_hamiltonian(temp, density, Pot, LzPhi, Atemp)
            temp(:,:)   = -iu*Atemp(:,:)*dt/(epsilon*i)
            Phi_new(:,:) = Phi_new(:,:) + temp(:,:)
        end do

        Phi(:,:) = Phi_new(:,:)
    end subroutine evolve

    ! Calculate HPhi (H:Hamiltonian, Phi:Wave function)
    subroutine apply_hamiltonian(Phi, density, Pot, LzPhi, HPhi)
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N), LzPhi(0:N,0:N)
        complex(kind(0d0)),intent(out) :: HPhi(0:N,0:N)
        integer                        :: i, j
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

                ! Cranking model
                HPhi(i,j) = HPhi(i,j) - OMEGA*LzPhi(i,j)/omega_x
            end do
        end do
    end subroutine

    ! Solve Energy Expected Value
    double precision function solve_energy(Phi, Pot, LzPhi)
        double precision,intent(in)    :: Pot(0:N,0:N)
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N), LzPhi(0:N,0:N)
        complex(kind(0d0))             :: HPhi(0:N,0:N)
        call apply_hamiltonian(Phi, abs(Phi)**2d0, Pot, LzPhi, HPhi)
        solve_energy = dble(integrate(conjg(Phi)*HPhi))
    end function

    ! Calculate Probability Current
    function calc_flux(Phi) result(Flux)
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N)
        double precision              :: Flux(0:N,0:N,1:2)
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
    end function

    ! Calculate rotation of flux vectors
    subroutine calc_rotation(Flux, Rot)
        double precision,intent(in)  :: Flux(0:N,0:N,1:2)
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
    subroutine calc_angular_momentum(Phi, LzPhi)
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N)
        complex(kind(0d0)),intent(out) :: LzPhi(0:N,0:N)
        integer                        :: i, j
        double precision               :: x, y
        complex(kind(0d0))             :: dPhi_dX(0:N,0:N), dPhi_dY(0:N,0:N)

        dPhi_dX(:,:) = dcmplx(0d0, 0d0)
        do j = 0, N
            do i = 0, N

                if (i > 1) then
                    dPhi_dX(i,j) = dPhi_dX(i,j) - 8d0*Phi(i-1,j)
                end if
                if (i > 2) then
                    dPhi_dX(i,j) = dPhi_dX(i,j) + Phi(i-2,j)
                end if
                if (i < N-1) then
                    dPhi_dX(i,j) = dPhi_dX(i,j) - Phi(i+2,j)
                end if
                if (i < N-2) then
                    dPhi_dX(i,j) = dPhi_dX(i,j) + 8d0*Phi(i+1,j)
                end if
            end do
        end do

        dPhi_dY(:,:) = dcmplx(0d0, 0d0)
        do j = 0, N
            do i = 0, N

                if (j > 1) then
                    dPhi_dY(i,j) = dPhi_dY(i,j) - 8d0*Phi(i,j-1)
                end if
                if (j > 2) then
                    dPhi_dY(i,j) = dPhi_dY(i,j) + Phi(i,j-2)
                end if
                if (j < N-1) then
                    dPhi_dY(i,j) = dPhi_dY(i,j) - Phi(i,j+2)
                end if
                if (j < N-2) then
                    dPhi_dY(i,j) = dPhi_dY(i,j) + 8d0*Phi(i,j+1)
                end if
            end do
        end do

        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i

                LzPhi(i,j) = -iu*(x*dPhi_dY(i,j) - y*dPhi_dX(i,j))/(12d0*dh)
            end do
        end do
    end subroutine

    ! Calculate Expected Angular Momentum Value (Excluding Plank constant)
    subroutine calc_angular_momentum_expected_value(Phi, LzPhi, Lz)
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N)
        complex(kind(0d0)),intent(out) :: LzPhi(0:N,0:N)
        double precision,intent(out)   :: Lz
        complex(kind(0d0))             :: temp
        temp = integrate(conjg(Phi)*LzPhi)
        ! Check wether sum is almost real here (Not implemented yet)
        Lz = dble(temp)
        if (aimag(temp) > 1d-6) then
            write (*, *) "Angular momentum is not properly calculated"
        end if
    end subroutine  
    
    subroutine calc_widths(Phi, sigma_x, sigma_y)
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