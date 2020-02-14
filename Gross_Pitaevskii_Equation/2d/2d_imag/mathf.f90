! Mathematical Procedures
module mathf
    use constants
    implicit none

    interface integrate
        module procedure integrate_real, integrate_complex
    end interface
    interface phase
        module procedure phase_coordinates, phase_complex
    end interface
contains
    ! Integration of real function
    function integrate_real(f) result(sum)
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
    end function
    ! Integration of complex function
    function integrate_complex(f) result(sum)
        complex(kind(0d0)),intent(in)  :: f(0:N, 0:N)
        complex(kind(0d0))             :: sum
        complex(kind(0d0))             :: sum_temp(0:N)
        integer                        :: i

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
    end function

    ! Normalize wave function
    subroutine normalize(Phi)
        complex(kind(0d0)),intent(inout) :: Phi(0:N, 0:N)
        Phi(:, :) = Phi(:, :) / sqrt(integrate(abs(Phi(:,:))**2d0))
    end subroutine normalize

    ! Imaginary-time propagation
    subroutine evolve(Phi, LzPhi, Pot)
        complex(kind(0d0)),intent(inout) :: Phi(0:N,0:N)
        complex(kind(0d0)),intent(in)    :: LzPhi(0:N, 0:N)
        double precision,intent(in)      :: Pot(0:N,0:N)
        complex(kind(0d0))               :: Phi_old(0:N,0:N), Phi_new(0:N,0:N)
        complex(kind(0d0))               :: temp(0:N,0:N), Atemp(0:N,0:N)
        double precision                 :: density(0:N,0:N)
        integer                          :: i

        Phi_old(:,:) = Phi(:,:)
        ! First term of Taylor expansion
        temp(:,:)    = Phi_old(:,:)
        Phi_new(:,:) = temp(:,:)
        ! Other terms of Taylor expansion
        do i = 1, 10
            Atemp = apply_hamiltonian(temp, abs(Phi_old)**2d0, LzPhi, Pot)
            temp(:,:)     = -Atemp(:,:)*dt/(epsilon*i)
            Phi_new(:,:) = Phi_new(:,:) + temp(:,:)
        end do

        ! Mix the previous density and calculated wave function's density
        density(:,:) = 0.7d0*abs(Phi_old(:,:))**2d0 + 0.3d0*abs(Phi_new(:,:))**2d0

        ! First term of Taylor expansion
        temp(:,:)    = Phi_old(:,:)
        Phi_new(:,:) = temp(:,:)
        ! Other terms of Taylor expansion
        do i = 1, 10
            Atemp = apply_hamiltonian(temp, density, LzPhi, Pot)
            temp(:,:)     = -Atemp(:,:)*dt/(epsilon*i)
            Phi_new(:,:) = Phi_new(:,:) + temp(:,:)
        end do

        Phi(:,:) = Phi_new(:,:)
        call normalize(Phi)
    end subroutine evolve

    ! Calculate HPhi (H:Hamiltonian, Phi:Wave function)
    function apply_hamiltonian(Phi, density, LzPhi, Pot) result(HPhi)
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N), LzPhi(0:N,0:N)
        double precision,intent(in)    :: Pot(0:N,0:N)
        complex(kind(0d0))             :: HPhi(0:N,0:N)
        double precision,intent(in)    :: density(0:N,0:N)
        integer                        :: i, j
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

                ! Cranking model
                HPhi(i,j) = HPhi(i,j) - OMEGA*LzPhi(i,j)/omega_x
            end do
        end do
    end function

    ! Solve Energy Expected Value
    function solve_energy(Phi, Pot, LzPhi) result(mu)
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N), LzPhi(0:N,0:N)
        double precision,intent(in)   :: Pot(0:N,0:N)
        double precision              :: mu
        complex(kind(0d0))            :: HPhi(0:N,0:N)

        HPhi = apply_hamiltonian(Phi, abs(Phi)**2d0, LzPhi, Pot)
        mu = dble(integrate(conjg(Phi)*HPhi))
    end function

    ! Make Quantized Vortex by changing the phase
    subroutine make_vortex(Phi)
        complex(kind(0d0)),intent(inout)  :: Phi(0:N,0:N)
        integer                           :: i, j
        double precision                  :: x, y

        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i
                
                Phi(i,j) = exp(iu*phase(y,x))*Phi(i,j)
            end do
        end do
    end subroutine

    ! Get phase of the complex
    pure elemental double precision function phase_complex(Z)
        complex(kind(0d0)),intent(in) :: Z
        phase_complex = atan2(aimag(Z), real(Z))
    end function

    ! Get phase information from given coordinates
    double precision function phase_coordinates(Y, X)
        double precision,intent(in) :: X, Y
        phase_coordinates = atan2(Y, X)
    end function

    ! Calculate LzPhi (Excluding Plank constant)
    function calc_angular_momentum(Phi) result(LzPhi)
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N)
        complex(kind(0d0))             :: LzPhi(0:N,0:N)
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
    end function

    ! Calculate Expected Angular Momentum Value (Excluding Plank constant)
    double precision function calc_angular_momentum_expected_value(Phi, LzPhi)
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N), LzPhi(0:N,0:N)
        integer                        :: i
        complex(kind(0d0))             :: sum_cmplx

        sum_cmplx = integrate(conjg(Phi)*LzPhi)

        ! Check if calculated Lz doesn't have imaginary part since the operator is Hermitian
        ! (<Lz> should be quantized, right? I don't have confidence so just leave it like this...)
        if (aimag(sum_cmplx) > 1d-6) then
            write (*, '(X, A, 2F13.10, A)') "Angular momentum is not properly calculated <Lz> = (", sum_cmplx, ")"
            stop
        end if

        ! ( If Lz should be quantized, check if it's integer and then store it into integer variable, Lz )
        calc_angular_momentum_expected_value = dble(sum_cmplx)
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
end module mathf
