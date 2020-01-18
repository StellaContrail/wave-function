! Mathematical Procedures
module mathf
    use constants
    implicit none
contains
    ! Integration of real function
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

    ! Normalize wave function
    subroutine normalize(Phi, N, dh)
        integer,intent(in)               :: N
        double precision,intent(in)      :: dh
        complex(kind(0d0)),intent(inout) :: Phi(0:N, 0:N)
        double precision                 :: sum
        call integrate(abs(Phi(:,:))**2d0, N, dh, sum)
        Phi(:, :) = Phi(:, :) / sqrt(sum)
    end subroutine normalize

    ! Imaginary-time propagation
    subroutine evolve(Phi_old, N, dt, dh, epsilon, kappa, density, Lz, Pot, Phi_next, OMEGA)
        integer,intent(in)             :: N
        double precision,intent(in)    :: dt, dh, epsilon, kappa, density(0:N,0:N), OMEGA
        complex(kind(0d0)),intent(in)  :: Phi_old(0:N,0:N), Lz(0:N, 0:N), Pot(0:N,0:N)
        complex(kind(0d0)),intent(out) :: Phi_next(0:N,0:N)
        complex(kind(0d0))             :: temp(0:N,0:N), Atemp(0:N,0:N)
        integer                        :: i
        ! First term of Taylor expansion
        temp(:,:)     = Phi_old(:,:)
        Phi_next(:,:) = temp(:,:)
        ! Other terms of Taylor expansion
        do i = 1, 10
            call apply_hamiltonian(temp, N, dh, dt, epsilon, kappa, density, Lz, Pot, Atemp, OMEGA)
            temp(:,:)     = -Atemp(:,:)*dt/(epsilon*i)
            Phi_next(:,:) = Phi_next(:,:) + temp(:,:)
        end do
    end subroutine evolve

    ! Calculate HPhi (H:Hamiltonian, Phi:Wave function)
    subroutine apply_hamiltonian(Phi, N, dh, dt, epsilon, kappa, density, LzPhi, Pot, HPhi, OMEGA)
        integer,intent(in)             :: N
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N), LzPhi(0:N,0:N), Pot(0:N,0:N)
        complex(kind(0d0)),intent(out) :: HPhi(0:N,0:N)
        double precision,intent(in)    :: dh, epsilon, kappa, dt
        double precision,intent(in)    :: OMEGA ! Rotation rate of the cranking model
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
                HPhi(i,j) = HPhi(i,j) - OMEGA*LzPhi(i,j)/(20d0*pi)
            end do
        end do
    end subroutine

    ! Solve Energy Expected Value
    subroutine solve_energy(Phi, Pot, Lz, N, epsilon, kappa, mu, dh, dt, OMEGA)
        integer,intent(in)            :: N
        double precision,intent(in)   :: dh, epsilon, kappa, dt, OMEGA
        complex(kind(0d0)),intent(in) :: Pot(0:N,0:N)
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N), Lz(0:N,0:N)
        double precision,intent(out)  :: mu
        complex(kind(0d0))            :: HPhi(0:N,0:N), sum_temp(0:N), sum_cmplx
        integer                       :: i
        mu = 0d0
        call apply_hamiltonian(Phi, N, dh, dt, epsilon, kappa, abs(Phi)**2d0, Lz, Pot, HPhi, OMEGA)
        sum_temp(:) = dcmplx(0d0, 0d0)
        do i = 0, N
            if (i == 0 .or. i == N) then
                sum_temp(:) = sum_temp(:) + 0.5d0*conjg(Phi(i,:))*HPhi(i,:)*dh
            else 
                sum_temp(:) = sum_temp(:) + conjg(Phi(i,:))*HPhi(i,:)*dh
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
        mu = dble(sum_cmplx)
    end subroutine

    ! Shift wave fuction's phase partially
    subroutine shift_phase(Phi_IN, N, x_start, x_end, y_start, y_end, Phi_OUT, angle)
        integer,intent(in)             :: x_start, x_end, y_start, y_end, N
        double precision,intent(in)    :: angle
        complex(kind(0d0)),intent(in)  :: Phi_IN(0:N,0:N)
        complex(kind(0d0)),intent(out) :: Phi_OUT(0:N,0:N)
        integer                        :: i
        Phi_OUT(:,:) = Phi_IN(:,:)
        do i = x_start, x_end
            Phi_OUT(i,y_start:y_end) = exp(iu*angle) * Phi_IN(i,y_start:y_end)
        end do
    end subroutine

    ! Make Quantized Vortex by changing the phase
    subroutine make_vortex(Phi, N, xmax, dh, Phi_phased)
        integer,intent(in)             :: N
        double precision,intent(in)    :: dh, xmax
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N)
        complex(kind(0d0)),intent(out) :: Phi_phased(0:N,0:N)
        integer                        :: i, j
        double precision               :: x, y, degree

        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i
                
                degree = atan2(y, x)
                Phi_phased(i,j) = exp(iu*degree)*Phi(i,j)
            end do
        end do
    end subroutine

    ! Get phase of the complex number
    function get_phase(Z) result(phase)
        complex(kind(0d0)),intent(in) :: Z
        double precision              :: phase
        
        phase = atan2(aimag(z), real(z))
    end function

    ! Get phase distribution of the field
    subroutine get_phase_field(Phi, N, Phase)
        integer,intent(in)            :: N
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N)
        double precision,intent(out)  :: Phase(0:N,0:N)
        integer                       :: i, j

        do j = 0, N
            do i = 0, N
                Phase(i,j) = get_phase(Phi(i,j))
            end do
        end do
    end subroutine

    ! Calculate LzPhi (Excluding Plank constant)
    subroutine calc_angular_momentum(Phi, N, xmax, dh, LzPhi)
        integer,intent(in)             :: N
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N)
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
    end subroutine

    ! Calculate Expected Angular Momentum Value (Excluding Plank constant)
    subroutine calc_angular_momentum_expected_value(Phi, N, dh, LzPhi, Lz)
        integer,intent(in)             :: N
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N)
        double precision,intent(in)    :: dh
        complex(kind(0d0)),intent(out) :: LzPhi(0:N,0:N)
        integer                        :: i
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

        ! Check if calculated Lz doesn't have imaginary part since the operator is Hermitian
        ! (<Lz> should be quantized, right? I don't have confidence so just leave it like this...)
        if (aimag(sum_cmplx) > 1d-6) then
            write (*, '(X, A, 2F13.10, A)') "Angular momentum is not properly calculated <Lz> = (", sum_cmplx, ")"
            stop
        end if

        ! ( If Lz should be quantized, check if it's integer and then store it into integer variable, Lz )
        Lz = dble(sum_cmplx)
    end subroutine  

    ! Calculate Probability Current
    subroutine calc_flux(Phi, N, mass, dh, Flux)
        integer,intent(in)            :: N
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N)
        double precision,intent(in)   :: dh, mass
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
end module mathf
