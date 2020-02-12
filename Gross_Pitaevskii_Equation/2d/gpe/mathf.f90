! Mathematical Procedures
! Maybe apply perfectly matched layer ? :https://arxiv.org/pdf/1312.1565.pdf

module mathf
    use constants
    implicit none

    ! Integrate given function with any precision
    interface integrate
        module procedure integrate_real, integrate_complex
    end interface
    ! Get phase of given complex number/position
    interface phase
        module procedure phase_coordinates, phase_complex
    end interface
    ! Get circulation of the system
    interface circulation
        module procedure circulation_flux, circulation_phase
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

    subroutine evolve(Phi, LzPhi, Pot, OMEGA, isImag, m)
        complex(kind(0d0)),intent(inout)     :: Phi(0:N,0:N)
        complex(kind(0d0)),intent(in)        :: LzPhi(0:N, 0:N)
        double precision,intent(in)          :: Pot(0:N,0:N), OMEGA
        double precision,intent(in),optional :: m
        logical,intent(in)                   :: isImag
        complex(kind(0d0))                   :: Phi_old(0:N,0:N), Phi_new(0:N,0:N)
        complex(kind(0d0))                   :: temp(0:N,0:N), Atemp(0:N,0:N)
        double precision                     :: density(0:N,0:N)
        integer                              :: i, j
        complex(kind(0d0))                   :: lambda
        double precision                     :: dt
        if (isImag .and. (.not.present(m))) then
            stop "parameter m is missing when calling evolve()"
        end if
        if (isImag) then
            lambda = 1d0
            dt = dt_imag
        else
            lambda = iu
            dt = dt_real
        end if
        Phi_old(:,:) = Phi(:,:)

        ! First term of Taylor expansion
        temp(:,:)    = Phi_old(:,:)
        Phi_new(:,:) = temp(:,:)
        ! Other terms of Taylor expansion
        do i = 1, 15
            Atemp        = apply_hamiltonian(temp, abs(Phi_old)**2d0, LzPhi, Pot, OMEGA)
            temp(:,:)    = -lambda*Atemp(:,:)*dt/(epsilon*i)
            Phi_new(:,:) = Phi_new(:,:) + temp(:,:)
        end do

        ! Mix the previous density and calculated wave function's density
        if (isImag) then
            density(:,:) = m*abs(Phi_old(:,:))**2d0 + (1d0-m)*abs(Phi_new(:,:))**2d0
        else
            density(:,:) = 0.5d0*abs(Phi_old(:,:))**2d0 + 0.5d0*abs(Phi_new(:,:))**2d0
        end if

        ! First term of Taylor expansion
        temp(:,:)    = Phi_old(:,:)
        Phi_new(:,:) = temp(:,:)
        ! Other terms of Taylor expansion
        do i = 1, 15
            Atemp        = apply_hamiltonian(temp, density, LzPhi, Pot, OMEGA)
            temp(:,:)    = -lambda*Atemp(:,:)*dt/(epsilon*i)
            Phi_new(:,:) = Phi_new(:,:) + temp(:,:)
        end do

        Phi(:,:) = Phi_new(:,:)
        if (isImag) then
            call normalize(Phi)
        end if
    end subroutine evolve

    ! Calculate HPhi (H:Hamiltonian, Phi:Wave function)
    function apply_hamiltonian(Phi, density, LzPhi, Pot, OMEGA) result(HPhi)
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N), LzPhi(0:N,0:N)
        double precision,intent(in)    :: Pot(0:N,0:N)
        complex(kind(0d0))             :: HPhi(0:N,0:N)
        double precision,intent(in)    :: density(0:N,0:N), OMEGA
        integer                        :: i, j
        HPhi(:,:) = 0d0
        ! Laplacian part (Five Point Stencil)
        do j = 0, N
            do i = 0, N
                if (i == 0) then
                    HPhi(i,j) = -Phi(i+2,j)+16d0*Phi(i+1,j)-30d0*Phi(i,j)-16d0*Phi(i+1,j)+Phi(i+2,j)
                end if
                if (i == 1) then
                    HPhi(i,j) = -Phi(i+2,j)+16d0*Phi(i+1,j)-30d0*Phi(i,j)+16d0*Phi(i-1,j)+Phi(i+2,j)
                end if
                if (i == N-1) then
                    HPhi(i,j) = Phi(i-2,j)+16d0*Phi(i+1,j)-30d0*Phi(i,j)+16d0*Phi(i-1,j)-Phi(i-2,j)
                end if
                if (i == N) then
                    HPhi(i,j) = Phi(i-2,j)-16d0*Phi(i-1,j)-30d0*Phi(i,j)+16d0*Phi(i-1,j)-Phi(i-2,j)
                end if
                if (j == 0) then
                    HPhi(i,j) = HPhi(i,j) -Phi(i,j+2)+16d0*Phi(i,j+1)-30d0*Phi(i,j)-16d0*Phi(i,j+1)+Phi(i,j+2)
                end if
                if (j == 1) then
                    HPhi(i,j) = HPhi(i,j) -Phi(i,j+2)+16d0*Phi(i,j+1)-30d0*Phi(i,j)+16d0*Phi(i,j-1)+Phi(i,j+2)
                end if
                if (j == N) then
                    HPhi(i,j) = HPhi(i,j) +Phi(i,j-2)-16d0*Phi(i,j-1)-30d0*Phi(i,j)+16d0*Phi(i,j-1)-Phi(i,j-2)
                end if
                if (j == N-1) then
                    HPhi(i,j) = HPhi(i,j) +Phi(i,j-2)+16d0*Phi(i,j+1)-30d0*Phi(i,j)+16d0*Phi(i,j-1)-Phi(i,j-2)
                end if

                if ( (1 < i .and. i < N-1) .and. (1 < j .and. j < N-1) ) then
                    HPhi(i,j) = -Phi(i+2,j)+16d0*Phi(i+1,j)-30d0*Phi(i,j)+16d0*Phi(i-1,j)-Phi(i-2,j)
                    HPhi(i,j) = HPhi(i,j) -Phi(i,j+2)+16d0*Phi(i,j+1)-30d0*Phi(i,j)+16d0*Phi(i,j-1)-Phi(i,j-2)
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
                HPhi(i,j) = HPhi(i,j) - epsilon*(OMEGA/omega_x)*LzPhi(i,j)
            end do
        end do
    end function

    ! Solve Energy Expected Value
    function solve_energy(Phi, Pot, LzPhi, OMEGA) result(mu)
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N), LzPhi(0:N,0:N)
        double precision,intent(in)   :: Pot(0:N,0:N), OMEGA
        double precision              :: mu
        complex(kind(0d0))            :: HPhi(0:N,0:N)

        HPhi = apply_hamiltonian(Phi, abs(Phi)**2d0, LzPhi, Pot, OMEGA)
        mu = dble(integrate(conjg(Phi)*HPhi))
    end function

    ! Make Quantized Vortex by changing the phase
    subroutine make_vortex(Phi, m, radius)
        complex(kind(0d0)),intent(inout)       :: Phi(0:N,0:N)
        integer,           intent(in)          :: m
        double precision,  intent(in),optional :: radius
        integer                                :: i, j
        double precision                       :: x, y

        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i
                if (present(radius)) then
                    if ((x-x0_vortex)**2d0 + (y-y0_vortex)**2d0 < radius**2d0) then
                        Phi(i,j) = Phi(i,j) * exp(iu*m*phase(y-y0_vortex, x-x0_vortex))
                    end if
                else
                    Phi(i,j) = Phi(i,j) * exp(iu*m*phase(y-y0_vortex, x-x0_vortex))
                end if
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
    function calc_LzPhi(Phi) result(LzPhi)
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N)
        complex(kind(0d0))             :: LzPhi(0:N,0:N)
        integer                        :: i, j
        double precision               :: x, y
        complex(kind(0d0))             :: dPhi_dX(0:N,0:N), dPhi_dY(0:N,0:N)

        dPhi_dX(:,:) = dcmplx(0d0, 0d0)
        do j = 0, N
            do i = 0, N

                if (i == 0) then
                    dPhi_dX(i,j) = -Phi(i+2,j)+8d0*Phi(i+1,j)+8d0*Phi(i+1,j)-Phi(i+2,j)
                end if
                if (i == 1) then
                    dPhi_dX(i,j) = -Phi(i+2,j)+8d0*Phi(i+1,j)-8d0*Phi(i-1,j)-Phi(i+2,j)
                end if
                if (i == N) then
                    dPhi_dX(i,j) = Phi(i-2,j)-8d0*Phi(i-1,j)-8d0*Phi(i-1,j)+Phi(i-2,j)
                end if
                if (i == N-1) then
                    dPhi_dX(i,j) = Phi(i-2,j)+8d0*Phi(i+1,j)-8d0*Phi(i-1,j)+Phi(i-2,j)
                end if

                if ( (1 < i .and. i < N-1) .and. (1 < j .and. j < N-1) ) then
                    dPhi_dX(i,j) = -Phi(i+2,j)+8d0*Phi(i+1,j)-8d0*Phi(i-1,j)+Phi(i-2,j)
                end if
            end do
        end do

        dPhi_dY(:,:) = dcmplx(0d0, 0d0)
        do j = 0, N
            do i = 0, N
                if (j == 0) then
                    dPhi_dY(i,j) = -Phi(i,j+2)+8d0*Phi(i,j+1)+8d0*Phi(i,j+1)-Phi(i,j+2)
                end if
                if (j == 1) then
                    dPhi_dY(i,j) = -Phi(i,j+2)+8d0*Phi(i,j+1)-8d0*Phi(i,j-1)-Phi(i,j+2)
                end if
                if (j == N) then
                    dPhi_dY(i,j) = Phi(i,j-2)-8d0*Phi(i,j-1)-8d0*Phi(i,j-1)+Phi(i,j-2)
                end if
                if (j == N-1) then
                    dPhi_dY(i,j) = Phi(i,j-2)+8d0*Phi(i,j+1)-8d0*Phi(i,j-1)+Phi(i,j-2)
                end if

                if ( (1 < i .and. i < N-1) .and. (1 < j .and. j < N-1) ) then
                    dPhi_dY(i,j) = -Phi(i,j+2)+8d0*Phi(i,j+1)-8d0*Phi(i,j-1)+Phi(i,j-2)
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
    ! flag : Handle errors?
    double precision function calc_Lz(Phi, LzPhi, flag)
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N), LzPhi(0:N,0:N)
        logical,intent(in)             :: flag
        complex(kind(0d0))             :: sum_cmplx

        sum_cmplx = integrate(conjg(Phi)*LzPhi)

        if (flag) then
            if (aimag(sum_cmplx) > 1d-6) then
                write (*, '(X, A, F0.10, X, F0.10, A)') "Angular momentum is not properly calculated <Lz> = (", sum_cmplx, ")"
                stop
            end if
        end if

        calc_Lz = dble(sum_cmplx)
    end function  

    ! Calculate Probability Current
    function calc_flux(Phi) result(Flux)
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N)
        double precision              :: Flux(0:N,0:N,1:2)
        integer                       :: i, j
        do j = 0, N
            do i = 0, N
                if (i == 0) then
                    Flux(i,j,1) = aimag(conjg(Phi(i,j))*(-Phi(i+2,j)+8d0*Phi(i+1,j)+8d0*Phi(i+1,j)-Phi(i+2,j))/(12d0*dh))
                end if
                if (j == 0) then
                    Flux(i,j,2) = aimag(conjg(Phi(i,j))*(-Phi(i,j+2)+8d0*Phi(i,j+1)+8d0*Phi(i,j+1)-Phi(i,j+2))/(12d0*dh))
                end if
                if (i == 1) then
                    Flux(i,j,1) = aimag(conjg(Phi(i,j))*(-Phi(i+2,j)+8d0*Phi(i+1,j)-8d0*Phi(i-1,j)-Phi(i+2,j))/(12d0*dh))
                end if
                if (j == 1) then
                    Flux(i,j,2) = aimag(conjg(Phi(i,j))*(-Phi(i,j+2)+8d0*Phi(i,j+1)-8d0*Phi(i,j-1)-Phi(i,j+2))/(12d0*dh))
                end if

                if (i == N-1) then
                    Flux(i,j,1) = aimag(conjg(Phi(i,j))*(Phi(i-2,j)+8d0*Phi(i+1,j)-8d0*Phi(i-1,j)+Phi(i-2,j))/(12d0*dh))
                end if
                if (j == N-1) then
                    Flux(i,j,2) = aimag(conjg(Phi(i,j))*(Phi(i,j-2)+8d0*Phi(i,j+1)-8d0*Phi(i,j-1)+Phi(i,j-2))/(12d0*dh))
                end if
                if (i == N) then
                    Flux(i,j,1) = aimag(conjg(Phi(i,j))*(Phi(i-2,j)-8d0*Phi(i-1,j)-8d0*Phi(i-1,j)+Phi(i-2,j))/(12d0*dh))
                end if
                if (j == N) then
                    Flux(i,j,2) = aimag(conjg(Phi(i,j))*(Phi(i,j-2)-8d0*Phi(i,j-1)-8d0*Phi(i,j-1)+Phi(i,j-2))/(12d0*dh))
                end if

                if ( (1 < i .and. i < N-1) .and. (1 < j .and. j < N-1) ) then
                    Flux(i,j,1) = aimag(conjg(Phi(i,j))*(-Phi(i+2,j)+8d0*Phi(i+1,j)-8d0*Phi(i-1,j)+Phi(i-2,j))/(12d0*dh))
                    Flux(i,j,2) = aimag(conjg(Phi(i,j))*(-Phi(i,j+2)+8d0*Phi(i,j+1)-8d0*Phi(i,j-1)+Phi(i,j-2))/(12d0*dh))
                end if
            end do
        end do

        !Flux(:,:,:) = Flux(:,:,:)*(hbar/mass)
    end function

    ! Calculate rotation of flux vectors
    subroutine calc_rotation(Flux, Rot)
        double precision,intent(in)  :: Flux(0:N,0:N,1:2)
        double precision,intent(out) :: Rot(0:N,0:N)
        double precision             :: dFluxY_dX(0:N,0:N), dFluxX_dY(0:N,0:N)
        integer                      :: i, j
        double precision             :: x, y

        do j = 0, N
            y = -xmax + dh*j
            do i = 0, N
                x = -xmax + dh*i

                if (i == 0) then
                    dFluxY_dX(i,j) = -Flux(i+2,j,2)+8d0*Flux(i+1,j,2)+8d0*Flux(i+1,j,2)-Flux(i+2,j,2)
                end if
                if (i == 1) then
                    dFluxY_dX(i,j) = -Flux(i+2,j,2)+8d0*Flux(i+1,j,2)-8d0*Flux(i-1,j,2)-Flux(i+2,j,2)
                end if
                if (i == N) then
                    dFluxY_dX(i,j) = Flux(i-2,j,2)-8d0*Flux(i-1,j,2)-8d0*Flux(i-1,j,2)+Flux(i-2,j,2)
                end if
                if (i == N-1) then
                    dFluxY_dX(i,j) = Flux(i-2,j,2)+8d0*Flux(i+1,j,2)-8d0*Flux(i-1,j,2)+Flux(i-2,j,2)
                end if

                if (j == 0) then
                    dFluxX_dY(i,j) = -Flux(i,j+2,1)+8d0*Flux(i,j+1,1)+8d0*Flux(i,j+1,1)-Flux(i,j+2,1)
                end if
                if (j == 1) then
                    dFluxX_dY(i,j) = -Flux(i,j+2,1)+8d0*Flux(i,j+1,1)-8d0*Flux(i,j-1,1)-Flux(i,j+2,1)
                end if
                if (j == N) then
                    dFluxX_dY(i,j) = Flux(i,j-2,1)-8d0*Flux(i,j-1,1)-8d0*Flux(i,j-1,1)+Flux(i,j-2,1)
                end if
                if (j == N-1) then
                    dFluxX_dY(i,j) = Flux(i,j-2,1)+8d0*Flux(i,j+1,1)-8d0*Flux(i,j-1,1)+Flux(i,j-2,1)
                end if

                if ( (1 < i .and. i < N-1) .and. (1 < j .and. j < N-1) ) then
                    dFluxY_dX(i,j) = -Flux(i+2,j,2)+8d0*Flux(i+1,j,2)-8d0*Flux(i-1,j,2)+Flux(i-2,j,2)
                    dFluxX_dY(i,j) = -Flux(i,j+2,1)+8d0*Flux(i,j+1,1)-8d0*Flux(i,j-1,1)+Flux(i,j-2,1)
                end if

                dFluxY_dX = dFluxY_dX / ( 12d0*dh )
                dFluxX_dY = dFluxX_dY / ( 12d0*dh )

                Rot(i,j) = dFluxY_dX(i,j) - dFluxX_dY(i,j)
            end do
        end do
    end subroutine

    double precision function circulation_flux(Phi, Flux)
        complex(kind(0d0)),intent(in) :: Phi(0:N, 0:N)
        double precision,  intent(in) :: Flux(0:N,0:N,1:2)
        integer                       :: i
        double precision              :: sum
        integer                       :: lb, ub
        integer,parameter             :: width = 10
        lb = floor(N/2d0) - width
        ub = ceiling(N/2d0) + width

        sum = 0d0
        ! C1
        do i = lb, ub
            sum = sum + v_dr_(Phi, Flux, i, lb, +1, 0)
        end do
        ! C2
        do i = lb, ub
            sum = sum + v_dr_(Phi, Flux, ub, i, 0, +1)
        end do
        ! C3
        do i = ub, lb, -1
            sum = sum + v_dr_(Phi, Flux, i, ub, -1, 0)
        end do
        ! C4
        do i = ub, lb, -1
            sum = sum + v_dr_(Phi, Flux, lb, i, 0, -1)
        end do
        circulation_flux = sum
    end function
    
    double precision function v_dr_(Phi, Flux, i, j, coe_dx, coe_dy)
        complex(kind(0d0)), intent(in) :: Phi(0:N, 0:N)
        double precision,   intent(in) :: Flux(0:N,0:N,1:2)
        integer,            intent(in) :: i, j, coe_dx, coe_dy

        v_dr_ = coe_dx*Flux(i,j,1)/abs(Phi(i,j))**2d0 + coe_dy*Flux(i,j,2)/abs(Phi(i,j))**2d0
        v_dr_ = v_dr_ * dh
    end function

    double precision function circulation_phase(Phi)
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N)
        integer                       :: i
        double precision              :: sum
        integer                       :: lb, ub
        integer,parameter             :: width = 10
        lb = floor(N/2d0) - width
        ub = ceiling(N/2d0) + width

        sum = 0d0
        ! C1
        do i = lb, ub
            if(phase(Phi(i+1,lb)) - phase(Phi(i-1,lb)) < 0) then
                sum = sum + 0.5d0*( (phase(Phi(i+1,lb))+2d0*pi) - phase(Phi(i-1,lb)) )
            else
                sum = sum + v_dr(Phi, i, lb, +1, 0)
            end if
        end do
        ! C2
        do i = lb, ub
            if(phase(Phi(ub,i+1)) - phase(Phi(ub,i-1)) < 0) then
                sum = sum + 0.5d0*( (phase(Phi(ub,i+1))+2d0*pi) - phase(Phi(ub,i-1)) )
            else
                sum = sum + v_dr(Phi, ub, i, 0, +1)
            end if
        end do
        ! C3
        do i = ub, lb, -1
            if(phase(Phi(i-1,ub)) - phase(Phi(i+1,ub)) < 0) then
                sum = sum - 0.5d0*( (phase(Phi(i+1,ub))) - (phase(Phi(i-1,ub))+2d0*pi) )
            else
                sum = sum + v_dr(Phi, i, ub, -1, 0)
            end if
        end do
        ! C4
        do i = ub, lb, -1
            if (phase(Phi(lb,i-1)) - phase(Phi(lb,i+1)) < 0) then
                sum = sum - 0.5d0*( (phase(Phi(lb, i+1))) - (phase(Phi(lb, i-1))+2d0*pi) )
            else
                sum = sum + v_dr(Phi, lb, i, 0, -1)
            end if
        end do
        circulation_phase = sum
    end function

    double precision function v_dr(Phi, i, j, coe_dx, coe_dy)
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N)
        integer,           intent(in) :: i, j, coe_dx, coe_dy

        v_dr = coe_dx*(phase(Phi(i+1, j))-phase(Phi(i-1,j)))
        v_dr = v_dr + coe_dy*(phase(Phi(i, j+1))-phase(Phi(i, j-1)))
        v_dr = 0.5d0 * v_dr
    end function
end module mathf
