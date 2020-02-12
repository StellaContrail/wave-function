! Mathematical Procedures
! TODO: 境界条件を自由に決められるようにする 具体的にはHamiltoinian内及びLzPhi内の微分演算子の扱いを変える (Dark Solitonを得るため)
! REFER -> https://www.ias.ac.in/article/fulltext/pram/077/05/0929-0947

module mathf
    use constants
    implicit none

    ! Integrate given function with any precision
    interface integrate
        module procedure integrate_real, integrate_complex
    end interface
    ! Get phase of given complex number/position
    interface phase
        module procedure phase_complex
    end interface
contains
    ! Integration of real function
    function integrate_real(f) result(sum)
        double precision,intent(in)  :: f(0:N)
        double precision             :: sum
        integer                      :: i

        sum = 0d0
        do i = 0, N
            if (i == 0 .or. i == N) then
                sum = sum + 0.5d0*f(i)*dh
            else
                sum = sum + f(i)*dh
            end if
        end do
    end function
    ! Integration of complex function
    function integrate_complex(f) result(sum)
        complex(kind(0d0)),intent(in)  :: f(0:N)
        complex(kind(0d0))             :: sum
        integer                        :: i

        sum = 0d0
        do i = 0, N
            if (i == 0 .or. i == N) then
                sum = sum + 0.5d0*f(i)*dh
            else
                sum = sum + f(i)*dh
            end if
        end do
    end function

    ! Normalize wave function
    subroutine normalize(Phi)
        complex(kind(0d0)),intent(inout) :: Phi(0:N)
        Phi(:) = Phi(:) / sqrt(integrate(abs(Phi(:))**2d0))
    end subroutine normalize

    subroutine evolve(Phi, Pot, isImag, m)
        complex(kind(0d0)),intent(inout)     :: Phi(0:N)
        double precision,intent(in)          :: Pot(0:N)
        double precision,intent(in),optional :: m
        logical,intent(in)                   :: isImag
        complex(kind(0d0))                   :: Phi_old(0:N), Phi_new(0:N)
        complex(kind(0d0))                   :: temp(0:N), Atemp(0:N)
        double precision                     :: density(0:N)
        integer                              :: i, j
        complex(kind(0d0))                   :: lambda
        if (isImag .and. (.not.present(m))) then
            stop "parameter m is missing when calling evolve()"
        end if
        if (isImag) then
            lambda = dcmplx(1d0, 0d0)
        else
            lambda = iu
        end if

        Phi_old(:) = Phi(:)
        ! First term of Taylor expansion
        temp(:)    = Phi_old(:)
        Phi_new(:) = temp(:)
        ! Other terms of Taylor expansion
        if (isImag) then
            do i = 1, 15
                Atemp      = apply_hamiltonian(temp, abs(Phi_old)**2d0, Pot)
                temp(:)    = -lambda*Atemp(:)*dt_imag/(epsilon*i)
                Phi_new(:) = Phi_new(:) + temp(:)
            end do
        else
            do i = 1, 15
                Atemp      = apply_hamiltonian(temp, abs(Phi_old)**2d0, Pot)
                temp(:)    = -lambda*Atemp(:)*dt_real/(epsilon*i)
                Phi_new(:) = Phi_new(:) + temp(:)
            end do
        end if

        ! Mix the previous density and calculated wave function's density
        if (isImag) then
            density(:) = m*abs(Phi_old(:))**2d0 + (1d0-m)*abs(Phi_new(:))**2d0
        else
            density(:) = 0.5d0*abs(Phi_old(:))**2d0 + 0.5d0*abs(Phi_new(:))**2d0
        end if

        ! First term of Taylor expansion
        temp(:)    = Phi_old(:)
        Phi_new(:) = temp(:)
        ! Other terms of Taylor expansion
        if (isImag) then
            do i = 1, 15
                Atemp      = apply_hamiltonian(temp, density, Pot)
                temp(:)    = -lambda*Atemp(:)*dt_imag/(epsilon*i)
                Phi_new(:) = Phi_new(:) + temp(:)
            end do
        else
            do i = 1, 15
                Atemp      = apply_hamiltonian(temp, density, Pot)
                temp(:)    = -lambda*Atemp(:)*dt_real/(epsilon*i)
                Phi_new(:) = Phi_new(:) + temp(:)
            end do
        end if

        Phi(:) = Phi_new(:)
        if (isImag) then
            call normalize(Phi)
        end if
    end subroutine evolve

    ! Calculate HPhi (H:Hamiltonian, Phi:Wave function)
    function apply_hamiltonian(Phi, density, Pot) result(HPhi)
        complex(kind(0d0)),intent(in)  :: Phi(0:N)
        double precision,intent(in)    :: Pot(0:N)
        complex(kind(0d0))             :: HPhi(0:N)
        double precision,intent(in)    :: density(0:N)
        integer                        :: i
        HPhi(:) = 0d0
        ! Laplacian part (Five Point Stencil)
        do i = 0, N
            if (1 < i .and. i < N-1) then
                HPhi(i) = -Phi(i+2)+16d0*Phi(i+1)-30d0*Phi(i)+16d0*Phi(i-1)-Phi(i-2)
            else
                if (i == 0) then
                    HPhi(i) = -Phi(i+2)+16d0*Phi(i+1)-30d0*Phi(i)+16d0*Phi(i+1)-Phi(i+2)
                end if
                if (i == 1) then
                    HPhi(i) = -Phi(i+2)+16d0*Phi(i+1)-30d0*Phi(i)+16d0*Phi(i-1)-Phi(i+2)
                end if
                if (i == N) then
                    HPhi(i) = -Phi(i-2)+16d0*Phi(i-1)-30d0*Phi(i)+16d0*Phi(i-1)-Phi(i-2)
                end if
                if (i == N-1) then
                    HPhi(i) = -Phi(i-2)+16d0*Phi(i+1)-30d0*Phi(i)+16d0*Phi(i-1)-Phi(i-2)
                end if
            end if

        end do
        HPhi(:) = HPhi(:) / (12d0*dh*dh)
        
        do i = 0, N
            ! Kinetic energy
            HPhi(i) = -0.5d0*epsilon*epsilon*HPhi(i)

            ! External potential
            HPhi(i) = HPhi(i) + Pot(i)*Phi(i)

            ! Nonlinear part
            HPhi(i) = HPhi(i) + kappa*density(i)*Phi(i)
        end do
    end function

    subroutine apply_boundary_condition(Phi)
        complex(kind(0d0)),intent(inout) :: Phi(0:N)
        Phi(0) = 0.2d0
        Phi(N) = 0.2d0
    end subroutine

    ! Solve Energy Expected Value
    function solve_energy(Phi, Pot) result(mu)
        complex(kind(0d0)),intent(in) :: Phi(0:N)
        double precision,intent(in)   :: Pot(0:N)
        double precision              :: mu
        complex(kind(0d0))            :: HPhi(0:N)

        HPhi = apply_hamiltonian(Phi, abs(Phi)**2d0, Pot)
        mu = dble(integrate(conjg(Phi)*HPhi))
    end function

    ! Change phase
    subroutine shift_phase(Phi, m)
        complex(kind(0d0)),intent(inout)       :: Phi(0:N)
        integer,           intent(in)          :: m
        integer                                :: i, j
        double precision                       :: x

        do i = 0, N
            x = -xmax + dh*i
            if (x0_vortex < x) then
                Phi(i) = Phi(i) * exp(iu*m*pi)
            end if
        end do 
    end subroutine

    ! Get phase of the complex
    pure elemental double precision function phase_complex(Z)
        complex(kind(0d0)),intent(in) :: Z
        phase_complex = atan2(aimag(Z), real(Z))
    end function

    ! Calculate Probability Current
    function calc_flux(Phi) result(Flux)
        complex(kind(0d0)),intent(in) :: Phi(0:N)
        double precision              :: Flux(0:N)
        integer                       :: i, j

        Flux(0) = 0d0
        Flux(N) = 0d0

        do i = 1, N-1
            Flux(i) = aimag(conjg(Phi(i))*(Phi(i+1)-Phi(i-1))/(2d0*dh))
        end do

        !Flux(:,:,:) = Flux(:,:,:)*(hbar/mass)
    end function

    double precision function circulation(Phi)
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N)
        integer                       :: i
        double precision              :: sum
        integer                       :: lb, ub
        integer,parameter             :: width = 5
        lb = floor(N/2d0) - width
        ub = ceiling(N/2d0) + width

        sum = 0d0
        do i = lb, ub
            sum = sum + v_dr(Phi, i, lb, +1, 0)
            sum = sum + v_dr(Phi, ub, i, 0, +1)
            sum = sum + v_dr(Phi, i, ub, -1, 0)
            sum = sum + v_dr(Phi, lb, i, 0, -1)
        end do
        circulation = sum
    end function

    double precision function v_dr(Phi, i, j, coe_dx, coe_dy)
        complex(kind(0d0)),intent(in) :: Phi(0:N,0:N)
        integer,           intent(in) :: i, j, coe_dx, coe_dy
        v_dr = 0.5d0*(coe_dx*(phase(Phi(i+1, j)) - phase(Phi(i-1, j))) + coe_dy*(phase(Phi(i, j+1)) - phase(Phi(i, j-1))))
        ! v_dr = v_dr * (hbar/mass)
    end function
end module mathf
