! Mathematical Procedures
module mathf
    implicit none
contains
    ! Integration of function f using Trapezoidal rule
    ! f   : Integrand array
    ! N   : Dimension of space excluding the first element
    ! dh  : Step distance of space
    ! sum : The result of the integration
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

    ! Normalize the given function f
    ! f   : Function to be normalized
    ! N   : Dimension of f
    ! dh  : step distance of space
    subroutine normalize(f, N, dh)
        integer,intent(in)                :: N
        double precision,intent(in)       :: dh
        complex(kind(0d0)),intent(inout)  :: f(0:N, 0:N)
        double precision                  :: sum
        call integrate(abs(f(:, :))**2d0, N, dh, sum)
        f(:, :) = f(:, :) / sqrt(sum)
    end subroutine normalize

    ! Calculate C := exp(Ax)
    ! A : REAL array having dimension of NxN
    ! x : Complex array having dimension of N
    ! C : Complex array having dimension of N
    subroutine evolve(Phi_old, N, dt, dh, epsilon, kappa, density, Pot, Phi_next)
        integer,intent(in)             :: N
        double precision,intent(in)    :: dt, dh, epsilon, kappa, density(0:N,0:N), Pot(0:N,0:N)
        complex(kind(0d0)),intent(in)  :: Phi_old(0:N,0:N)
        complex(kind(0d0)),intent(out) :: Phi_next(0:N,0:N)
        integer                        :: i
        complex(kind(0d0))             :: temp(0:N,0:N), Atemp(0:N,0:N)

        ! First term of Taylor expansion
        temp(:,:)     = Phi_old(:,:)
        Phi_next(:,:) = temp(:,:)

        ! Other terms of Taylor expansion
        do i = 1, 4
            call apply_hamiltonian(temp, N, dh, epsilon, kappa, density, Pot, Atemp)
            ! Atemp = H*LastTerm
            temp(:,:)   = -Atemp(:,:)*dt/(epsilon*i)
            Phi_next(:,:) = Phi_next(:,:) + temp(:,:)
        end do
    end subroutine evolve

    ! Calculate C := HPhi
    ! Where H is Hamiltonian of the system
    ! and Phi is 2 dimensional wave function
    subroutine apply_hamiltonian(Phi, N, dh, epsilon, kappa, density, Pot, HPhi)
        integer,intent(in)             :: N
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N)
        complex(kind(0d0)),intent(out) :: HPhi(0:N,0:N)
        integer                        :: i, j
        double precision,intent(in)    :: dh, epsilon, kappa
        double precision,intent(in)    :: Pot(0:N,0:N), density(0:N,0:N)
        
        HPhi(:,:) = dcmplx(0d0, 0d0)
        ! Laplacian part
        do j = 0, N
            do i = 0, N
                HPhi(i,j) = -4d0*Phi(i,j)
                if (0 < i) then
                    HPhi(i,j) = HPhi(i,j) + Phi(i-1,j)
                end if
                if (0 < j) then
                    HPhi(i,j) = HPhi(i,j) + Phi(i,j-1)
                end if
                if (i < N) then
                    HPhi(i,j) = HPhi(i,j) + Phi(i+1,j)
                end if
                if (j < N) then
                    HPhi(i,j) = HPhi(i,j) + Phi(i,j+1)
                end if
            end do
        end do
        HPhi(:,:) = HPhi(:,:) / (dh*dh)

        ! Kinetic energy
        HPhi(:,:) = -0.5d0*epsilon*epsilon*HPhi(:,:)

        ! External potential
        HPhi(:,:) = HPhi(:,:) + Pot(:,:)

        ! Nonlinear part
        HPhi(:,:) = HPhi(:,:) + kappa*density(:,:)
    end subroutine

    ! Solve Energy Expected Value
    subroutine solve_energy(Phi, Pot, N, epsilon, kappa, mu, dh)
        integer,intent(in)             :: N
        double precision,intent(in)    :: dh, epsilon, kappa
        double precision,intent(in)    :: Pot(0:N,0:N)
        complex(kind(0d0)),intent(in)  :: Phi(0:N,0:N)
        double precision,intent(out)   :: mu
        complex(kind(0d0))             :: HPhi(0:N,0:N), PhiHPhi(0:N,0:N)
        double precision               :: f(0:N, 0:N)
        mu = 0d0

        call apply_hamiltonian(Phi, N, dh, epsilon, kappa, abs(Phi)**2d0, Pot, HPhi)
        PhiHPhi(:, :) = matmul(conjg(Phi(:,:)), HPhi(:,:))
        f(:, :) = dble(PhiHPhi(:, :))
        call integrate(f, N, dh, mu)
    end subroutine
end module mathf