! Imaginary time-development of Gross-Pitaevskii Equation in 2 dimensional space
! The Equations and constants are to be referenced to DOI:10.1016/S0021-9991(03)00102-5
! TODO: Store initial values from setting file into variables

program main
    use io
    use setting
    use mathf
    use constants
    implicit none

    ! Physical values
    integer                        :: N                ! Number of space steps in a direction
    complex(kind(0d0)),allocatable :: Phi_phased(:, :) ! Wave function phased by given angle
    complex(kind(0d0)),allocatable :: Phi_next(:, :)   ! Wave function at next step
    complex(kind(0d0)),allocatable :: Phi_prev(:, :)   ! Wave function at previous step
    complex(kind(0d0)),allocatable :: Pot(:, :)        ! Potential in Laboratory frame
    double precision               :: dh               ! Step of distance in the x-direction
    double precision               :: dt               ! Step of time     in the t-direction
    double precision               :: xmax             ! largest x position (Boundary position)
    double precision               :: mass             ! mass of a boson
    double precision               :: omega_x          ! angular velocity of harmonic potential in x direction
    double precision               :: omega_y          ! angular velocity of harmonic potential in y direction
    double precision               :: gamma            ! The ratio of angular velocity in y direction to the one in x direction
    integer                        :: ParticleCount    ! number of bose particles
    double precision               :: ScatteringLength ! s-wave scattering length (>0:Repulsive <0:Attractive)
    double precision               :: mu               ! chemical potential
    double precision,allocatable   :: j(:, :)          ! probability current
    double precision,allocatable   :: Phase_field(:,:) ! Phase field
    complex(kind(0d0)),allocatable :: LzPhi(:,:)       ! Angular momentum operator on wave function
    double precision               :: Lz               ! Angular momentum itself
    double precision,allocatable   :: Flux(:,:,:)      ! Probability current
    double precision               :: OMEGA            ! Angular velocity on Cranking model

    ! Coefficients and variables (not user defined)
    double precision               :: Azero            ! length of the harmonic oscillator ground state
    double precision               :: Xs               ! characteristic length of the condensate
    double precision               :: epsilon          ! squared ratio of Azero to Xs
    double precision               :: kappa            ! coeffient of the nonlinear term
    integer                        :: i                ! Loop variable
    double precision               :: mu_old           ! Chemical potential at previous step

    ! Output File Path
    character(*),parameter         :: fn_initial     = "data_initial.txt"         ! (Output) Trial wave function
    character(*),parameter         :: fn_potential   = "data_potential.txt"       ! (Output) Potential
    character(*),parameter         :: fn_result      = "data.txt"                 ! (Output) Result wave function
    character(*),parameter         :: fn_phased      = "data_phased.txt"          ! (Output) Phased wave function
    character(*),parameter         :: fn_phase_field = "data_phase_field.txt"     ! (Output) Phase distribution of wave function

    ! Definition of physical values (this could be replaced with I/O)
    ! These values are referenced from
    ! 'Numerical Solution of the Gross-Pitaevskii Equation for Bose-Einstein Condensation'
    ! by Weizhu Bao et al. (2003)
    mass             = 1.4d-25
    omega_x          = 20d0 * pi
    omega_y          = omega_x
    gamma            = omega_y / omega_x
    ParticleCount    = 1000
    ScatteringLength = 5.1d-9
    
    ! Number of steps in a direction
    N                = 100 - 1
    ! Angular velocity on Cranking model
    OMEGA            = 1000d0
    ! Variable allocation
    allocate (Phi_next(0:N,0:N), Phi_prev(0:N,0:N))
    allocate (Phi_phased(0:N, 0:N), LzPhi(0:N,0:N))
    allocate (Pot(0:N,0:N), j(0:N,0:N), Phase_field(0:N,0:N), Flux(0:N,0:N,1:2))

    ! Calculate/Set other variables
    xmax    = 15d0
    Azero   = sqrt(hbar/(omega_x*mass))
    Xs      = Azero
    epsilon = (Azero/Xs)**2d0
    kappa   = (4d0*pi*ScatteringLength*ParticleCount/Azero)*(Azero/Xs)**5d0
    dh      = xmax / (N/2 + 0.5d0)
    dt      = 0.01d0*dh*dh

    ! Display settings
    print *, "Physical constants of the system----------------------------------"
    print *, "<Fundamental Physical Constants>"
    print *, "m  (Mass of the bose particle)   [kg] = ", mass
    print *, "omegax(Angular velocity of HO)[rad/s] = ", omega_x
    print *, "omegay(Angular velocity of HO)[rad/s] = ", omega_y
    print *, "N (Particle Count)            [count] = ", ParticleCount
    print *, "a (ScatteringLength)              [m] = ", ScatteringLength
    print *, "n (Dimension of the space)    [count] = ", N + 1
    print *, "A0 (Length of the HO Ground State)[m] = ", Azero
    print *, "Xs (Characteristic Length)        [m] = ", Xs
    print *, "dh (Step of distance)                 = ", dh
    print *, "dt (Step of time)                     = ", dt
    print *, "<Coefficients of NLSE terms>"
    print *, "Epsilon (A0/Xs)^2                     = ", epsilon
    print *, "Kappa (Coefficient of NL term)        = ", kappa
    print *, "<Other Configuration Values>"
    print *, "Delta (4*pi*a*N/a_0)                  = ", (4d0*pi*ScatteringLength*ParticleCount)/Azero
    print *, "Healing length (8*pi*|a|*N/Xs^3)^-0.5 = ", ((8d0*pi*abs(ScatteringLength)*ParticleCount)/(Xs**3d0))**(-0.5d0)
    print *, "------------------------------------------------------------------"
    write (*, *)

    print *, "Press Enter to start calculation..."
    read (*, *)

    ! Initialization
    call initialize(Phi_next, Phi_prev, Pot, N, dh, xmax, gamma)
    call normalize(Phi_prev, N, dh)
    write (*, *) "- Initialized wave function and potential"

    ! Save initial wave function
    open(10, file=fn_initial)
    call output_complex(10, Phi_prev, N, dh, xmax)
    close(10)
    write (*, *) "- Initial Trial Wave Function => ", fn_initial

    ! Save potential form
    open(10, file=fn_potential)
    call output_potential(10, Pot, N, dh, xmax)
    close(10)
    write (*, *) "- Given Potential Form => ", fn_potential
    
    ! Save probability current
    call calc_flux(Phi_prev, N, mass, dh, Flux)
    open(11, file="data_initial_flux.txt")
    call output_flux(11, Flux, N, dh, xmax)
    close(11)
    print *, "- Probability current  => ", "data_initial_flux.txt"

    ! Calculate chemical potential of initial state
    call solve_energy(Phi_prev, Pot, LzPhi, N, epsilon, kappa, mu, dh, dt, OMEGA)

    ! Calculate z-component angular momentum expected value
    call calc_angular_momentum(Phi_prev, N, xmax, dh, LzPhi)
    call calc_angular_momentum_expected_value(Phi_prev, N, dh, LzPhi, Lz)
    write (*, '(X, A, F13.10, X, A)') "- Angular momentum <Lz> = ", Lz, "hbar"
    write (*, *) "- OMEGA = ", OMEGA

    ! Start solving 2D GPE
    do i = 1, 50000
        ! Calculate LzPhi
        call calc_angular_momentum(Phi_prev, N, xmax, dh, LzPhi)
        ! Evolve the system
        call evolve(Phi_prev, N, dt, dh, epsilon, kappa, abs(Phi_prev)**2d0, LzPhi, Pot, Phi_next, OMEGA)
        ! Mix the previous density and calculated wave function's density
        Phi_next(:,:) = sqrt((0.7d0*abs(Phi_prev(:,:))**2d0 + 0.3d0*abs(Phi_next(:,:))**2d0))
        call normalize(Phi_next, N, dh)

        ! Calculate chemical potential including cranking term
        mu_old = mu
        call solve_energy(Phi_next, Pot, LzPhi, N, epsilon, kappa, mu, dh, dt, OMEGA)

        ! Display how many iterations have been processed
        if (mod(i, 1000) == 0) then
            write (*, '(X, A, I0, A)') "* ", i, " calculations have been done"
        end if

        ! Substitution to calculate self-consistent equation again
        Phi_prev(:,:) = Phi_next(:,:)

        ! Check if chemical potential has been converged or not
        if (abs(mu_old - mu) < 1d-6) then
            print *, "- Chemical potential has been converged"
            write (*, '(X, A, I0, A)') "- Calculation successfully completed with ", i, " iterations"
            exit
        end if
    end do
    ! Warning of calculation divergence
    if (i >= 50000) then
        write (*, *) "* Calculation has been exceeded its iteration limit. Incorrect result is expected."
    end if
    write (*, *)

    ! Save calculation result
    open(10, file=fn_result)
    call output_complex(10, Phi_prev, N, dh, xmax)
    close (10)
    write (*, *) "- Result Wave Function => ", fn_result

    ! Make a vortex by shifting phase by pi continually
    call make_vortex(Phi_prev, N, xmax, dh, Phi_phased)
    open(11, file=fn_phased)
    call output_complex(11, Phi_phased, N, dh, xmax)
    close(11)
    print *, "- Phased Wave Function => ", fn_phased
    
    ! Save phase distribution of given wave function
    call get_phase_field(Phi_prev, N, Phase_field)
    open(11, file=fn_phase_field)
    call output_real(11, Phase_field, N, dh, xmax)
    close(11)
    print *, "- Phase field          => ", fn_phase_field

    ! Save probability current
    call calc_flux(Phi_prev, N, mass, dh, Flux)
    open(11, file="data_flux.txt")
    call output_flux(11, Flux, N, dh, xmax)
    close(11)
    print *, "- Probability current  => ", "data_flux.txt"

    ! Save angular momentum information
    open(10, file="angular_momentum.txt")
    ! Calculate angular momentum in z-direction
    call calc_angular_momentum(Phi_prev, N, xmax, dh, LzPhi)
    call output_complex(10, LzPhi, N, dh, xmax)
    close(10)
    ! Calculate z-component angular momentum expected value
    call calc_angular_momentum_expected_value(Phi_prev, N, dh, LzPhi, Lz)
    write (*, '(X, A, F13.10, X, A)') "- Angular momentum <Lz> = ", Lz, "hbar"

    ! [Temporary] External Potential + Self-consistent term
    open(10, file="total_potential.txt")
    call output_potential(10, Pot+kappa*abs(Phi_prev)**2d0, N, dh, xmax)
    close(10)

    ! Free allocated variables from memory
    deallocate (Phi_next, Phi_prev, Phi_phased, Pot, j, Phase_field, LzPhi, Flux)
end program 
