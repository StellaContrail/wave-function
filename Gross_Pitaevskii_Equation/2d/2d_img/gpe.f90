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
    complex(kind(0d0)),allocatable :: Phi(:, :) ! Wave function phased by given angle
    complex(kind(0d0)),allocatable :: LzPhi(:,:)       ! Angular momentum operator on wave function
    double precision,allocatable   :: Pot(:, :)        ! Potential in Laboratory frame
    double precision               :: mu               ! chemical potential
    double precision               :: Lz               ! Angular momentum itself
    double precision               :: t1, t2           ! Calculation time variables

    ! Coefficients and variables (not user defined)
    integer                        :: i,j              ! Loop variable
    double precision               :: mu_old           ! Chemical potential at previous step
    double precision               :: x, y

    ! Output File Path
    character(*),parameter         :: fn_initial     = "data_initial.txt"         ! (Output) Trial wave function
    character(*),parameter         :: fn_potential   = "data_potential.txt"       ! (Output) Potential
    character(*),parameter         :: fn_result      = "data.txt"                 ! (Output) Result wave function
    character(*),parameter         :: fn_phased      = "data_phased.txt"          ! (Output) Phased wave function
    character(*),parameter         :: fn_phase_field = "data_phase_field.txt"     ! (Output) Phase distribution of wave function
    
    ! Variable allocation
    allocate (Phi(0:N,0:N), LzPhi(0:N,0:N), Pot(0:N,0:N))

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
    print *, "xmax (Absolute value of maximum x)[m] = ", xmax
    print *, "<Coefficients of NLSE terms>"
    print *, "Epsilon (A0/Xs)^2                     = ", epsilon
    print *, "Kappa (Coefficient of NL term)        = ", kappa
    print *, "<Other Configuration Values>"
    print *, "Delta (4*pi*a*N/a_0)                  = ", (4d0*pi*ScatteringLength*ParticleCount)/Azero
    print *, "Healing length (8*pi*|a|*N/Xs^3)^-0.5 = ", ((8d0*pi*abs(ScatteringLength)*ParticleCount)/(Xs**3d0))**(-0.5d0)
    print *, "------------------------------------------------------------------"
    write (*, *)

    ! Initialization
    call initialize(Phi, Pot)
    write (*, *) "- Initialized wave function and potential"

    ! Save initial wave function
    call output(fn_initial, Phi)
    write (*, *) "- Initial Trial Wave Function => ", fn_initial

    ! Save potential form
    call output(fn_potential, Pot)
    write (*, *) "- Given Potential Form => ", fn_potential
    
    ! Save probability current
    call output("data_initial_flux.txt", calc_flux(Phi))
    print *, "- Probability current  => ", "data_initial_flux.txt"

    ! Calculate chemical potential of initial state
    mu = solve_energy(Phi, Pot, LzPhi)

    ! Calculate z-component angular momentum expected value
    LzPhi = calc_angular_momentum(Phi)
    Lz    = calc_angular_momentum_expected_value(Phi, LzPhi)
    write (*, '(X, A, F13.10, X, A)') "- Angular momentum <Lz> = ", Lz, "hbar"
    write (*, *) "- OMEGA = ", OMEGA

    write (*, *)
    print *, "Press Enter to start calculation..."
    read (*, *)
    write (*, '(X, A)') "* Calculation has been initiated"

    ! Start solving 2D GPE
    call cpu_time(t1)
    do i = 1, 50000
        ! Calculate LzPhi
        LzPhi = calc_angular_momentum(Phi)
        ! Evolve the system
        call evolve(Phi, LzPhi, Pot)

        ! Calculate chemical potential including cranking term
        mu_old = mu
        mu = solve_energy(Phi, Pot, LzPhi)

        ! Display how many iterations have been processed
        if (mod(i, 1000) == 0) then
            write (*, '(X, A, I0, A)') "* ", i, " calculations have been done"
        end if

        ! Check if chemical potential has been converged or not
        if (abs(mu_old - mu) < 1d-6) then
            print *, "- Chemical potential has been converged"
            write (*, '(X, A, I0, A)') "- Calculation successfully completed with ", i, " iterations"
            exit
        end if
    end do
    call cpu_time(t2)
    ! Warning of calculation divergence
    if (i >= 50000) then
        write (*, *) "* Calculation has been exceeded its iteration limit. Incorrect result is expected."
    end if
    write (*, *)
    
    write (*, '(X, A, F5.1, A)') "- Calculation time took", t2 - t1, " seconds"
    
    

    ! Save calculation result
    call output(fn_result, Phi)
    write (*, *) "- Result Wave Function => ", fn_result
    
    ! Save phase distribution of given wave function
    call output(fn_phase_field, phase(Phi))
    print *, "- Phase field          => ", fn_phase_field

    ! Save probability current
    call output("data_flux.txt", calc_flux(Phi))
    print *, "- Probability current  => ", "data_flux.txt"

    ! Save angular momentum information
    LzPhi = calc_angular_momentum(Phi)
    call output("angular_momentum.txt", LzPhi)
    ! Calculate z-component angular momentum expected value
    Lz = calc_angular_momentum_expected_value(Phi, LzPhi)
    write (*, '(X, A, F13.10, X, A)') "- Angular momentum <Lz> = ", Lz, "hbar"

    ! [Temporary] External Potential + Self-consistent term
    call output("total_potential.txt", Pot+kappa*abs(Phi)**2d0)

    ! Make a vortex by shifting phase by pi continually
    call make_vortex(Phi)
    call output(fn_phased, Phi)
    print *, "- Phased Wave Function => ", fn_phased

    ! Free allocated variables from memory
    deallocate (Phi, Pot, LzPhi)
end program 
