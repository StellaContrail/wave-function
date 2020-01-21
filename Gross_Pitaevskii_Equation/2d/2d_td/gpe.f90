! Time-development of Gross-Pitaevskii Equation in 2 dimensional space
! For a reason of Hamiltonian being not dense matrix, we will use just multiple dimensional array here

program main
    use constants
    use io
    use setting
    use mathf
    implicit none
    complex(kind(0d0)),allocatable :: Phi(:, :)   ! Wave function at next step
    double precision,allocatable   :: Pot(:, :)        ! Potential
    double precision,allocatable   :: Pot_TD(:, :)     ! Time-dependent Potential
    double precision               :: mu               ! chemical potential

    ! Coefficients and variables (not user defined)
    integer                        :: i                ! Loop variable
    double precision               :: prob
    double precision,allocatable   :: Flux(:,:,:)
    double precision,allocatable   :: Rot(:,:)
    double precision               :: width_x, width_y
    integer                        :: total_iterations ! How many iterations to calculate real time propagation
    complex(kind(0d0)),allocatable :: LzPhi(:,:)       ! Angular momentum operator on wave function
    double precision               :: Lz               ! Angular momentum itself
    double precision               :: t1, t2

    ! Output File Path
    character(*),parameter         :: fn_initial   = "data_initial.txt"
    character(*),parameter         :: fn_potential = "data_potential.txt"
    character(*),parameter         :: fn_result    = "data.txt"
    character(*),parameter         :: fn_flux      = "data_flux.txt"
    character(*),parameter         :: fn_rotation  = "data_rotation.txt"
    character(*),parameter         :: fn_phased    = "data_phased.txt"

    ! Allocation of variables
    allocate (Phi(0:N,0:N), Pot(0:N,0:N), Pot_TD(0:N,0:N))
    allocate (Flux(0:N,0:N,1:2), Rot(0:N,0:N), LzPhi(0:N,0:N))

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
    Pot_TD(:, :) = Pot(:, :)
    write (*, *) "- Loaded Wave Function <= ", "data_input.txt"
    write (*, *) "- Initialized wave functions and potential"

    ! Save initial wave function
    open(10, file=fn_initial)
    call output(10, Phi)
    close(10)
    write (*, *) "- Initial Wave Function => ", fn_initial
    
    ! Calculate z-component angular momentum expected value
    call calc_angular_momentum(Phi, LzPhi)
    call calc_angular_momentum_expected_value(Phi, LzPhi, Lz)
    write (*, '(X, A, F13.10, X, A)') "- Angular momentum <Lz> = ", Lz, "hbar"
    write (*, *) "- OMEGA = ", OMEGA

    print *, "Press Enter key to start calculation..."
    read (*, *)
    write (*, '(X, A)') "* Calculation has been initiated"

    open(10, file=fn_result)
    open(20, file=fn_potential)
    open(30, file=fn_flux)
    open(40, file=fn_rotation)
    !open(50, file="widths.txt")
    total_iterations = 1000
    call cpu_time(t1)
    do i = 1, total_iterations
        ! Evolve the system
        call evolve(Phi, Pot_TD, LzPhi)

        ! Save wave function, potential, flux, rotation, and widths every 50th step
        if (mod(i, 10) == 0) then
            ! Save wave function
            call output(10, Phi)

            ! Save potential form
            call output(20, Pot_TD)

            ! Save flux distribution
            Flux = calc_flux(Phi)
            call output(30, Flux)

            ! Save rotation of flux
            call calc_rotation(Flux, Rot)
            call output_rotation(40, Rot)
            
            ! Save width of the condensate
            !call calc_widths(Phi_next, width_x, width_y)
            !call output_widths(50, i, width_x, width_y)
        end if

        ! Refresh the iteration display every 500th step
        if (mod(i, 50) == 0) then
            if (i > 50) then
                write (*, '(A)', advance='no') char(13)
            end if
            write (*, '(X, A, I5, A, I5, A)', advance='no') "- ", i, "/", total_iterations, " completed"
        end if

        ! Vary potential form depending on time
        !call vary_potential(Pot, Pot_TD, i, total_iterations)
    end do
    call cpu_time(t2)
    write (*, *)
    ! Check if the probability is conserved
    !call integrate(abs(Phi)**2d0, prob)
    close (10)
    close (20)
    close (30)
    close (40)
    !close (50)
    write (*, '(X, A, F5.1, A)') "- Calculation took ", t2 - t1, " sec"

    call calc_angular_momentum(Phi, LzPhi)
    call calc_angular_momentum_expected_value(Phi, LzPhi, Lz)
    write (*, '(X, A, F13.10, A)') "Angular Momentum = ", Lz, " hbar"
    write (*, *) "- Result Wave Function => ", fn_result
end program 
