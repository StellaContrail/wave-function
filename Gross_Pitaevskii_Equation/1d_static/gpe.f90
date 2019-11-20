program main
    use io
    use setting
    use mathf
    implicit none
    ! Mathematical constants
    double precision,parameter     :: pi   = acos(-1d0)       ! PI
    complex(kind(0d0)),parameter   :: iu   = dcmplx(0d0, 1d0) ! Imaginary unit
    ! Physical constants
    double precision,parameter     :: hbar = 1.05d-34         ! Reduced Plank constant
    ! Physical values
    integer                        :: N                ! Number of division in space
    integer                        :: DIM              ! Dimension of discretized wave function array
    complex(kind(0d0)),allocatable :: Phi_temp(:)      ! Temporary wave function used by zhbev
    complex(kind(0d0)),allocatable :: Phi_next(:)      ! Wave function at next step
    complex(kind(0d0)),allocatable :: Phi_prev(:)      ! Wave function at previous step
    double precision,allocatable   :: Pot(:)           ! Potential
    double precision,allocatable   :: H(:,:)
    double precision               :: dh               ! Step of distance in the x-direction
    double precision               :: xmax             ! largest x position (Boundary position)
    double precision               :: mass             ! mass of a boson
    double precision               :: omega            ! angular velocity of harmonic potential
    integer                        :: ParticleCount    ! number of bose particles
    double precision               :: ScatteringLength ! s-wave scattering length
    double precision               :: mu               ! chemical potential
    ! Coefficients and variables (not user defined)
    double precision               :: Azero            ! length of the harmonic oscillator ground state
    double precision               :: Xs               ! characteristic length of the condensate
    double precision               :: epsilon          ! squared ratio of Azero to Xs
    double precision               :: kappa            ! coeffient of the nonlinear term
    double precision               :: mu0              ! chemical potential returned by zhbev
    integer                        :: i                ! Loop variable
    logical                        :: loop_end_flag    ! Loop end flag
    character(len=1)               :: yn               ! Yes/No flag
    ! Definition of physical values (this could be replaced with I/O)
    ! These values are referenced from
    ! 'Numerical Solution of the Gross-Pitaevskii Equation for Bose-Einstein Condensation'
    ! by Weizhu Bao et al. (2003)
    mass             = 1.44d-25
    omega            = 20d0 * pi
    ParticleCount    = 100
    ScatteringLength = 5.1d-9
    N                = 2**8 - 1!129   ! N must be an odd number
    DIM              = N + 1 ! Include n=0 point
    allocate (Phi_next(0:N), Phi_prev(0:N), Pot(0:N))
    allocate (Phi_temp(0:N), H(0:N,0:N))
    ! Calculation of coefficients and variables using defined physical values
    xmax    = 10d0
    Azero   = sqrt(hbar/(omega*mass))
    Xs      = Azero   ! Usually chosen to be Azero for a weak/moderate interaction
    epsilon = (Azero/Xs)**2d0
    kappa   = (4d0*pi*ScatteringLength*ParticleCount/Azero)*(Azero/Xs)**5d0
    dh      = xmax / (n/2 + 0.5d0)

    ! Show configuration of fundamental physical constants
    print *, "Physical constants of the system----------------------------------"
    print *, "<Fundamental Physical Constants>"
    print *, "m  (Mass of the bose particle)   [kg] = ", mass
    print *, "omega (Angular velocity of HO)[rad/s] = ", omega
    print *, "N (Particle Count)            [count] = ", ParticleCount
    print *, "a (ScatteringLength)              [m] = ", ScatteringLength
    print *, "n (Dimension of the space)    [count] = ", N
    print *, "A0 (Length of the HO Ground State)[m] = ", Azero
    print *, "Xs (Characteristic Length)        [m] = ", Xs
    print *, "dh (Step of distance)             [m] = ", dh
    print *, "<Coefficients of NLSE terms>"
    print *, "Epsilon (A0/Xs)^2                     = ", epsilon
    print *, "Kappa (Coefficient of NL term)        = ", kappa
    print *, "------------------------------------------------------------------"
    write (*, *)
    
    print *, "Calculation Start-----------------------------------------"
    ! Initialization of wave function and potential
    call initialize(Phi_next, Phi_prev, Pot, N, dh, xmax, Azero)
    write (*, *) "- Initialized the wave function and potential function"

    ! Normalization of wave function
    call normalize(Phi_prev, N, dh)
    
    ! Output the initial wave function to file
    open(10, file="data_initial.txt")
    call output(10, Phi_prev, N, dh, xmax)
    close(10)
    write (*, *) "- Initial wave function has been saved into ", "data_initial.txt"
    print *, ""

    ! Start I/O Procedure
    open(10, file="data.txt")
    ! Set the loop exit flag to be false
    loop_end_flag = .false.
    ! Solve the inconsistent equation until the chemical potential converges
    do i = 1, 50
        write (*, '(A, I4, A)') "#", i, " step calculation began -------------------------------------------"
        ! Construct the hamiltonian using present wave function data
        call hamiltonian(H, Pot, abs(Phi_prev)**2d0, N, dh, epsilon, kappa)
        print *, "- New hamiltonian has been reconstructed with initial assumed density |"
        ! Solve Nonlinear Schroedinger Equation Using the Assumed Wave function
        call solve_eigen2(H, Phi_temp, mu0, DIM)
        print *, "- NLSE has been successfully calculated with initial assumed density  |"
        ! Take out the first (lowest energy) eigenvector
        Phi_next(:) = Phi_temp(:)
        ! Take an average between t and t+dt and make new probability whose time step is at t+0.5*dt approximately
        Phi_temp(:) = 0.5d0*(Phi_next(:) + Phi_prev(:))
        ! Construct the hamiltonian using present wave function data
        call hamiltonian(H, Pot, abs(Phi_temp)**2d0, N, dh, epsilon, kappa)
        print *, "- New hamiltonian has been reconstructed with averaged density        |"
        ! Solve Nonlinear Schroedinger Equation with averaged probability density calculated from t and t+dt wave functions
        call solve_eigen2(H, Phi_next, mu0, DIM)
        print *, "- NLSE has been successfully calculated with averaged density         |"

        ! Normalize the wave function
        call normalize(Phi_next, N, dh)
        print *, "- Wave function has been normalized                                   |"

        ! Check if chemical potential has been converged or not
        if (abs(mu0 - mu) < 1d-6) then
            print *, "* Chemical potential has been converged!                              |"
            loop_end_flag = .true.
        end if

        ! Substitute chemical potential
        mu = mu0
        print '(X, A, F9.5, A)', "- New Chemical potential : ", mu, " [J]                              |"
        print *, "- Chemical potential has been updated                                 |"

        Phi_prev = Phi_next
        print *, "- Finished                                                            |"

        if (loop_end_flag) then
            exit
        end if
    end do
    print *, "-----------------------------------------------------------------------"
    print *, ""
    write (*, *) "- All calculation procedures have been finished"
    call output(10, Phi_prev, N, dh, xmax)
    close (10)
    write (*, *) "- Calculation result has been saved into ", "data.txt"
    print *, "----------------------------------------------------------"
    write (*, *)
    print *, "Result of the calculation ----------------------------------------"
    print '(X, A, F9.5)', "mu (Chemical Potential) [J] = ", mu
    write (*, *)
    print *, "Wave function half of whose phase is changed by pi is saved into a file => ", "data_shifted.txt"
    call apply_phase_shift(Phi_prev(floor(N/2d0):N), N-floor(N/2d0)+1, iu, pi, Phi_prev(floor(N/2d0):N))
    open(10, file="data_shifted.txt")
    call output(10, Phi_prev, N, dh, xmax)
    close(10)
end program 
