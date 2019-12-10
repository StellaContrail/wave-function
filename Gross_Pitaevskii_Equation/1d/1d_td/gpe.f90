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
    complex(kind(0d0)),allocatable :: Phi_temp(:)      ! Temporary wave function used by zhbev
    complex(kind(0d0)),allocatable :: Phi_next(:)      ! Wave function at next step
    complex(kind(0d0)),allocatable :: Phi_prev(:)      ! Wave function at previous step
    double precision,allocatable   :: Pot(:)           ! Potential
    double precision,allocatable   :: Pot_TD(:)        ! Time-Dependant Potential
    double precision,allocatable   :: H(:,:)
    double precision               :: dh               ! Step of distance in the x-direction
    double precision               :: dt               ! Step of time     in the t-direction
    double precision               :: xmax             ! largest x position (Boundary position)
    double precision               :: mass             ! mass of a boson
    double precision               :: omega            ! angular velocity of harmonic potential
    integer                        :: ParticleCount    ! number of bose particles
    double precision               :: ScatteringLength ! s-wave scattering length
    double precision               :: mu               ! chemical potential
    double precision,allocatable   :: j(:)             ! probability current
    ! Coefficients and variables (not user defined)
    double precision               :: Azero            ! length of the harmonic oscillator ground state
    double precision               :: Xs               ! characteristic length of the condensate
    double precision               :: epsilon          ! squared ratio of Azero to Xs
    double precision               :: kappa            ! coeffient of the nonlinear term
    double precision,allocatable   :: mus(:)           ! multiple chemical potentials returned by zhbev
    integer                        :: i                ! Loop variable
    double precision               :: prob             ! Probability of the system (should be 1)
    character(:),allocatable       :: string           ! ouput string
    logical                        :: enable           ! enable output
    integer                        :: iter_interval    ! output every iter_interval
    integer                        :: output_count     ! Number of processed data
    ! Definition of physical values (this could be replaced with I/O)
    ! These values are referenced from
    ! 'Numerical Solution of the Gross-Pitaevskii Equation for Bose-Einstein Condensation'
    ! by Weizhu Bao et al. (2003)
    mass             = 1.4d-25
    omega            = 20d0 * pi
    ParticleCount    = 100
    ScatteringLength = 5.1d-9
    N                = 2**8 - 1
    allocate (Phi_next(0:N), Phi_prev(0:N), Pot(0:N), mus(0:N), j(0:N), Pot_TD(0:N))
    allocate (Phi_temp(0:N), H(0:N,0:N))
    ! Calculation of coefficients and variables using defined physical values
    xmax    = 10d0
    Azero   = sqrt(hbar/(omega*mass))
    Xs      = Azero   ! Usually chosen to be Azero for a weak/moderate interaction
    epsilon = (Azero/Xs)**2d0
    kappa   = (4d0*pi*ScatteringLength*ParticleCount/Azero)*(Azero/Xs)**5d0
    dh      = xmax / (n/2 + 0.5d0)
    dt      = 0.4d0*dh*dh

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
    print *, "dt (Step of time)                 [s] = ", dt
    print *, "<Coefficients of NLSE terms>"
    print *, "Epsilon (A0/Xs)^2                     = ", epsilon
    print *, "Kappa (Coefficient of NL term)        = ", kappa
    print *, "<Other Configuration Values>"
    print *, "Delta (4*pi*a*N/a_0)                  = ", (4d0*pi*ScatteringLength*ParticleCount)/Azero
    print *, "Healing length (8*pi*|a|*N/Xs^3)^-0.5 = ", ((8d0*pi*abs(ScatteringLength)*ParticleCount)/(Xs**3d0))**(-0.5d0)
    print *, "------------------------------------------------------------------"
    write (*, *)
    print *, "Press any key to start calculation..."
    read (*, *)
    print *, "Calculation Start-----------------------------------------"
    ! Initialization of wave function and potential
    call initialize(Phi_next, Phi_prev, Pot, N, dh, xmax)
    write (*, *) "- Initialized the wave function and potential function"
    Pot_TD(:) = Pot(:)

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
    open(11, file="data_current.txt")
    open(30, file="data_potential.txt")
    allocate(character(len=80) :: string)
    enable = .true.
    iter_interval = 50
    output_count = 0
    ! Solve the inconsistent equation until the wave function converges
    do i = 1, 50000
        write (string, '(A, I6, A)') "#", i, " step calculation began ---------------------------------------"
        call print_ex(string, enable, 'E', iter_interval, i)
        ! Construct the hamiltonian using present wave function data
        call hamiltonian(H, Pot_TD, abs(Phi_prev)**2d0, N, dh, epsilon, kappa)
        write (string, '(X, A)') "- New hamiltonian has been reconstructed with first assumed density |"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Solve Nonlinear Schroedinger Equation Using the Assumed Wave function
        call exp_mat(H, Phi_prev, N, dt, epsilon, iu, Phi_next)
        write (string, '(X, A)') "- NLSE has been successfully calculated with first assumed density  |"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Calculate the average wave function which are located at two different time points
        Phi_temp(:) = sqrt(0.5d0*(abs(Phi_next(:))**2d0 + abs(Phi_prev(:))**2d0))

        ! Construct the hamiltonian using present wave function data
        call hamiltonian(H, Pot_TD, abs(Phi_temp)**2d0, N, dh, epsilon, kappa)
        write (string, '(X, A)') "- New hamiltonian has been reconstructed with averaged density      |"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Solve Nonlinear Schroedinger Equation Using the Assumed Wave function
        call exp_mat(H, Phi_prev, N, dt, epsilon, iu, Phi_next)
        write (string, '(X, A)') "- NLSE has been successfully calculated with averaged density       |"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Check if the probability is conserved
        call integrate(abs(Phi_next)**2d0, N, dh, prob)
        write (string, '(X, A, F15.10, A)') "- Probability       = ", prob, "                               |"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Calculate chemical potential
        call expected_value_symm(Phi_next, H, N, mu, dh)
        write (string, '(X, A, F10.5, A)') "- Chemical Potential = ", mu, "                                   |"
        call print_ex(string, enable, 'E', iter_interval, i)

        if (mod(i, 50) == 0) then
            call output(10, Phi_next, N, dh, xmax)
            write (string, '(X, A)') "- Wave function has been saved into file                            |"
            call print_ex(string, enable, 'E', iter_interval, i)
            call output_real(30, Pot_TD, N, dh, xmax)
            output_count = output_count + 1

            ! Probability current calculation and output
            call calc_current(Phi_next, N, dh, hbar, mass, j)
            call output_current(11, j, N, dh, xmax)
            write (string, '(X, A)') "- Flux has been saved into file                                     |"
            call print_ex(string, enable, 'E', iter_interval, i)
        end if

        ! Substitute Phi_next into Phi_prev to calculate the TDGPE
        Phi_prev = Phi_next

        write (string, '(X, A)') "- Finished                                                          |"
        call print_ex(string, enable, 'E', iter_interval, i)

        call change_potential(Pot, Pot_TD, N, i)
    end do
    print *, "---------------------------------------------------------------------"
    print *, ""
    write (*, *) "- All calculation procedures have been finished"
    close (10)
    close (11)
    close (30)
    write (*, *) "- Calculation result has been saved into ", "data.txt"
    write (*, *) output_count, " groups of data"
    print *, "----------------------------------------------------------"
    write (*, *)
    print *, "Result of the calculation ----------------------------------------"
    call expected_value_symm(Phi_prev, H, N, mu, dh)
    print '(X, A, F10.5)', "mu (Chemical Potential) [J] = ", mu
    write (*, *)
end program 