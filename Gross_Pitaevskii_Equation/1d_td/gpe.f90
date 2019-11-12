program main
    use io
    use setting
    use mathf
    implicit none
    ! Mathematical constants
    double precision,parameter     :: pi   = acos(-1d0)       ! PI
    complex(kind(0d0)),parameter   :: iu   = dcmplx(0d0, 1d0) ! Imaginary unit
    ! Physical constants
    double precision,parameter     :: hbar = 1d0              ! Reduced Plank constant
    ! Physical values
    integer                        :: N                ! Dimension of Space
    complex(kind(0d0)),allocatable :: Phi_temp(:)      ! Temporary wave function used by zhbev
    complex(kind(0d0)),allocatable :: Phi_next(:)      ! Wave function at next step
    complex(kind(0d0)),allocatable :: Phi_prev(:)      ! Wave function at previous step
    double precision,allocatable   :: Pot(:)           ! Potential
    double precision,allocatable   :: H(:,:)
    double precision               :: dh               ! Step of distance in the x-direction
    double precision               :: dt               ! Step of time     in the t-direction
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
    double precision,allocatable   :: mus(:)           ! multiple chemical potentials returned by zhbev
    integer                        :: i                ! Loop variable
    double precision               :: prob             ! Probability of the system (should be 1)
    character(:),allocatable       :: string           ! ouput string
    logical                        :: enable           ! enable output
    integer                        :: iter_interval    ! output every iter_interval
    ! Definition of physical values (this could be replaced with I/O)
    ! These values are referenced from
    ! 'Numerical Solution of the Gross-Pitaevskii Equation for Bose-Einstein Condensation'
    ! by Weizhu Bao et al. (2003)
    mass             = 1d0
    omega            = 1d0
    ParticleCount    = 100
    ScatteringLength = 5.1d-9
    N                = 129
    allocate (Phi_next(1:N), Phi_prev(1:N), Pot(1:N), mus(1:N))
    allocate (Phi_temp(1:N), H(1:N,1:N))
    ! Calculation of coefficients and variables using defined physical values
    xmax    = 10d0
    Azero   = sqrt(hbar/(omega*mass))
    Xs      = Azero   ! Usually chosen to be Azero for a weak/moderate interaction
    epsilon = (Azero/Xs)**2d0
    kappa   = (4d0*pi*ScatteringLength*ParticleCount/Azero)*(Azero/Xs)**5d0
    dh      = (2d0*xmax)/N
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
    open(11, file="data_flux.txt")
    allocate(character(len=50) :: string)
    enable = .true.
    iter_interval = 50
    ! Solve the inconsistent equation until the wave function converges
    do i = 1, 50000
        write (string, '(A, I6, A)') "#", i, " step calculation began ------------"
        call print_ex(string, enable, 'E', iter_interval, i)
        ! Construct the hamiltonian using present wave function data
        call hamiltonian(H, Pot, N, dh, epsilon, kappa, Phi_prev)
        write (string, '(X, A)') "- New hamiltonian has been reconstructed |"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Solve Nonlinear Schroedinger Equation Using the Assumed Wave function
        call exp_mat(H, Phi_prev, N, dt, epsilon, iu, Phi_temp)
        write (string, '(X, A)') "- NLSE has been successfully calculated  |"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Update the wave function
        Phi_next(:) = Phi_temp(:)
        write (string, '(X, A)') "- Wave function has been updated         |"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Check if the probability is conserved
        call integrate(abs(Phi_next)**2d0, N, dh, prob)
        write (string, '(X, A, F15.10, A)') "- Probability       = ", prob, "    |"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Calculate chemical potential
        call expected_value_symm(Phi_next, H, N, mu)
        write (string, '(X, A, F10.5, A)') "- Chemical Potential = ", mu, "        |"
        call print_ex(string, enable, 'E', iter_interval, i)

        if (mod(i, 50) == 0) then
            call output(10, Phi_next, N, dh, xmax)
            write (string, '(X, A)') "- Wave function has been saved into file |"
            call print_ex(string, enable, 'E', iter_interval, i)

            ! j(x=0) = div[Phi(x=0)] = ∂_x[Phi(x=0)]
            write (11, *) dt*i, (Phi_next(66)-Phi_next(64))/(2d0*dh)
            write (string, '(X, A)') "- Flux has been saved into file          |"
            call print_ex(string, enable, 'E', iter_interval, i)
        end if

        ! Substitute Phi_next into Phi_prev to calculate the TDGPE
        Phi_prev = Phi_next

        write (string, '(X, A)') "- Finished                               |"
        call print_ex(string, enable, 'E', iter_interval, i)
    end do
    print *, "------------------------------------------"
    print *, ""
    write (*, *) "- All calculation procedures have been finished"
    close (10)
    close (11)
    write (*, *) "- Calculation result has been saved into ", "data.txt"
    print *, "----------------------------------------------------------"
    write (*, *)
    print *, "Result of the calculation ----------------------------------------"
    call expected_value_symm(Phi_prev, H, N, mu)
    print '(X, A, F10.5)', "mu (Chemical Potential) [J] = ", mu
    write (*, *)
end program 
