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
    integer                        :: N                ! Dimension of Space
    complex(kind(0d0)),allocatable :: Phi_temp(:)   ! Temporary wave function used by zhbev
    complex(kind(0d0)),allocatable :: Phi_next(:)      ! Wave function at next step
    complex(kind(0d0)),allocatable :: Phi_prev(:)      ! Wave function at previous step
    double precision,allocatable   :: Pot(:)           ! Potential
    double precision,allocatable   :: H(:,:)
    double precision               :: dh               ! Step of distance in the x-direction
    double precision               :: dt               ! Step of time     in the t-direction
    double precision               :: xmax             ! largest x position (Boundary position)
    double precision               :: mass             ! mass of a boson
    double precision               :: omega            ! angular velocity of harmonic potential
    double precision               :: ParticleCount    ! number of bose particles
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
    ! Definition of physical values (this could be replaced with I/O)
    ! These values are referenced from
    ! 'Numerical Solution of the Gross-Pitaevskii Equation for Bose-Einstein Condensation'
    ! by Weizhu Bao et al. (2003)
    mass             = 1.44d-25
    omega            = 20d0 * pi
    ParticleCount    = 1d2
    ScatteringLength = 5.1d-9
    N                = 128
    allocate (Phi_next(1:N), Phi_prev(1:N), Pot(1:N), mus(1:N))
    allocate (Phi_temp(1:N), H(1:N,1:N))
    ! Calculation of coefficients and variables using defined physical values
    xmax    = 10d0
    Azero   = sqrt(hbar/(omega*mass))
    Xs      = Azero   ! Usually chosen to be Azero for a weak/moderate interaction
    epsilon = (Azero/Xs)**2d0
    kappa   = (4d0*pi*ScatteringLength*ParticleCount/Azero)*(Azero/Xs)**5d0
    dh      = (2d0*xmax)/N
    dt      = 0.1d0*dh*dh

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
    print *, "dh (A step of distance)           [m] = ", dh
    print *, "dt (A step of time)               [s] = ", dt
    print *, "<Coefficients of NLSE terms>"
    print *, "Epsilon (A0/Xs)^2                     = ", epsilon
    print *, "Kappa (Coefficient of NL term)        = ", kappa
    print *, "------------------------------------------------------------------"
    write (*, *)
    
    print *, "Calculation Start-----------------------------------------"
    ! Initialization of wave function and potential
    call initialize(Phi_next, Phi_prev, Pot, N, dh, xmax)
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
    ! Solve the inconsistent equation until the wave function converges
    do i = 1, 5000
        write (*, '(A, I4, A)') "#", i, " step calculation began --------------"
        ! Construct the hamiltonian using present wave function data
        call hamiltonian(H, Pot, N, dh, epsilon, kappa, Phi_prev)
        print *, "- New hamiltonian has been reconstructed |"

        ! Solve Nonlinear Schroedinger Equation Using the Assumed Wave function
        call exp_mat(H, Phi_prev, N, dt, epsilon, iu, Phi_temp)
        print *, "- NLSE has been successfully calculated  |"

        ! Update the wave function
        Phi_next(:) = Phi_temp(:)
        print *, "- Wave function has been updated         |"

        ! Check if the probability is conserved
        call integrate(abs(Phi_next)**2d0, N, dh, prob)
        print '(X, A, F15.10, A)', "- Probability = ", prob, "          |"

        Phi_prev = Phi_next
        if (mod(i, 50) == 0) then
            call output(10, Phi_prev, N, dh, xmax)
            print *, "- Wave function has been saved into file |"
        end if

        print *, "- Finished                               |"
    end do
    print *, "------------------------------------------"
    print *, ""
    write (*, *) "- All calculation procedures have been finished"
    close (10)
    write (*, *) "- Calculation result has been saved into ", "data.txt"
    print *, "----------------------------------------------------------"
    write (*, *)
    print *, "Result of the calculation ----------------------------------------"
    print *, "mu (Chemical Potential) [J] = "
    write (*, *)
end program 
