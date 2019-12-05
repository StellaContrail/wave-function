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
    double precision,allocatable   :: Phi_temp(:)      ! Temporary wave function used by zhbev
    double precision,allocatable   :: Phi_next(:)      ! Wave function at next step
    double precision,allocatable   :: Phi_prev(:)      ! Wave function at previous step
    double precision,allocatable   :: Pot(:)           ! Potential
    double precision,allocatable   :: Pot_TD(:)        ! Time-Dependent Potential
    double precision,allocatable   :: H(:,:)           ! Hamiltonian
    double precision,allocatable   :: H_temp(:,:)      ! Temporary Hamiltonian used by zhbev
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
    integer                        :: i,j,k            ! Loop variable
    logical                        :: loop_end_flag    ! Loop end flag
    double precision               :: Degree           ! Degree of asymmetry of wave function
    ! Definition of physical values (this could be replaced with I/O)
    ! These values are referenced from
    ! 'Numerical Solution of the Gross-Pitaevskii Equation for Bose-Einstein Condensation'
    ! by Weizhu Bao et al. (2003)
    mass             = 1.44d-25
    omega            = 20d0 * pi
    ParticleCount    = 1000 ! Until 120 or so is the threshold for precise calculation
    ScatteringLength = 5.1d-9
    !N                = 2**4 - 1!129   ! N must be an odd number
    ! Calculation of coefficients and variables using defined physical values
    xmax    = 5d0
    Azero   = sqrt(hbar/(omega*mass))
    Xs      = 2d0*Azero   ! Usually chosen to be Azero for a weak/moderate interaction
    epsilon = (Azero/Xs)**2d0
    kappa   = (4d0*pi*ScatteringLength*ParticleCount/Azero)*(Azero/Xs)**5d0
    !dh      = xmax / (n/2 + 0.5d0)
    dh      = 0.8d0
    !xmax    = (N+0.5d0)*dh
    N       = int(2d0*xmax/dh)
    if (mod(N, 2) == 0) then
        N = N + 1
    end if
    DIM              = N + 1 ! Include n=0 point
    allocate (Phi_next(0:N), Phi_prev(0:N), Pot(0:N), H(0:N,0:N), Pot_TD(0:N))
    allocate (Phi_temp(1:DIM), H_temp(1:DIM,1:DIM))

    ! Show configuration of fundamental physical constants
    print *, "Physical constants of the system----------------------------------"
    print *, "<Fundamental Physical Constants>"
    print *, "m  (Mass of the bose particle)   [kg] = ", mass
    print *, "omega (Angular velocity of HO)[rad/s] = ", omega
    print *, "N  (Particle Count)           [count] = ", ParticleCount
    print *, "a  (ScatteringLength)             [m] = ", ScatteringLength
    print *, "n  (Dimension of the space)   [count] = ", N
    print *, "A0 (Length of the HO Ground State)[m] = ", Azero
    print *, "Xs (Characteristic Length)        [m] = ", Xs
    print *, "dh (Step of distance)             [m] = ", dh*Xs
    print *, "<Coefficients of NLSE terms>"
    print *, "Epsilon (A0/Xs)^2                     = ", epsilon
    print *, "Kappa (Coefficient of NL term)        = ", kappa
    print *, "<Other Configuration Values>"
    print *, "Delta (4*pi*a*N/a_0)                  = ", (4d0*pi*ScatteringLength*ParticleCount)/Azero
    print *, "Healing length (8*pi*|a|*N/Xs^3)^-0.5 = ", ((8d0*pi*abs(ScatteringLength)*ParticleCount)/(Xs**3d0))**(-0.5d0)
    print *, "------------------------------------------------------------------"
    write (*, *)
    
    print *, "Calculation Start-----------------------------------------"
    ! Initialization of wave function and potential
    call initialize(Phi_next, Phi_prev, Pot, N, dh, xmax)
    write (*, *) "- Initialized the wave function and potential function"

    Pot_TD(:) = Pot(:)

    ! For development usefulness (TO BE DELETED)
    open(20, file="data_asymmetry_degree.txt")

    open(30, file="data_potential.txt")
    call output_potential(30, Pot, N, dh, xmax)
    close(30)

    ! Normalization of wave function
    call normalize(Phi_prev, N, dh)
    
    ! Output the initial wave function to file
    open(10, file="data_initial.txt")
    call output_real(10, Phi_prev, N, dh, xmax)
    close(10)
    write (*, *) "- Initial wave function has been saved into ", "data_initial.txt"

    ! Start I/O Procedure
    open(10, file="data.txt")
    open(30, file="data_td_potential.txt")
    ! Set the loop exit flag to be false
    loop_end_flag = .false.
    do i = 1, 1000
        if (mod(i, 50) == 0) then
            write (*, *) "ITERATION COUNT = ", i
        end if

        call hamiltonian(H, Pot, abs(Phi_prev)**2d0, N, dh, epsilon, kappa)
        do j = 0, N
            do k = 0, N
                H_temp(j+1, k+1) = H(j,k)
            end do
        end do

        Pot_TD(:) =  Pot(:) + kappa*abs(Phi_prev(:))**2d0
        call output_potential(30, Pot_TD, N, dh, xmax)

        call solve_eigen(H_temp(1:DIM, 1:DIM), Phi_temp(1:DIM), mu0, DIM)
        do j = 0, N
            Phi_next(j) = Phi_temp(j+1)
        end do
        call normalize(Phi_next, N, dh)
        Phi_prev(0:N) = Phi_next(0:N)

        ! For development usefulness : TO BE DELETED
        call calc_asymmetry_degree(Phi_prev, N, Degree)
        write (20, *) i, Degree
        write (20, *)

        if (abs(mu0 - mu) < 1d-6) then
            loop_end_flag = .true.
        end if

        mu = mu0

        if (loop_end_flag) then
            exit
        end if
    end do
    close(20)
    write (*, *) "- All calculation procedures have been finished"
    call output_real(10, Phi_prev, N, dh, xmax)
    close (10)
    call calc_asymmetry_degree(Phi_prev, N, Degree)
    write (*, '(X, A, F15.10)') "- Degree of system's asymmetry = ", Degree
    write (*, *) "- Calculation result has been saved into ", "data.txt"
    print *, "----------------------------------------------------------"
    write (*, *)
    print *, "Result of the calculation ----------------------------------------"
    print '(X, A, F9.5)', "mu (Chemical Potential)     = ", mu
    write (*, *)
end program 
