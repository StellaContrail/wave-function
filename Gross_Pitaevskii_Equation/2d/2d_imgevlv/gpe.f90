! Imaginary time-development of Gross-Pitaevskii Equation in 2 dimensional space
! For a reason of Hamiltonian being not dense matrix, we will use just multiple dimensional array here

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
    complex(kind(0d0)),allocatable :: Phi_temp(:, :)   ! Temporary wave function used by zhbev
    complex(kind(0d0)),allocatable :: Phi_next(:, :)   ! Wave function at next step
    complex(kind(0d0)),allocatable :: Phi_prev(:, :)   ! Wave function at previous step
    double precision,allocatable   :: Pot(:, :)        ! Potential
    double precision               :: dh               ! Step of distance in the x-direction
    double precision               :: dt               ! Step of time     in the t-direction
    double precision               :: xmax             ! largest x position (Boundary position)
    double precision               :: mass             ! mass of a boson
    double precision               :: omega_x          ! angular velocity of harmonic potential in x direction
    double precision               :: omega_y          ! angular velocity of harmonic potential in y direction
    double precision               :: gamma            ! The ratio of angular velocity in y direction to the one in x direction
    integer                        :: ParticleCount    ! number of bose particles
    double precision               :: ScatteringLength ! s-wave scattering length
    double precision               :: mu               ! chemical potential
    double precision,allocatable   :: j(:, :)          ! probability current
    ! Coefficients and variables (not user defined)
    double precision               :: Azero            ! length of the harmonic oscillator ground state
    double precision               :: Xs               ! characteristic length of the condensate
    double precision               :: epsilon          ! squared ratio of Azero to Xs
    double precision               :: kappa            ! coeffient of the nonlinear term
    integer                        :: i                ! Loop variable
    character(:),allocatable       :: string           ! ouput string
    logical                        :: enable           ! enable output
    integer                        :: iter_interval    ! output every iter_interval
    logical                        :: loop_end_flag
    double precision               :: mu0
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
    N                = 50 - 1
    allocate (Phi_next(0:N,0:N), Phi_prev(0:N,0:N), Pot(0:N,0:N), j(0:N,0:N))
    allocate (Phi_temp(0:N, 0:N))
    ! Calculation of coefficients and variables using defined physical values
    xmax    = 15d0
    Azero   = sqrt(hbar/(omega_x*mass))
    Xs      = Azero   ! Usually chosen to be Azero for a weak/moderate interaction
    epsilon = (Azero/Xs)**2d0
    kappa   = (4d0*pi*ScatteringLength*ParticleCount/Azero)*(Azero/Xs)**5d0
    dh      = xmax / (n/2 + 0.5d0)
    dt      = 0.1d0*dh*dh
    loop_end_flag = .false.

    ! Show configuration of fundamental physical constants
    print *, "Physical constants of the system----------------------------------"
    print *, "<Fundamental Physical Constants>"
    print *, "m  (Mass of the bose particle)   [kg] = ", mass
    print *, "omegax(Angular velocity of HO)[rad/s] = ", omega_x
    print *, "omegay(Angular velocity of HO)[rad/s] = ", omega_y
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
    ! Initialization of wave functions and potential
    call initialize(Phi_next, Phi_prev, Pot, N, dh, xmax, gamma)
    write (*, *) "- Initialized the wave function and potential function"

    ! Normalization of wave function
    call normalize(Phi_prev, N, dh)

    ! Output the initial wave function to file
    open(10, file="data_initial.txt")
    call output(10, Phi_prev, N, dh, xmax)
    close(10)
    write (*, *) "- Initial wave function has been saved into ", "data_initial.txt"

    ! Output the form of potential to file
    open(10, file="data_pot.txt")
    call output_real(10, Pot, N, dh, xmax)
    close(10)
    write (*, *) "- Potential form has been saved into ", "data_pot.txt"

    ! Start I/O Procedure
    open(10, file="data.txt")
    open(11, file="data_current.txt")
    allocate(character(len=80) :: string)
    enable = .true.
    iter_interval = 50
    ! Solve the inconsistent equation until the wave function converges
    do i = 1, 50000
        write (string, '(A, I6, A)') "#", i, " step calculation began ---------------------------------------"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Evolve the system
        call evolve(Phi_prev, N, dt, dh, epsilon, kappa, abs(Phi_prev)**2d0, Pot, Phi_next)
        write (string, '(X, A)') "- The imaginary time evolution has been carried out                 |"
        call print_ex(string, enable, 'E', iter_interval, i)
        call normalize(Phi_next, N, dh)

        ! Calculate chemical potential
        mu0 = mu
        call solve_energy(Phi_next, Pot, N, epsilon, kappa, mu, dh)
        write (string, '(X, A, F10.5, A)') "- Chemical Potential = ", mu, "                                   |"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Check if chemical potential has been converged or not
        if (abs(mu0 - mu) < 1d-6) then
            print *, "* Chemical potential has been converged!                            |"
            loop_end_flag = .true.
        end if

        ! Substitute Phi_next into Phi_prev to calculate the TDGPE
        Phi_prev = Phi_next

        write (string, '(X, A)') "- Finished                                                          |"
        call print_ex(string, enable, 'E', iter_interval, i)

        if (loop_end_flag) then
            exit
        end if
    end do
    print *, "---------------------------------------------------------------------"
    print *, ""
    write (*, *) "- All calculation procedures have been finished"
    call output(10, Phi_prev, N, dh, xmax)
    close (10)
    write (*, *) "- Calculation result has been saved into ", "data.txt"
    print *, "----------------------------------------------------------"
    write (*, *)
    print *, "Result of the calculation ----------------------------------------"
    print '(X, A, F10.5)', "mu (Chemical Potential) = ", mu
    write (*, *)
    print *, "Wave function half of whose phase is changed by pi is saved into a file => ", "data_shifted.txt"
    call shift_phase(Phi_prev, N, 0, N, ceiling(N/2d0), N, Phi_temp, iu, pi)
    open(11, file="data_shifted.txt")
    call output(11, Phi_temp, N, dh, xmax)
    close(11)
end program 
