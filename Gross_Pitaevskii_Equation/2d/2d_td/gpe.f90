! Time-development of Gross-Pitaevskii Equation in 2 dimensional space
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
    double precision,allocatable   :: Pot_TD(:, :)     ! Time-dependent Potential
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
    double precision               :: prob
    double precision,allocatable   :: Flux(:,:,:)
    double precision,allocatable   :: Rot(:,:)
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
    allocate (Phi_next(0:N,0:N), Phi_prev(0:N,0:N), Pot(0:N,0:N), j(0:N,0:N), Pot_TD(0:N,0:N))
    allocate (Phi_temp(0:N,0:N))
    allocate (Flux(0:N,0:N,1:2), Rot(0:N,0:N))
    ! Calculation of coefficients and variables using defined physical values
    xmax    = 15d0
    Azero   = sqrt(hbar/(omega_x*mass))
    Xs      = Azero   ! Usually chosen to be Azero for a weak/moderate interaction
    epsilon = (Azero/Xs)**2d0
    kappa   = (4d0*pi*ScatteringLength*ParticleCount/Azero)*(Azero/Xs)**5d0
    dh      = xmax / (n/2 + 0.5d0)
    dt      = 0.01d0*dh*dh
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
    Pot_TD(:, :) = Pot(:, :)

    ! Output the initial wave function to file
    open(10, file="data_initial.txt")
    call output(10, Phi_prev, N, dh, xmax)
    close(10)
    write (*, *) "- Initial wave function has been saved into ", "data_initial.txt"

    ! Start I/O Procedure
    open(10, file="data.txt")
    open(30, file="data_potential.txt")
    open(40, file="data_flux.txt")
    open(50, file="data_rotation.txt")
    !open(11, file="data_current.txt")
    allocate(character(len=80) :: string)
    enable = .true.
    iter_interval = 50
    ! Solve the inconsistent equation until the wave function converges
    do i = 1, 10000
        write (string, '(A, I6, A)') "#", i, " step calculation began ---------------------------------------"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Evolve the system
        call evolve(Phi_prev, N, dt, dh, epsilon, kappa, iu, abs(Phi_prev)**2d0, Pot_TD, Phi_next)
        write (string, '(X, A)') "- The real time evolution has been carried out                      |"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Check if the probability is conserved
        call integrate(abs(Phi_next)**2d0, N, dh, prob)
        write (string, '(X, A, F15.10, A)') "- Probability       = ", prob, "                               |"
        call print_ex(string, enable, 'E', iter_interval, i)

        ! Calculate chemical potential
        mu0 = mu
        call solve_energy(Phi_next, Pot_TD, N, epsilon, kappa, mu, dh)
        write (string, '(X, A, F10.5, A)') "- Chemical Potential = ", mu, "                                   |"
        call print_ex(string, enable, 'E', iter_interval, i)

        if (mod(i, 50) == 0) then
            call output(10, Phi_next, N, dh, xmax)
            write (string, '(X, A)') "- Wave function has been saved into file                            |"
            call print_ex(string, enable, 'E', iter_interval, i)

            ! Output the present potential form into a file
            call output_potential(30, Pot_TD, N, dh, xmax)

            ! Calculate the flux of every discretized point
            call calc_flux(Phi_next, N, mass, dh, hbar, Flux)
            call output_flux(40, Flux, N, dh, xmax)

            ! Calculate the rotation of every discretized point
            call calc_rotation(Flux, N, dh, xmax, Rot)
            call output_rotation(50, Rot, N, dh, xmax)
        end if

        ! Substitute Phi_next into Phi_prev to calculate the TDGPE
        Phi_prev = Phi_next

        write (string, '(X, A)') "- Finished                                                          |"
        call print_ex(string, enable, 'E', iter_interval, i)

        call vary_potential(Pot, Pot_TD, N, dh, i, xmax)

        if (loop_end_flag) then
            exit
        end if
    end do
    print *, "---------------------------------------------------------------------"
    print *, ""
    write (*, *) "- All calculation procedures have been finished"
    close (10)
    !close (11)
    close (30)
    close (40)
    close (50)
    write (*, *) "- Calculation result has been saved into ", "data.txt"
    print *, "----------------------------------------------------------"
    write (*, *)
    print *, "Result of the calculation ----------------------------------------"
    !call expected_value_symm(Phi_prev, H, N, mu)
    print '(X, A, F10.5)', "mu (Chemical Potential) = ", mu
    write (*, *)
    print *, "Wave function half of whose phase is changed by pi is saved into a file => ", "data_shifted.txt"
    !call apply_phase_shift(Phi_prev(floor(N/2d0):N), N-floor(N/2d0)+1, iu, pi, Phi_prev(floor(N/2d0):N))
    !open(11, file="data_shifted.txt")
    !call output(11, Phi_prev, N, dh, xmax)
    !close(11)
end program 
