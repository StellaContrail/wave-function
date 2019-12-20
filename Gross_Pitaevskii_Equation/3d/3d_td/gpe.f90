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
    complex(kind(0d0)),allocatable :: Phi_phased(:, :, :)! Temporary wave function used by zhbev
    complex(kind(0d0)),allocatable :: Phi_next(:, :, :)! Wave function at next step
    complex(kind(0d0)),allocatable :: Phi_prev(:, :, :)! Wave function at previous step
    double precision,allocatable   :: Pot(:, :, :)     ! Potential
    double precision,allocatable   :: Pot_TD(:,:,:)    ! Time-dependent potential
    double precision,allocatable   :: Flux(:,:,:,:)    ! Probability current
    double precision,allocatable   :: Rot(:,:,:,:)     ! Rotation of flux
    double precision               :: dh               ! Step of distance in the x-direction
    double precision               :: dt               ! Step of time     in the t-direction
    double precision               :: xmax             ! largest x position (Boundary position)
    double precision               :: mass             ! mass of a boson
    double precision               :: omega_x          ! angular velocity of harmonic potential in x direction
    double precision               :: omega_y          ! angular velocity of harmonic potential in y direction
    double precision               :: omega_z          ! angular velocity of harmonic potential in z direction
    double precision               :: gamma_y          ! The ratio of angular velocity in y direction to the one in x direction
    double precision               :: gamma_z          ! The ratio of angular velocity in z direction to the one in z direction
    integer                        :: ParticleCount    ! number of bose particles
    double precision               :: ScatteringLength ! s-wave scattering length
    double precision               :: mu               ! chemical potential
    double precision,allocatable   :: j(:, :, :)       ! probability current
    ! Coefficients and variables (not user defined)
    double precision               :: Azero            ! length of the harmonic oscillator ground state
    double precision               :: Xs               ! characteristic length of the condensate
    double precision               :: epsilon          ! squared ratio of Azero to Xs
    double precision               :: kappa            ! coeffient of the nonlinear term
    integer                        :: i                ! Loop variable
    double precision               :: prob             ! Probability
    double precision               :: prob_old         ! Probability of previous step
    integer(kind=8)                :: t1, t2           ! Processed time
    double precision               :: delta_t          ! delta_t := t2 - t1
    integer(kind=8)                :: tr               ! Time ticks rate
    integer                        :: time_left        ! Estimated time left
    integer                        :: iterations       ! Count of Iterations
    integer                        :: display_output_interval  ! Output progress after every display_output_interval iterations
    
    ! Output File Path
    character(*),parameter         :: fn_initial_proj     = "data_initial_projection.txt"
    character(*),parameter         :: fn_initial_raw      = "data_initial_raw.txt"

    character(*),parameter         :: fn_potential_cutout = "data_potential_cutout.txt"
    character(*),parameter         :: fn_potential_raw    = "data_potential_raw.txt"

    character(*),parameter         :: fn_result_proj      = "data_projection.txt"
    character(*),parameter         :: fn_result_raw       = "data_raw.txt"

    character(*),parameter         :: fn_flux             = "data_flux.txt"
    character(*),parameter         :: fn_rotation         = "data_rotation.txt"

    ! Definition of physical values (this could be replaced with I/O)
    ! These values are referenced from
    ! 'Numerical Solution of the Gross-Pitaevskii Equation for Bose-Einstein Condensation'
    ! by Weizhu Bao et al. (2003)
    mass             = 1.4d-25
    omega_x          = 20d0 * pi
    omega_y          = omega_x
    omega_z          = omega_x
    gamma_y          = omega_y / omega_x
    gamma_z          = omega_z / omega_x
    ParticleCount    = 1000
    ScatteringLength = 5.1d-9

    ! Number of steps in a direction
    N                = 30 - 1
    prob_old         = 1d0
    allocate (Phi_next(0:N,0:N,0:N), Phi_prev(0:N,0:N,0:N), Pot(0:N,0:N,0:N), j(0:N,0:N,0:N))
    allocate (Phi_phased(0:N,0:N,0:N), Pot_TD(0:N,0:N,0:N), Flux(0:N,0:N,0:N,1:3), Rot(0:N,0:N,0:N,1:3))
    
    ! Other variables for setup
    xmax    = 5d0
    Azero   = sqrt(hbar/(omega_x*mass))
    Xs      = Azero
    epsilon = (Azero/Xs)**2d0
    kappa   = (4d0*pi*ScatteringLength*ParticleCount/Azero)*(Azero/Xs)**5d0
    dh      = xmax / (n/2 + 0.5d0)
    dt      = 0.01d0*dh*dh

    ! Display settings
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

    print *, "Press Enter key to start calculation..."
    read (*, *)
    
    ! Initialization
    call initialize(Phi_next, Phi_prev, Pot, N, dh, xmax, gamma_y, gamma_z)
    Pot_TD(:,:,:) = Pot(:,:,:)
    write (*, *) "- Initialized wave functions and potential"

    ! Save initial wave function
    open(10, file=fn_initial_proj)
    call output_projection(10, Phi_prev, N, dh, xmax)
    close(10)
    open(10, file=fn_initial_raw)
    call output_complex_raw(10, Phi_prev, N, dh, xmax)
    close(10)
    write (*, *) "- Initial Wave Function (PROJECTION) => ", fn_initial_proj
    write (*, *) "- Initial Wave Function        (RAW) => ", fn_initial_raw

    ! Calculate probability
    call integrate(abs(Phi_prev)**2d0, N, dh, prob)
    write (*, *) "- Initial Probability        : ", prob
    ! Calculate chemical potential
    call solve_energy(Phi_prev, Pot, N, epsilon, kappa, mu, dh)
    write (*, *) "- Initial chemical potential : ", mu

    open(10, file=fn_result_proj)
    open(15, file=fn_result_raw)
    open(20, file=fn_potential_cutout)
    open(25, file=fn_potential_raw)
    open(30, file=fn_flux)
    open(40, file=fn_rotation)
    call system_clock(t1)
    iterations = 15000
    display_output_interval = 500
    do i = 1, iterations
        ! Evolve the system
        call evolve(Phi_prev, N, dt, dh, epsilon, kappa, iu, abs(Phi_prev)**2d0, Pot_TD, Phi_next)
        Phi_next(:,:,:) = sqrt(0.5d0*abs(Phi_prev(:,:,:))**2d0 + 0.5d0*abs(Phi_next(:,:,:))**2d0)

        if (mod(i, 50) == 0) then
            ! Save wave function
            call output_projection(10, Phi_next, N, dh, xmax)
            call output_complex_raw(15, Phi_next, N, dh, xmax)
            ! Save potential form
            call output_potential_cutout(20, Pot_TD, N, dh, xmax, 25, 25, 25)
            call output_real_raw(25, Pot_TD, N, dh, xmax)
            ! Save flux distribution
            call calc_flux(Phi_next, N, mass, dh, hbar, Flux)
            call output_flux(30, Flux, N, dh, xmax)
            ! Save rotation of flux
            call calc_rotation(Flux, N, dh, xmax, Rot)
            call output_rotation(40, Rot, N, dh, xmax)
        end if

        if (mod(i, display_output_interval) == 0) then
            if (i > display_output_interval) then
                write (*, '(A)', advance='no') char(13)
            end if
            call system_clock(t2, tr)
            delta_t = (t2 - t1)/dble(tr)
            time_left = int(delta_t*((iterations-i)/dble(display_output_interval)))
            write (*, '(X, A, I5, A, I5, A, I5, A)', advance='no') "* ", i, "/", iterations, &
                                                            " completed (", time_left," seconds left)"
            t1 = t2
        end if

        ! Substitution to step forward in time
        Phi_prev(:,:,:) = Phi_next(:,:,:)
        ! Vary potential form depending on time
        call vary_potential(Pot, Pot_TD, N, dh, pi, i, iterations, xmax)
    end do
    write (*, *)
    write (*, *)
    ! Calculate probability
    call integrate(abs(Phi_next)**2d0, N, dh, prob)
    write (*, *) "- Final Probability          : ", prob
    ! Calculate chemical potential
    call solve_energy(Phi_next, Pot_TD, N, epsilon, kappa, mu, dh)
    write (*, *) "- Final chemical potential   : ", mu

    close (10)
    close (20)
    close (30)
    close (40)
    write (*, *) "- Result Wave Function (PROJECTION) => ", fn_result_proj
end program 
