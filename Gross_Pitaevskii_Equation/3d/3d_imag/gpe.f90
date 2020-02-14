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
    integer                        :: N                   ! Number of division in space
    complex(kind(0d0)),allocatable :: Phi_phased(:, :, :) ! Temporary wave function used by zhbev
    double precision,allocatable   :: Phi_next(:, :, :)   ! Wave function at next step
    double precision,allocatable   :: Phi_prev(:, :, :)   ! Wave function at previous step
    double precision,allocatable   :: Phi_temp(:, :, :)
    double precision,allocatable   :: Pot(:, :, :)        ! Potential
    double precision               :: dh                  ! Step of distance in the x-direction
    double precision               :: dt                  ! Step of time     in the t-direction
    double precision               :: xmax                ! largest x position (Boundary position)
    double precision               :: mass                ! mass of a boson
    double precision               :: omega_x             ! angular velocity of harmonic potential in x direction
    double precision               :: omega_y             ! angular velocity of harmonic potential in y direction
    double precision               :: omega_z             ! angular velocity of harmonic potential in z direction
    double precision               :: gamma_y             ! The ratio of angular velocity in y direction to the one in x direction
    double precision               :: gamma_z             ! The ratio of angular velocity in z direction to the one in x direction
    integer                        :: ParticleCount       ! number of bose particles
    double precision               :: ScatteringLength    ! s-wave scattering length
    double precision               :: mu                  ! chemical potential
    double precision,allocatable   :: j(:, :, :)          ! probability current
    ! Coefficients and variables (not user defined)
    double precision               :: Azero               ! length of the harmonic oscillator ground state
    double precision               :: Xs                  ! characteristic length of the condensate
    double precision               :: epsilon             ! squared ratio of Azero to Xs
    double precision               :: kappa               ! coeffient of the nonlinear term
    integer                        :: i                   ! Loop variable
    double precision               :: mu_old              ! Chemical potential
    double precision,allocatable   :: Phase_field(:,:,:)
    
    ! Output File Path
    character(*),parameter         :: fn_initial           = "data_initial.txt"
    character(*),parameter         :: fn_potential_raw     = "data_potential_raw.txt"
    character(*),parameter         :: fn_potential_cutout  = "data_potential_cutout.txt"
    character(*),parameter         :: fn_result_raw        = "data_raw.txt"
    character(*),parameter         :: fn_result_projection = "data_projection.txt"
    character(*),parameter         :: fn_phased            = "data_phased.txt"

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
    ! Allocation of variables
    allocate (Phi_next(0:N,0:N,0:N), Phi_prev(0:N,0:N,0:N), Pot(0:N,0:N,0:N), j(0:N,0:N,0:N))
    allocate (Phi_phased(0:N, 0:N,0:N), Phi_temp(0:N, 0:N, 0:N), Phase_field(0:N,0:N,0:N))

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
    call normalize(Phi_prev, N, dh)
    write (*, *) "- Initialized wave functions and potential"

    ! Save initial wave function
    open(10, file=fn_initial)
    call output_projection(10, Phi_prev, N, dh, xmax)
    close(10)
    write (*, *) "- Initial Trial Wave Function    => ", fn_initial

    ! Save potential form
    open(10, file=fn_potential_cutout)
    call output_potential_cutout(10, Pot, N, dh, xmax, 25, 25, 25)
    close(10)
    open(10, file=fn_potential_raw)
    call output_potential_raw(10, Pot, N, dh, xmax)
    close(10)
    write (*, *) "- Given Potential Form  (CUTOUT) => ", fn_potential_cutout
    write (*, *) "- Given Potential Form     (RAW) => ", fn_potential_raw

    ! Calculate chemical potential of initial state
    call solve_energy(Phi_prev, Pot, N, epsilon, kappa, mu_old, dh)

    ! Start solving 2D GPE
    do i = 1, 50000
        ! Evolve the system
        call evolve(Phi_prev, N, dt, dh, epsilon, kappa, abs(Phi_prev)**2d0, Pot, Phi_temp)
        Phi_next(:,:,:) = sqrt(0.7d0*abs(Phi_prev(:,:,:))**2d0 + 0.3d0*abs(Phi_next(:,:,:))**2d0)
        call normalize(Phi_next, N, dh)

        ! Calculate chemical potential
        mu_old = mu
        call solve_energy(Phi_next, Pot, N, epsilon, kappa, mu, dh)

        if (mod(i, 100) == 0) then
            write (*, '(X, A, I0, A)') "* ", i, " calculations have done"
        end if

        ! Substitution to calculate self-consistent equation again
        Phi_prev(:,:,:) = Phi_next(:,:,:)

        ! Check if chemical potential has been converged or not
        if (abs(mu_old - mu) < 1d-6) then
            print *, "* Chemical potential has been converged!"
            exit
        end if
    end do
    write (*, '(X, A, I0, A)') "- Calculation successfully completed with ", i, " iterations"
    write (*, *)

    ! Save calculation result
    open(10, file=fn_result_projection)
    call output_projection(10, Phi_prev, N, dh, xmax)
    close (10)
    write (*, *) "- Result Wave Function (PROJECTION) => ", fn_result_projection
    open(11, file=fn_result_raw)
    call output_real(11, Phi_prev, N, dh, xmax)
    close (11)
    write (*, *) "- Result Wave Function (RAW)        => ", fn_result_raw
    
    ! Make a vortex
    call make_vortex(Phi_prev, N, xmax, dh, iu, Phi_phased)
    open(11, file=fn_phased)
    call output_complex(11, Phi_phased, N, dh, xmax)
    close(11)
    write (*, *) "- Phased Wave Function              => ", fn_phased

    call get_phase_field(Phi_phased, N, Phase_field)
    open(11, file="phase.txt")
    call output_real(11, Phase_field, N, dh, xmax)
    close(11)
    write (*, *) "- Phase                             => ", "phase.txt"

    ! Free allocated variables from memory
    deallocate (Phi_next, Phi_prev, Phi_phased, Pot, j)
end program 
