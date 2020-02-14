! Imaginary time-development of Gross-Pitaevskii Equation in 2 dimensional space
! The Equations and constants are to be referenced to DOI:10.1016/S0021-9991(03)00102-5

program main
    use io
    use setting
    use mathf
    use constants
    implicit none
    ! Physical values
    complex(kind(0d0)),allocatable :: Phi(:, :)        ! Wave function phased by given angle
    complex(kind(0d0)),allocatable :: LzPhi(:,:)       ! Angular momentum operator on wave function
    double precision  ,allocatable :: Pot(:, :)        ! Potential in Laboratory frame
    double precision  ,allocatable :: Flux(:,:,:)
    double precision  ,allocatable :: Rot(:,:)
    double precision               :: mu               ! chemical potential
    double precision               :: Lz               ! Angular momentum itself
    double precision               :: t1, t2           ! Calculation time variables
    double precision               :: prob
    integer                        :: limit_iterations
    integer                        :: time_iterations
    integer                        :: time_iterations_interval
    double precision               :: n_flux, n_phase
    integer                        :: n_vortex, imag_pot_mode, real_pot_mode

    ! Coefficients and variables (not user defined)
    integer                        :: i,j              ! Loop variable
    double precision               :: mu_old           ! Chemical potential at previous step

    ! File names
    character(*),parameter       :: fn_potential_imaginary                   = "potential_imaginary.txt"
    character(*),parameter       :: fn_potential_real                        = "potential_real.txt"
    character(*),parameter       :: fn_wavefunction_imaginary_initial        = "wavefunction_imaginary_initial.txt"
    character(*),parameter       :: fn_wavefunction_imaginary_result         = "wavefunction_imaginary_final.txt"
    character(*),parameter       :: fn_wavefunction_real_result              = "wavefunction_real_time_development.txt"
    character(*),parameter       :: fn_current_imag_result                   = "probability_current_imag_result.txt"
    character(*),parameter       :: fn_current_real_result                   = "probability_current_real_time_development.txt"
    character(*),parameter       :: fn_rotation_real_result                  = "rotation_real_time_development.txt"
    character(*),parameter       :: fn_phase_distribution_imag_result        = "phase_distribution_imaginary.txt"
    character(*),parameter       :: fn_phase_distribution_real_result        = "phase_distribution_real_time_development.txt"
    character(*),parameter       :: fn_energy_iteration_dependance_imaginary = "energy_iteration_dependence_imaginary.txt"
    character(*),parameter       :: fn_energy_iteration_dependance_real      = "energy_iteration_dependence_real.txt"

    ! Allocate variables
    allocate (Phi(0:N,0:N), LzPhi(0:N,0:N), Pot(0:N,0:N))
    allocate (Flux(0:N,0:N,1:2), Rot(0:N,0:N))

    ! Write out present settings
    print *, "Physical constants of the system----------------------------------"
    print *, "<Fundamental Physical Constants>"
    print *, "n         (Dimension of the space) = ", N + 1
    print *, "dh              (Step of distance) = ", dh
    print *, "dt_imag             (Step of time) = ", dt_imag
    print *, "dt_real             (Step of time) = ", dt_real
    print *, "N                 (Particle Count) = ", ParticleCount
    print *, "a               (ScatteringLength) = ", ScatteringLength
    print *, "Xs         (Characteristic Length) = ", Xs
    print *, "hbar      (Reduced Plank Constant) = ", hbar
    print *, "omega_x   (Angular velocity of HO) = ", omega_x
    print *, "omega_y   (Angular velocity of HO) = ", omega_y
    print *, "m      (Mass of the bose particle) = ", mass
    print *, "xmax (Absolute value of maximum x) = ", xmax
    print *, "A0 (Length of the HO Ground State) = ", Azero
    print *, ""
    print *, "<Coefficients of NLSE terms>"
    print *, "Epsilon                  (A0/Xs)^2 = ", epsilon
    print *, "Kappa     (Coefficient of NL term) = ", kappa
    print *, ""
    print *, "<Other Useful Values>"
    print *, "4*pi*a*N/A0                (delta) = ", (4d0*pi*ScatteringLength*ParticleCount)/Azero
    print *, "1/sqrt(8*pi*|a|*N/Xs^3)      (x_h) = ", ((8d0*pi*abs(ScatteringLength)*ParticleCount)/(Xs**3d0))**(-0.5d0)
    print *, "m*omega_x^2*(Xs)^0.5               = ", mass*(omega_x**2d0)*sqrt(Xs)
    print *, "epsilon*(OMEGA_imag/omega_x)       = ", epsilon*(OMEGA_imag/omega_x)
    print *, "epsilon*(OMEGA_real/omega_x)       = ", epsilon*(OMEGA_real/omega_x)
    print *, "------------------------------------------------------------------"

    ! .TRUE.  => Solve Quantities
    ! .FALSE. => Solve time development
    n_vortex = 1
    imag_pot_mode = 6
    real_pot_mode = 5
    if (.false.) then
        write (*, *) "Calculation initiated."
        open (12, file="quantities.txt", action="write")
        open (11, file=fn_energy_iteration_dependance_imaginary)
        do j = 0, 40
            call initialize(Pot, imag_pot_mode, Phi)
            call make_vortex(Phi, n_vortex)
            LzPhi = calc_LzPhi(Phi)
            Lz    = calc_Lz(Phi, LzPhi, .true.)
            mu = solve_energy(Phi, Pot, LzPhi, 0.25d0*j)
            write (11, *) 0, mu
            limit_iterations = 500000
            do i = 1, limit_iterations
                LzPhi = calc_LzPhi(Phi)
                call evolve(Phi, LzPhi, Pot, 0.25d0*j, .true., 0.7d0)
                mu_old = mu
                mu = solve_energy(Phi, Pot, LzPhi, 0.25d0*j)
                write (11, *) i, mu
                if (abs(mu_old - mu) < 1d-10) then
                    write (*, '(A)', advance='no') char(13)
                    write (*, '(F0.10, X, F0.11)', advance='no') 0.25d0*j, abs(mu_old - mu)
                    write (*, *)
                    write (*, '(X, A, I0, A)') "- Calculation successfully completed with ", i, " iterations"
                    exit
                else
                    if (mod(i,100) == 0) then
                        write (*, '(A)', advance='no') char(13)
                        write (*, '(F0.10, X, F0.10)', advance='no') 0.25d0*j, abs(mu_old - mu)
                    end if
                end if
            end do
            if (i >= limit_iterations) then
                write (*, *)
                write (*, *) "* Calculation has been exceeded its iteration limit. Incorrect result is expected."
                close (11)
                close (12)
            end if
            LzPhi = calc_LzPhi(Phi)
            Lz = calc_Lz(Phi, LzPhi, .true.)
            Flux = calc_flux(Phi)
            n_flux  = circulation(Phi, Flux) / (2d0*pi)
            n_phase = circulation(Phi) / (2d0*pi)
            write (*, '(F0.10, X, 5(F0.5, X))') 0.25d0*j, Lz, mu, mu*ENERGY_UNIT_IN_DIMENSIONLESS_GPE/(1.602d-19), n_phase, n_flux

            write (12, *) 0.25d0*j, Lz, mu, mu*ENERGY_UNIT_IN_DIMENSIONLESS_GPE/(1.602d-19), n_phase, n_flux
        end do
        close (11)
        close (12)

    else
        write (*, *) "Press Enter to initiate calculation"
        read (*, *)
        call initialize(Pot, imag_pot_mode, Phi)
        call make_vortex(Phi, n_vortex)
        write (*, '(X, A, F0.3, A, F0.3, A)') "- Phase Shifted at (",x0_vortex, ", ", y0_vortex, ")"
        write (*, '(X, A, F0.10)') "- OMEGA = ", OMEGA_imag
        call output(fn_wavefunction_imaginary_initial, Phi)
        call output_potential(fn_potential_imaginary, Pot)
        open(11, file=fn_energy_iteration_dependance_imaginary)
        LzPhi = calc_LzPhi(Phi)
        Flux = calc_flux(Phi)
        mu = solve_energy(Phi, Pot, LzPhi, OMEGA_imag)
        write (11, *) 0, mu, mu*ENERGY_UNIT_IN_DIMENSIONLESS_GPE
        Lz    = calc_Lz(Phi, LzPhi, .true.)
        write (*, '(X, A, F0.10, X, A)') "- <Lz> = ", Lz, "hbar"
        write (*, '(X, A, F0.15, A)') "- mu   = ", mu*ENERGY_UNIT_IN_DIMENSIONLESS_GPE/(1.602d-19), " [eV]" 
        write (*, '(X, A)') "* Calculating 2D GPE Imaginary-time development for ground state"
        n_flux = circulation(Phi, Flux) / (2d0*pi)
        n_phase = circulation(Phi) / (2d0*pi)
        write (*, '(X, A, F0.10, A, F0.10)') "- n    = ",  n_flux, " ERROR =>", abs(n_vortex-n_flux)
        write (*, '(X, A, F0.10, A, F0.10)') "- n    = ", n_phase, " ERROR =>", abs(n_vortex-n_phase)
        call cpu_time(t1)
        limit_iterations = 500000
        do i = 1, limit_iterations
            LzPhi = calc_LzPhi(Phi)
            call evolve(Phi, LzPhi, Pot, OMEGA_imag, .true., 0.7d0)
            mu_old = mu
            mu = solve_energy(Phi, Pot, LzPhi, OMEGA_imag)
            write (11, *) i, mu, mu*ENERGY_UNIT_IN_DIMENSIONLESS_GPE
            if (mod(i, 500) == 0) then
                if (i > 500) then
                    write (*, '(A)', advance='no') char(13)
                end if
                write (*, '(X, A, I7, A, F0.9)', advance='no') "- ", i, " calculations have been done ", abs(mu_old-mu)
            end if
            if (abs(mu_old - mu) < 1d-10) then
                write (*, *)
                write (*, '(X, A, I0, A)') "- Calculation successfully completed with ", i, " iterations"
                exit
            end if
        end do
        close(11)
        if (i >= limit_iterations) then
            write (*, *)
            stop "* Calculation has been exceeded its iteration limit. Incorrect result is expected."
        end if
        call output(fn_wavefunction_imaginary_result, Phi)
        call output(fn_phase_distribution_imag_result, phase(Phi))
        LzPhi = calc_LzPhi(Phi)
        Lz = calc_Lz(Phi, LzPhi, .true.)
        Flux = calc_flux(Phi)
        call output(fn_current_imag_result, Flux)
        n_flux  = circulation(Phi, Flux) / (2d0*pi)
        n_phase = circulation(Phi) / (2d0*pi)
        write (*, '(X, A, F0.10, A)')    "- <Lz>    = ", Lz, " hbar"
        write (*, '(X, A, F0.15, A)')    "- mu      = ", mu*ENERGY_UNIT_IN_DIMENSIONLESS_GPE/(1.602d-19), " [eV]"
        write (*, '(X, A, F0.10, X, A)') "- n_flux  = ", circulation(Phi, Flux) / (2d0*pi)
        write (*, '(X, A, F0.10, X, A)') "- n_phase = ", circulation(Phi) / (2d0*pi)
        write (*, *)
        
        !call make_vortex(Phi, n_vortex)
    end if

    !----------------------- REAL TIME CALCULATION FROM HERE -------------------------------------------
    if (.true.) then
    call initialize(Pot, real_pot_mode)
    call output_potential(fn_potential_real, Pot)
    write (*, '(X, A)') "* Calculating 2D GPE real-time evoluton of calculated wave function"
    write (*, '(X, A, F0.10)') "- OMEGA = ", OMEGA_real
    open(10, file=fn_wavefunction_real_result)
    open(20, file=fn_current_real_result)
    open(30, file=fn_rotation_real_result)
    open(40, file=fn_phase_distribution_real_result)
    open(50, file=fn_energy_iteration_dependance_real)
    time_iterations = 5000
    time_iterations_interval = 50
    do i = 1, time_iterations
        LzPhi = calc_LzPhi(Phi)
        mu_old = mu
        call evolve(Phi, LzPhi, Pot, OMEGA_real, .false.)
        mu = solve_energy(Phi, Pot, LzPhi, OMEGA_real)
        write (50, *) i, mu, mu*ENERGY_UNIT_IN_DIMENSIONLESS_GPE
        if (mod(i, time_iterations_interval) == 0) then
            call output_complex_unit(10, Phi)
            Flux = calc_flux(Phi)
            call output(20, Flux)
            call calc_rotation(Flux, Rot)
            call output(30, Rot)
            if (onRotating) then
                call output(40, phase(Phi*exp(iu*mod(mu_old*dt_real*i/epsilon,2d0*pi))))
            else
                call output(40, phase(Phi))
            end if
        end if
        if (mod(i, 100) == 0) then
            if (i > 100) then
                write (*, '(A)', advance='no') char(13)
            end if
            write (*, '(X, A, I5, A, I5, A)', advance='no') "- ", i, "/", time_iterations, " completed"
        end if
    end do
    call cpu_time(t2)
    write (*, *)
    write (*, '(X, A, F0.10)') "- Probability = ", integrate(abs(Phi)**2d0)
    close (10)
    close (20)
    close (30)
    close (40)
    close (50)
    write (*, '(X, A, F5.1, A)') "- Calculation took ", t2 - t1, " sec"
    LzPhi = calc_LzPhi(Phi)
    Lz = calc_Lz(Phi, LzPhi, .false.)
    write (*, '(X, A, F0.10, A)')    "- <Lz>    = ", Lz, " hbar"
    write (*, '(X, A, F0.15, A)')    "- mu      = ", mu*ENERGY_UNIT_IN_DIMENSIONLESS_GPE/(1.602d-19), " [eV]"
    write (*, '(X, A, F0.10, X, A)') "- n_flux  = ", circulation(Phi, Flux) / (2d0*pi)
    write (*, '(X, A, F0.10, X, A)') "- n_phase = ", circulation(Phi) / (2d0*pi)
    end if
    ! Output variables data
    open(10, file="gnuplot_vars.txt")
    write (10, '(A)') "# A particle's mass"
    write (10, '(A, X, F0.10)') "mass =", mass
    write (10, '(A)') "# Angular velocity of Harmonic Oscillator Potential"
    write (10, '(A, X, F0.10)') "omega_x =", omega_x
    write (10, '(A, X, F0.10)') "omega_y =", omega_y
    write (10, '(A)') "# Number of particles in superfluid"
    write (10, '(A, X, F0.10)') "ParticleCount =", ParticleCount
    write (10, '(A)') "# S-wave scattering length"
    write (10, '(A, X, F0.10)') "ScatteringLength =", ScatteringLength
    write (10, '(A)') "# Characteristic length"
    write (10, '(A, X, F0.10)') "CharacteristicLength =", Xs
    write (10, '(A)') "# Space step"
    write (10, '(A, X, F0.10)') "dh =", dh
    write (10, '(A)') "# Time step"
    write (10, '(A, X, F0.10)') "dt =", dt_real
    write (10, '(A)') "# Dimension subtracted by one"
    write (10, '(A, X, I0)')    "N =", N
    write (10, '(A)') "# Pinning Location"
    write (10, '(A, X, F0.10)') "x0 =", x0
    write (10, '(A, X, F0.10)') "y0 =", y0
    write (10, '(A)') "# Number of iterations to calculate real-time development"
    write (10, '(A, X, I0)')    "iter =", time_iterations
    write (10, '(A)') "# Number of iterations to skip in real-time development"
    write (10, '(A, X, I0)')    "iter_output =", time_iterations_interval
    write (10, '(A)') "# Distance of how far the system is defined"
    write (10, '(A, X, F0.10)') "xmax =", xmax
    write (10, '(A)') "# File names of saved results"
    write (10, '(3A)')          "fn_potential_real ='", fn_potential_real, "'" 
    write (10, '(3A)')          "fn_potential_imaginary = '", fn_potential_imaginary, "'"
    write (10, '(3A)')          "fn_wavefunction_imaginary_initial = '", fn_wavefunction_imaginary_initial, "'"
    write (10, '(3A)')          "fn_wavefunction_imaginary_result = '", fn_wavefunction_imaginary_result, "'"
    write (10, '(3A)')          "fn_wavefunction_real_result = '", fn_wavefunction_real_result, "'"
    write (10, '(3A)')          "fn_current_imag_result = '", fn_current_imag_result, "'"
    write (10, '(3A)')          "fn_current_real_result = '", fn_current_real_result, "'"
    write (10, '(3A)')          "fn_rotation_real_result = '", fn_rotation_real_result, "'"
    write (10, '(3A)')          "fn_phase_distribution_imag_result = '", fn_phase_distribution_imag_result, "'"
    write (10, '(3A)')          "fn_phase_distribution_real_result = '", fn_phase_distribution_real_result, "'"
    write (10, '(A)') "# Angular velocity of Cranking model in each (imaginary/real)-time development"
    write (10, '(A, X, F0.10)') "OMEGA_imag =", OMEGA_imag
    write (10, '(A, X, F0.10)') "OMEGA_real =", OMEGA_real
    write (10, '(A)') "# Circularly-Symmetric Trap Potential Settings"
    write (10, '(A)') "# - Height"
    write (10, '(A, X, F0.10)') "Vmax =", Vmax
    write (10, '(A)') "# - Radius"
    write (10, '(A, X, F0.10)') "R0 =", R0
    write (10, '(A)') "# Pinning Potential Setings"
    write (10, '(A)') "# - Depth"
    write (10, '(A, X, F0.10)') "V0 =", V0
    write (10, '(A)') "# - How steep the pinning's slope is"
    write (10, '(A, X, F0.10)') "delta =", delta
    write (10, '(A)') "# - How wide the pinning is"
    write (10, '(A, X, F0.10)') "alpha =", alpha
    write (10, '(A)') "# Vortex location"
    write (10, '(A, X, F0.10)') "x0_vortex =", x0_vortex
    write (10, '(A, X, F0.10)') "y0_vortex =", y0_vortex
    close(10)

    ! Plot wave function time-lapse (Gnuplot command must be enabled in terminal)
    write (*, *) "Current animation"
    call execute_command_line('gnuplot "plot_current.plt"')
    !write (*, *) "Phase animation"
    call execute_command_line('gnuplot "plot_phase.plt"')
end program 
