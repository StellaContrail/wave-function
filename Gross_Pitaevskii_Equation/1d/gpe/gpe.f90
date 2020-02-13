! Imaginary time-development of Gross-Pitaevskii Equation in 2 dimensional space
! The Equations and constants are to be referenced to DOI:10.1016/S0021-9991(03)00102-5
! TODO: Store initial values from setting file into variables
! TODO: Load them and if they matches the present settings, skip the imaginary-time calculation part
program main
    use io
    use setting
    use mathf
    use constants
    implicit none
    ! Physical values
    complex(kind(0d0)),allocatable :: Phi(:)        ! Wave function phased by given angle
    double precision  ,allocatable :: Pot(:)        ! Potential in Laboratory frame
    double precision  ,allocatable :: Flux(:)
    double precision               :: mu               ! chemical potential
    double precision               :: t1, t2           ! Calculation time variables
    double precision               :: prob
    integer                        :: limit_iterations
    integer                        :: time_iterations
    integer                        :: time_iterations_interval

    ! Coefficients and variables (not user defined)
    integer                        :: i                ! Loop variable
    double precision               :: mu_old           ! Chemical potential at previous step

    ! File names
    character(*),parameter       :: fn_potential_imaginary                   = "potential_imaginary.txt"
    character(*),parameter       :: fn_potential_real                        = "potential_real.txt"
    character(*),parameter       :: fn_wavefunction_imaginary_initial        = "wavefunction_imaginary_initial.txt"
    character(*),parameter       :: fn_wavefunction_imaginary_result         = "wavefunction_imaginary_final.txt"
    character(*),parameter       :: fn_wavefunction_real_result              = "wavefunction_real_time_development.txt"
    character(*),parameter       :: fn_current_real_result                   = "probability_current_real_time_development.txt"
    character(*),parameter       :: fn_phase_distribution_imag_initial       = "phase_distribution_imag_initial.txt"
    character(*),parameter       :: fn_phase_distribution_imag_result        = "phase_distribution_imag_final.txt"
    character(*),parameter       :: fn_phase_distribution_real_result        = "phase_distribution_real_time_development.txt"
    character(*),parameter       :: fn_energy_iteration_dependance_imaginary = "energy_iteration_dependence_imaginary.txt"
    character(*),parameter       :: fn_energy_iteration_dependance_real      = "energy_iteration_dependence_real.txt"

    ! Allocate variables
    allocate (Phi(0:N), Pot(0:N), Flux(0:N))

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
    print *, "------------------------------------------------------------------"

    write (*, *) "Press Enter to initiate calculation"
    read (*, *)

    call initialize(Pot, 0, Phi)
    call shift_phase(Phi, 1)
    write (*, '(X, A, F0.3, A, F0.3)') "- Phase Shifted at x = ", x0_vortex
    call output(fn_wavefunction_imaginary_initial, Phi)
    call output_potential(fn_potential_imaginary, Pot)
    call output(fn_phase_distribution_imag_initial, phase(Phi))
    open(11, file=fn_energy_iteration_dependance_imaginary)
    mu = solve_energy(Phi, Pot)
    write (11, *) 0, mu, mu*ENERGY_UNIT_IN_DIMENSIONLESS_GPE
    write (*, '(X, A, F0.15, A)') "- mu   = ", mu*ENERGY_UNIT_IN_DIMENSIONLESS_GPE/(1.602d-19), " [eV]" 
    write (*, '(X, A)') "* Calculating 2D GPE Imaginary-time development for ground state"
    call cpu_time(t1)
    limit_iterations = 500000
    do i = 1, limit_iterations
        call evolve(Phi, Pot, .true., 0.7d0)
        mu_old = mu
        mu = solve_energy(Phi, Pot)
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
    call output(fn_phase_distribution_imag_result, phase(Phi)/pi)
    write (*, '(X, A, F0.15, A)') "- mu   = ", mu*ENERGY_UNIT_IN_DIMENSIONLESS_GPE/(1.602d-19), " [eV]" 
    write (*, *) mu
    write (*, *)

    !----------------------- REAL TIME CALCULATION FROM HERE -------------------------------------------
    call initialize(Pot, 0)
    call output_potential(fn_potential_real, Pot)
    write (*, '(X, A)') "* Calculating 2D GPE real-time evoluton of calculated wave function"
    open(10, file=fn_wavefunction_real_result)
    open(20, file=fn_current_real_result)
    open(40, file=fn_phase_distribution_real_result)
    open(50, file=fn_energy_iteration_dependance_real)
    time_iterations = 5000
    time_iterations_interval = 50
    do i = 1, time_iterations
        mu_old = mu
        call evolve(Phi, Pot, .false.)
        mu = solve_energy(Phi, Pot)
        write (50, *) i, mu, mu*ENERGY_UNIT_IN_DIMENSIONLESS_GPE
        if (mod(i, time_iterations_interval) == 0) then
            call output_complex_unit(10, Phi)
            Flux = calc_flux(Phi)
            call output(20, Flux)
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
    close (40)
    close (50)
    write (*, '(X, A, F5.1, A)') "- Calculation took ", t2 - t1, " sec"
    write (*, '(X, A, F0.15, A)') "- mu   = ", mu*ENERGY_UNIT_IN_DIMENSIONLESS_GPE/(1.602d-19), " [eV]"

    ! Output variables data
    open(10, file="gnuplot_vars.txt")
    write (10, '(A)') "# A particle's mass"
    write (10, '(A, X, F0.10)') "mass =", mass
    write (10, '(A)') "# Angular velocity of Harmonic Oscillator Potential"
    write (10, '(A, X, F0.10)') "omega_x =", omega_x
    write (10, '(A)') "# Number of particles in superfluid"
    write (10, '(A, X, I0)') "ParticleCount =", ParticleCount
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
    write (10, '(A)') "# Number of iterations to calculate real-time development"
    write (10, '(A, X, I0)')    "iter =", time_iterations
    write (10, '(A)') "# Number of iterations to skip in real-time development"
    write (10, '(A, X, I0)')    "iter_output =", time_iterations_interval
    write (10, '(A)') "# Distance of how far the system is defined"
    write (10, '(A, X, F0.10)') "xmax =", xmax
    write (10, '(A)') "# File names of saved results"
    write (10, '(3A)')          "fn_potential_imaginary = '", fn_potential_imaginary, "'"
    write (10, '(3A)')          "fn_wavefunction_imaginary_initial = '", fn_wavefunction_imaginary_initial, "'"
    write (10, '(3A)')          "fn_wavefunction_imaginary_result = '", fn_wavefunction_imaginary_result, "'"
    write (10, '(3A)')          "fn_wavefunction_real_result = '", fn_wavefunction_real_result, "'"
    write (10, '(3A)')          "fn_current_real_result = '", fn_current_real_result, "'"
    write (10, '(3A)')          "fn_phase_distribution_real_result = '", fn_phase_distribution_real_result, "'"
    write (10, '(A)') "# Vortex location"
    write (10, '(A, X, F0.10)') "x0_vortex =", x0_vortex
    close(10)

    ! Plot wave function time-lapse (Gnuplot command must be enabled in terminal)
    !write (*, *) "Current animation"
    !call execute_command_line('gnuplot "plot_current.plt"')
    !write (*, *) "Phase animation"
    !call execute_command_line('gnuplot "plot_phase.plt"')
end program 
