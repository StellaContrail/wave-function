! Imaginary time-development of Gross-Pitaevskii Equation in 2 dimensional space
! The Equations and constants are to be referenced to DOI:10.1016/S0021-9991(03)00102-5
! TODO: Store initial values from setting file into variables

program main
    use io
    use setting
    use mathf
    use constants
    implicit none
    ! Physical values
    complex(kind(0d0)),allocatable :: Phi(:, :)        ! Wave function phased by given angle
    complex(kind(0d0)),allocatable :: LzPhi(:,:)       ! Angular momentum operator on wave function
    double precision,allocatable   :: Pot(:, :)        ! Potential in Laboratory frame
    double precision,allocatable   :: Flux(:,:,:)
    double precision,allocatable   :: Rot(:,:)
    double precision               :: mu               ! chemical potential
    double precision               :: Lz               ! Angular momentum itself
    double precision               :: t1, t2           ! Calculation time variables
    double precision               :: prob
    integer                        :: time_iterations
    integer                        :: time_iterations_interval

    ! Coefficients and variables (not user defined)
    integer                        :: i                ! Loop variable
    double precision               :: mu_old           ! Chemical potential at previous step

    ! File names
    character(*),parameter       :: fn_potential_imaginary = "potential_imaginary.txt"
    character(*),parameter       :: fn_wavefunction_imaginary_initial = "wavefunction_imaginary_initial.txt"
    character(*),parameter       :: fn_wavefunction_imaginary_result = "wavefunction_imaginary_final.txt"
    character(*),parameter       :: fn_wavefunction_real_result = "wavefunction_real_time_development.txt"
    character(*),parameter       :: fn_current_real_result = "probability_current_real_time_development.txt"
    character(*),parameter       :: fn_rotation_real_result = "rotation_real_time_development.txt"

    ! Allocate variables
    allocate (Phi(0:N,0:N), LzPhi(0:N,0:N), Pot(0:N,0:N))
    allocate (Flux(0:N,0:N,1:2), Rot(0:N,0:N))

    ! Write out present settings
    print *, "Physical constants of the system----------------------------------"
    print *, "<Fundamental Physical Constants>"
    print *, "m  (Mass of the bose particle)   [kg] = ", mass
    print *, "omegax(Angular velocity of HO)[rad/s] = ", omega_x
    print *, "omegay(Angular velocity of HO)[rad/s] = ", omega_y
    print *, "N (Particle Count)            [count] = ", ParticleCount
    print *, "a (ScatteringLength)              [m] = ", ScatteringLength
    print *, "n (Dimension of the space)    [count] = ", N + 1
    print *, "A0 (Length of the HO Ground State)[m] = ", Azero
    print *, "Xs (Characteristic Length)        [m] = ", Xs
    print *, "dh (Step of distance)                 = ", dh
    print *, "dt (Step of time)                     = ", dt
    print *, "xmax (Absolute value of maximum x)[m] = ", xmax
    print *, "<Coefficients of NLSE terms>"
    print *, "Epsilon (A0/Xs)^2                     = ", epsilon
    print *, "Kappa (Coefficient of NL term)        = ", kappa
    print *, "<Other Configuration Values>"
    print *, "Delta (4*pi*a*N/a_0)                  = ", (4d0*pi*ScatteringLength*ParticleCount)/Azero
    print *, "Healing length (8*pi*|a|*N/Xs^3)^-0.5 = ", ((8d0*pi*abs(ScatteringLength)*ParticleCount)/(Xs**3d0))**(-0.5d0)
    print *, "------------------------------------------------------------------"
    call initialize(Pot, 5, Phi)
    !call make_vortex(Phi, 1d0)
    call output(fn_wavefunction_imaginary_initial, Phi)
    call output(fn_potential_imaginary, Pot)
    mu = solve_energy(Phi, Pot, LzPhi, OMEGA_imag)
    LzPhi = calc_LzPhi(Phi)
    Lz    = calc_Lz(Phi, LzPhi, .true.)
    write (*, '(X, A, F0.10, X, A)') "- <Lz> = ", Lz, "hbar"
    write (*, '(X, A)') "* Calculating 2D GPE Imaginary-time development for ground state"
    write (*, '(X, A, F0.10)') "- OMEGA = ", OMEGA_imag
    call cpu_time(t1)
    do i = 1, 50000
        LzPhi = calc_LzPhi(Phi)
        call evolve(Phi, LzPhi, Pot, OMEGA_imag, .true., 0.7d0)
        mu_old = mu
        mu = solve_energy(Phi, Pot, LzPhi, OMEGA_imag)
        if (mod(i, 500) == 0) then
            if (i > 500) then
                write (*, '(A)', advance='no') char(13)
            end if
            write (*, '(X, A, I5, A)', advance='no') "- ", i, " calculations have been done"
        end if
        if (abs(mu_old - mu) < 1d-10) then
            write (*, *)
            write (*, '(X, A, I0, A)') "- Calculation successfully completed with ", i, " iterations"
            exit
        end if
    end do
    if (i >= 50000) then
        stop "* Calculation has been exceeded its iteration limit. Incorrect result is expected."
    end if
    call make_vortex(Phi, 1.5d0)
    write (*, '(X, A)') "- Made a vortex"
    call output(fn_wavefunction_imaginary_result, Phi)
    LzPhi = calc_LzPhi(Phi)
    Lz = calc_Lz(Phi, LzPhi, .true.)
    write (*, '(X, A, F0.10, X, A)') "- <Lz> = ", Lz, "hbar"
    write (*, *)

    !----------------------- REAL TIME CALCULATION FROM HERE -------------------------------------------
    call initialize(Pot, 5)
    write (*, '(X, A)') "* Calculating 2D GPE real-time evoluton of calculated wave function"
    write (*, '(X, A, F0.5)') "- OMEGA = ", OMEGA_real
    open(10, file=fn_wavefunction_real_result)
    open(20, file=fn_current_real_result)
    open(30, file=fn_rotation_real_result)
    time_iterations = 1000
    time_iterations_interval = 10
    do i = 1, time_iterations
        LzPhi = calc_LzPhi(Phi)
        call evolve(Phi, LzPhi, Pot, OMEGA_real, .false.)
        if (mod(i, time_iterations_interval) == 0) then
            call output_complex_unit(10, Phi)
            Flux = calc_flux(Phi)
            call output(30, Flux)
            call calc_rotation(Flux, Rot)
            call output(30, Rot)
        end if
        if (mod(i, 50) == 0) then
            if (i > 50) then
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
    write (*, '(X, A, F5.1, A)') "- Calculation took ", t2 - t1, " sec"
    LzPhi = calc_LzPhi(Phi)
    Lz = calc_Lz(Phi, LzPhi, .false.)
    write (*, '(X, A, F0.10, A)') "- <Lz> = ", Lz, " hbar"

    ! Output variables data
    open(10, file="gnuplot_vars.txt")
    write (10, '(A, X, F0.10)') "dh =", dh
    write (10, '(A, X, F0.10)') "dt =", dt
    write (10, '(A, X, I0)')    "N =", N
    write (10, '(A, X, F0.10)') "x0 =", x0
    write (10, '(A, X, F0.10)') "y0 =", y0
    write (10, '(A, X, I0)')    "iter =", time_iterations
    write (10, '(A, X, I0)')    "iter_output =", time_iterations_interval
    write (10, '(A, X, F0.10)') "xmax =", xmax
    write (10, '(3A)')          "fn_potential_imaginary = '", fn_potential_imaginary, "'"
    write (10, '(3A)')          "fn_wavefunction_imaginary_initial = '", fn_wavefunction_imaginary_initial, "'"
    write (10, '(3A)')          "fn_wavefunction_imaginary_result = '", fn_wavefunction_imaginary_result, "'"
    write (10, '(3A)')          "fn_wavefunction_real_result = '", fn_wavefunction_real_result, "'"
    write (10, '(3A)')          "fn_current_real_result = '", fn_current_real_result, "'"
    write (10, '(3A)')          "fn_rotation_real_result = '", fn_rotation_real_result, "'"
    close(10)

    ! Plot wave function time-lapse
    write (*, *) "Generating animation of wave function's time evolution"
    call execute_command_line('gnuplot "plot.plt"')
end program 
