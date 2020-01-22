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
    complex(kind(0d0)),allocatable :: Phi(:, :) ! Wave function phased by given angle
    complex(kind(0d0)),allocatable :: LzPhi(:,:)       ! Angular momentum operator on wave function
    double precision,allocatable   :: Pot(:, :)        ! Potential in Laboratory frame
    double precision,allocatable   :: Flux(:,:,:)
    double precision,allocatable   :: Rot(:,:)
    double precision               :: mu               ! chemical potential
    double precision               :: Lz               ! Angular momentum itself
    double precision               :: t1, t2           ! Calculation time variables
    double precision               :: prob
    ! Coefficients and variables (not user defined)
    integer                        :: i              ! Loop variable
    double precision               :: mu_old           ! Chemical potential at previous step
    allocate (Phi(0:N,0:N), LzPhi(0:N,0:N), Pot(0:N,0:N))
    allocate (Flux(0:N,0:N,1:2), Rot(0:N,0:N))
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
    call output("wavefunction_imaginary_initial.txt", Phi)
    call output("potential_imaginary.txt", Pot)
    mu = solve_energy(Phi, Pot, LzPhi, OMEGA_imag)
    LzPhi = calc_LzPhi(Phi)
    Lz    = calc_Lz(Phi, LzPhi)
    write (*, '(X, A, F13.10, X, A)') "- Angular momentum <Lz> = ", Lz, "hbar"
    write (*, *) "- OMEGA = ", OMEGA_imag
    write (*, '(X, A)') "* Calculation has been initiated"
    call cpu_time(t1)
    do i = 1, 50000
        call evolve(Phi, LzPhi, Pot, OMEGA_imag, 0.5d0, .true.)
        LzPhi = calc_LzPhi(Phi)
        mu_old = mu
        mu = solve_energy(Phi, Pot, LzPhi, OMEGA_imag)
        if (mod(i, 1000) == 0) then
            write (*, '(X, A, I0, A)') "* ", i, " calculations have been done"
        end if
        if (abs(mu_old - mu) < 1d-6) then
            write (*, '(X, A, I0, A)') "- Calculation successfully completed with ", i, " iterations"
            exit
        end if
    end do
    if (i >= 50000) then
        stop "* Calculation has been exceeded its iteration limit. Incorrect result is expected."
    end if
    call output("wavefunction_imaginary_final.txt", Phi)
    Lz = calc_Lz(Phi, LzPhi)
    write (*, '(X, A, F13.10, X, A)') "- Angular momentum <Lz> = ", Lz, "hbar"

    !----------------------- REAL TIME CALCULATION FROM HERE -------------------------------------------
    call initialize(Pot, 5)
    write (*, *) "- OMEGA = ", OMEGA_real
    write (*, '(X, A)') "* Calculation has been initiated"
    open(10, file="wavefunction_real_time_development.txt")
    open(20, file="probability_current_real_time_development.txt")
    open(30, file="rotation_real_time_development.txt")
    do i = 1, 1000
        LzPhi = calc_LzPhi(Phi)
        call evolve(Phi, LzPhi, Pot, OMEGA_real, 0.7d0, .false.)
        if (mod(i, 10) == 0) then
            call output_complex_unit(10, Phi)
            Flux = calc_flux(Phi)
            call output(30, Flux)
            call calc_rotation(Flux, Rot)
            call output_rotation(40, Rot)
        end if
        if (mod(i, 50) == 0) then
            if (i > 50) then
                write (*, '(A)', advance='no') char(13)
            end if
            write (*, '(X, A, I5, A, I5, A)', advance='no') "- ", i, "/", 1000, " completed"
        end if
    end do
    call cpu_time(t2)
    prob = integrate(abs(Phi)**2d0)
    write (*, *)
    write (*, *) "Probability : ", prob
    close (10)
    close (20)
    close (30)
    write (*, '(X, A, F5.1, A)') "- Calculation took ", t2 - t1, " sec"
    LzPhi = calc_LzPhi(Phi)
    Lz = calc_Lz(Phi, LzPhi)
    write (*, '(X, A, F13.10, A)') "Angular Momentum = ", Lz, " hbar"
end program 
