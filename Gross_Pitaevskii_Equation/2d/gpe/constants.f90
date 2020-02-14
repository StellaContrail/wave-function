module constants
    implicit none
    ! Mathematical constants
    double precision,  parameter :: pi   = acos(-1d0)       ! PI
    complex(kind(0d0)),parameter :: iu   = dcmplx(0d0, 1d0) ! Imaginary unit
    
    ! Hartree's Atomic Units in SI Units
    double precision,  parameter :: AU_SI_hbar     = 1.05d-34 
    double precision,  parameter :: AU_SI_MASS     = 9.109d-31
    double precision,  parameter :: AU_SI_LENGTH   = 5.291d-11
    double precision,  parameter :: AU_SI_ENERGY   = 4.359d-18
    double precision,  parameter :: AU_SI_TIME     = 1.05d-34 /AU_SI_ENERGY
    double precision,  parameter :: AU_SI_VELOCITY = 2.188d6

    ! Calculation settings
    double precision,  parameter :: hbar             = AU_SI_hbar
    double precision,  parameter :: mass             = 1.45d-25
    double precision,  parameter :: omega_x          = 20d0*pi ! 20d0*pi ! 674.20d0
    double precision,  parameter :: omega_y          = 20d0*pi ! 20d0*pi ! 674.20d0
    double precision,  parameter :: ParticleCount    = 1000
    double precision,  parameter :: ScatteringLength = 109d0*0.529d-10 !109d0*0.529d-10
    ! Dimension
    integer         ,  parameter :: N                = 100 - 1
    ! Space and time step
    double precision,  parameter :: dh               = 0.2d0
    double precision,  parameter :: dt_imag          = 0.01d0
    double precision,  parameter :: dt_real          = 0.01d0
    ! Cranking model's angular velocity (Dimension is the same with omega_x or omega_y)
    double precision,  parameter :: OMEGA_imag       = 2.5d0
    double precision,  parameter :: OMEGA_real       = 0d0
    ! Height of trap potential
    double precision,  parameter :: Vmax             = 200d0
    ! Depth/Height of pinning potential (Dimensionless)
    double precision,  parameter :: V0               = 90d0
    ! Radius of circulary trap potential
    double precision,  parameter :: R0               = 9d0
    ! Other settings about pinning potential
    double precision,  parameter :: delta            = 4d0
    double precision,  parameter :: alpha            = 0d0
    ! Phase is shifted around this location
    double precision,  parameter :: x0_vortex        = 0d0
    double precision,  parameter :: y0_vortex        = 0d0
    ! Pinning site is located at this location
    double precision,  parameter :: x0               = 1.5d0
    double precision,  parameter :: y0               = 1.5d0
    ! whether to be in co-rotating frame with changing phase in phase animation
    logical,           parameter :: onRotating       = .true.

    ! Global constants
    double precision,  parameter :: xmax    = (N/2 + 0.5d0) * dh
    double precision,  parameter :: Azero   = sqrt(hbar/(omega_x*mass))
    double precision,  parameter :: Xs      = Azero
    double precision,  parameter :: gamma   = omega_y / omega_x
    double precision,  parameter :: epsilon = (Azero/Xs)**2d0
    double precision,  parameter :: kappa   = (4d0*pi*ScatteringLength*ParticleCount/Azero)*(Azero/Xs)**5d0
    ! Energy Unit in Dimensionless GPE
    double precision, parameter  :: ENERGY_UNIT_IN_DIMENSIONLESS_GPE = mass*omega_x**2d0*Xs**2d0
end module