module constants
    implicit none
    ! Mathematical constants
    double precision,parameter     :: pi   = acos(-1d0)       ! PI
    complex(kind(0d0)),parameter   :: iu   = dcmplx(0d0, 1d0) ! Imaginary unit
    
    ! Physical constants
    double precision,parameter     :: hbar = 1.05d-34         ! Reduced Plank constant

    ! Calculation settings
    double precision,parameter     :: mass             = 1.4d-25
    double precision,parameter     :: omega_x          = 20d0 * pi
    double precision,parameter     :: omega_y          = omega_x
    double precision,parameter     :: ParticleCount    = 1000
    double precision,parameter     :: ScatteringLength = 5.1d-9
    integer         ,parameter     :: N                = 100 - 1
    double precision,parameter     :: dh               = 0.16d0
    double precision,parameter     :: dt               = 0.002d0
    double precision,parameter     :: OMEGA            = 100d0 !(2*pi/(dt*5000d0))*10d0

    ! Global constants
    double precision,parameter     :: xmax    = (N/2 + 0.5d0) * dh
    double precision,parameter     :: Azero   = sqrt(hbar/(omega_x*mass))
    double precision,parameter     :: Xs      = Azero
    double precision,parameter     :: gamma   = omega_y / omega_x
    double precision,parameter     :: epsilon = (Azero/Xs)**2d0
    double precision,parameter     :: kappa   = (4d0*pi*ScatteringLength*ParticleCount/Azero)*(Azero/Xs)**5d0
end module