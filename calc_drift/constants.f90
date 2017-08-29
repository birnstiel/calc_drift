! ___________________________________________________________________________________________
! set all the constants
! ___________________________________________________________________________________________
MODULE constants
  IMPLICIT NONE
  INTEGER :: n_r
  doubleprecision :: pi,k_b,mu,m_p,Grav,M_sun,AU,sig_sb,sig_h2

  ! first some general constants
  PARAMETER(pi       = 3.141593)        ! PI
  PARAMETER(k_b      = 1.380658d-16)    ! Boltzmann constant in erg/K
  PARAMETER(mu       = 2.3d0)           ! mean molecular mass in proton masses
  PARAMETER(m_p      = 1.6726231d-24)   ! proton mass in g
  PARAMETER(Grav     = 6.67259d-8)      ! gravitational constant in cm^3 g^-1 s^-2
  PARAMETER(M_sun    = 1.989d33)        ! mass of the sun in g
  PARAMETER(AU       = 1.496d13)        ! astronomical unit in cm
  PARAMETER(sig_sb   = 5.670512d-5)     ! Stefan-Boltzmann constant in g s^-3 K^-4
  parameter(sig_h2    = 2d-15)          ! cross section of H2


  ! now the constants of the simulation
  PARAMETER(n_r      = 100)             ! the number of space cells

END MODULE constants
! ===========================================================================================
