MODULE variables
  USE constants,ONLY:M_sun,pi
  IMPLICIT NONE
  SAVE

  ! note: everything in cgs-units
  doubleprecision :: M_star             = M_sun   ! the mass of the star
  doubleprecision :: Sc                 = 1d0     ! Schmidt number
  doubleprecision :: rho_s              = 1.2d0   ! internal density of the grains
  doubleprecision :: stokes_factor      = pi/2d0
  INTEGER         :: peak_position(1)   = 1

END MODULE variables
