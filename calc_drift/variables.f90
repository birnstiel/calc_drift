MODULE variables
  USE constants,ONLY:n_r,M_sun,pi
  IMPLICIT NONE
  SAVE

  ! note: everything in cgs-units
  doubleprecision :: M_star             = M_sun   ! the mass of the star
  doubleprecision :: alpha(1:n_r)       = 0.01    ! the array for the alpha parameter
  doubleprecision :: dust_floor         = 1d-100  ! minimum value of the dust
  doubleprecision :: gas2dust_ratio     = 100d0   ! gas to dust ratio
  doubleprecision :: Sc                 = 1d0     ! Schmidt number
  doubleprecision :: stokes_factor      = pi/2d0
  INTEGER         :: peak_position(1)   = 1

END MODULE variables
