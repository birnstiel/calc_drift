module switches
	implicit none
	save

    integer	        :: dust_diffusion     = 0d0
    integer	        :: dust_radialdrift   = 1d0
    integer	        :: dust_drag          = 0d0
    doubleprecision :: drift_fudge_factor = 1d0 ! tune the strength of radial drift, 1=normal

end module switches
