! _____________________________________________________________________________
! this subroutine does the implicit timestep for a dust layer
! it uses a implicit donor cell advection-diffusion scheme with piecewise
! constant "interpolation"
!
! INPUT:         dt               = the length of the timestep
!                rho_s            = solid state density of the dust
!                grainsize        = size of the dust grains
!                x()              = radial grid
!                x05()            = radial half step grid
!                sigma_g()        = gas surface denisty
!                nu()             = viscosity
!                T()              = midplane temperature
!                sigma_dot()      = infall rate
!                sigma1()         = old dust surface density
!
! RETURNS:       sigma2()         = the updated dust surface density
!                dust_accretion   = dust accretion rate onto the star
!                dust_flux_o      = dust "accretion" to universe outside the grid
!                dust_flux()      = flux at each grid point
! _____________________________________________________________________________
subroutine duststep_donorcell(dt,grainsize,grainmass,x,x05,sigma_g,T,sigma_dot,sigma1,sigma2,dust_accretion,dust_flux_o &
                                &,v_drift,v_05,flim,Diff,v_gas,dust_flux,&
                                &coagulation_method,A,B,C,D)
use constants, ONLY:n_r,pi,Grav,mu,m_p,k_b
use variables, ONLY:M_star,alpha,Sc,peak_position
use switches,  ONLY:dust_diffusion,dust_radialdrift,dust_drag,drift_fudge_factor
implicit none
doubleprecision, intent(in)    :: dt,grainsize,grainmass
doubleprecision, intent(in)    :: x(1:n_r),x05(1:n_r),sigma_g(1:n_r),T(1:n_r),sigma_dot(1:n_r),v_gas(1:n_r)
doubleprecision, intent(inout) :: sigma1(1:n_r)
doubleprecision, intent(out)   :: sigma2(1:n_r),dust_accretion,dust_flux_o
doubleprecision, intent(out)   :: Diff(1:n_r),v_drift(1:n_r),v_05(1:n_r),flim(1:n_r),dust_flux(1:n_r)
doubleprecision, intent(out)   :: A(1:n_r),B(1:n_r),C(1:n_r),D(1:n_r)
integer, intent(in)            :: coagulation_method

doubleprecision :: T_05(1:n_r), sigma_g_05(1:n_r),cs_05(1:n_r),v_total ! interface values
doubleprecision :: u(1:n_r),u_old(1:n_r),g(1:n_r),h(1:n_r),k(1:n_r),l(1:n_r)
integer         :: n
doubleprecision :: St,v_dr,v_dr_old
doubleprecision :: get_St      ! function: calculates the stokes number

! first, calculate all the half-step arrays:
! temperature, gas surface density, dust diffusion coefficient and drift velocity
! calculate at grid-borders
do n = 2,n_r
    T_05(n) = 0.5d0*(T(n)+T(n-1))

    ! calculate speed of sound
    cs_05(n) = sqrt(k_b*T_05(n)/(mu*m_p))

    sigma_g_05(n) = 0.5d0*(sigma_g(n)+sigma_g(n-1))

    !St        =    grainsize*rho_s/sigma_g_05(n)*stokes_factor !ST
    St        =    get_ST(grainsize,T_05(n),sigma_g_05(n),x(n)) !ST: new method

    ! ----- dust diffusion
    if (dust_diffusion==1) then
        Diff(n) =   alpha(n) * k_b / (sqrt(Grav*M_star)*mu*m_p) * T(n) * x(n)**1.5d0 * &
                    & 1.d0/(1.d0 + St**2.d0)/Sc
    else
        Diff(n) = 0.d0
    endif

    ! ----- radial drift (headwind term)
    v_05(n) =    0.d0
    if ( (dust_radialdrift==1).and.(n>peak_position(1)) ) then !!! added ".and.(n>2)" condition because d/dr( Sig_g ) can diverge at the boundary
        if (grainsize>0.d0) then
            v_dr =   drift_fudge_factor/(St + 1.d0/St) * &
                    & min(cs_05(n), & !!!
                        & 1.d0/sqrt(Grav*M_star) * k_b/(mu*m_p) * sqrt(T_05(n)) * x05(n)**3.d0 / sigma_g_05(n) * &
                        & 1.d0/(x(n)-x(n-1)) * ( sigma_g(n) * sqrt(T(n)/x(n)**3.d0) - sigma_g(n-1) * sqrt(T(n-1)/x(n-1)**3.d0) ) &
                    & )
!            v_dr =   drift_fudge_factor/(St + 1.d0/St) * &
!                    & min(cs_05(n), & !!!
!                        & 2.d0/sqrt(Grav*M_star) * k_b/(mu*m_p) / (x(n-1)**(-1.5d0)+x(n)**(-1.5d0)) * &
!                        & 1.d0/(  sigma_g(n-1)/sqrt(T(n-1)*x(n-1)**3.d0) + sigma_g(n)/sqrt(T(n)*x(n)**3.d0)  ) * &
!                        & 1.d0/(x(n)-x(n-1)) * ( sigma_g(n) * sqrt(T(n)/x(n)**3.d0) - sigma_g(n-1) * sqrt(T(n-1)/x(n-1)**3.d0) ) &
!                    & )
            v_05(n) =    v_05(n) + v_dr
        else
            v_dr = 0.d0
        endif
        ! save this velocity
        if (n==peak_position(1)+1) then
            v_dr_old = v_dr
        endif
    endif

!!!EX
if (grainsize>0.d0) then
    v_05(1:peak_position(1)) = v_05(1:peak_position(1)) + v_dr_old
endif
!!!EX

    ! ----- gas-drag
    if ((dust_drag==1).and.(n>1).and.(n<=n_r)) then
        ! use only gas velocities smaller than sound speed
        if((abs(v_gas(n))<1d7)) then
            v_05(n)     = v_05(n) + v_gas(n)/(1d0+St**2.d0)
        else
            if ((n>1).and.(n<n_r)) then
                write(*,*) 'ERROR: TOO HIGH V_GAS AT i_r = ',n,':'
                write(*,*) '       v_gas  = ',v_gas(n)
                write(*,*) '       v_dust = ',v_05(n)
                write(*,*) '       cs     = ',cs_05(n)
                v_05(n)     = v_05(n) + v_gas(n)/(1d0+St**2.d0)
                !!!stop 58284
            endif
        endif
    endif
enddo


cs_05(1) = sqrt(k_b*T(1)/mu/m_p)
v_05(1)  = v_05(2)
v_05(n_r)  = v_05(n_r-1)



Diff(1)  = Diff(2)
Diff(n_r)    = 0d0
Diff(n_r-1)  = 0d0

if (coagulation_method==2) then!!!
    Diff(1:2)       = Diff(3)
    Diff(n_r-1:n_r) = Diff(n_r-2)
    v_05(1)  = v_05(2)
    v_05(n_r-1:n_r)  = 0d0
endif
!if ((coagulation_method==2).or.(coagulation_method==0)) then
!    Diff(n_r-2:n_r) = 0d0!!!
!    Diff(1:2)  = Diff(3)
!endif

!
! limit the radial velocity
!
forall (n=1:n_r,abs(v_05(n))>0.5d0*cs_05(n))
    ! next line means: v_05 = abs(0.5*c_s)*sign(v_05)
    v_05(n) = sign(0.5d0*cs_05(n),v_05(n))
end forall
!
! calculate the total radial velocity, that is
! advective velocity + mixing velocity, and from that
! the flux limiting factor
!
flim = 1.d0
do n = 2,n_r
    v_total = 0.5d0*(alpha(n-1)+alpha(n)) * cs_05(n)**2.d0/(sqrt(Grav*M_star/(0.5d0*(x(n-1)+x(n)))**3.d0))/Sc &
            & * (sigma_g(n-1)+sigma_g(n))/(sigma1(n-1)+sigma1(n)+1d-100) &
            & * 1.d0/(x(n)-x(n-1)) * (sigma1(n)/sigma_g(n) - sigma1(n-1)/sigma_g(n-1)) &
            & + v_05(n)
    flim(n) = min(1.d0,cs_05(n)/(  abs(v_total)+1d-100  ))
enddo
flim(1)   = 1.d0
flim(n_r) = 0.d0!!!

! ----- do the timestep
sigma2 = sigma1
u_old  = sigma1*x/grainmass
u      = u_old
g      = 1.d0
h      = sigma_g*x
K      = sigma_dot*x/grainmass
L      = 0.d0

! in: g2d ratio gradient = 0, outside: g2d ratio gradient = 0
call impl_donorcell_adv_diff_delta(n_r,x,Diff,v_05,g,h,K,L,flim,u,dt,1.d0,1.d0,0.d0,0.d0,0.d0,0.d0,&
&coagulation_method,A,B,C,D)


! update the surface density
sigma2 = u/x*grainmass
sigma1 = sigma2

! enforce the boundary condition
if (abs(sigma1(1)/sigma_g(1)/(sigma1(2)/sigma_g(2)))-1d0>1d-10) then
    sigma1(1)=sigma1(2)/sigma_g(2)*sigma_g(1)
    sigma2(1)=sigma1(1)
endif

do n=2,n_r-1
    v_drift(n) = 0.5d0*(v_05(n)+v_05(n+1))
enddo
Diff(1)    = Diff(2)
Diff(n_r)  = 0.d0
v_drift(1) = v_05(1)
v_drift(n_r) = v_05(n_r)

! get the exact implicit fluxes
do n = 2,n_r
    dust_flux(n) = u(n-1)*max(0.d0,v_05(n)) + u(n) * min(0.d0,v_05(n))&
                   & - flim(n) * 0.25d0*(Diff(n)+Diff(n-1)) * (h(n)+h(n-1)) *(g(n)/h(n)*u(n)-g(n-1)/h(n-1)*u(n-1)) / (x(n)-x(n-1))
enddo

! calculate dust accretion
n=0
dust_flux(1)   = dust_flux(n+2) + (  (u(n+1)-u_old(n+1))/dt - (k(n+1)+l(n+1)*u(n+1))) * (x05(n+2)-x05(n+1))
dust_accretion = -2d0*pi*dust_flux(1)*grainmass

! calculate dust outflow through the outer boundary
n=0
dust_flux_o    = dust_flux(n_r-n) + ( - (u(n_r-n)-u_old(n_r-n))/dt + (k(n_r-n)+l(n_r-n)*u(n_r-n))) * 2.d0*(x(n_r-n)-x05(n_r-n))
dust_flux_o    = 2d0*pi*dust_flux_o*grainmass

dust_flux = dust_flux*grainmass

end subroutine duststep_donorcell
! =============================================================================


! _____________________________________________________________________________
! Implicit donor cell advection-diffusion scheme with piecewise constant values
!
!     Perform one time step for the following PDE:
!
!        du    d /     \    d /              d  /     u  \ \
!        -- + -- | u v | - -- | h(x) Diff(x) -- |g(x)----| | = K + L u
!        dt   dx \     /   dx \              dx \    h(x)/ /
!
!     with boundary conditions
!
!         dgu/h |            |
!       p ----- |      + q u |       = r
!          dx   |x=xbc       |x=xbc
! INPUT:
!       n_x     = # of grid points
!       x()     = the grid
!       Diff()  = value of Diff @ cell center
!       v()     = the values for v @ interface (array(i) = value @ i-1/2)
!       g()     = the values for g(x)
!       h()     = the values for h(x)
!       K()     = the values for K(x)
!       L()     = the values for L(x)
!       flim()  = diffusion flux limiting factor at interfaces
!       u()     = the current values of u(x)
!       dt      = the time step
!
! OUTPUT:
!       u()     = the updated values of u(x) after timestep dt
!
! NEEDS:
!       subroutine  tridag(a,b,c,r,u,n)         to invert tridiagonal matrix
! _____________________________________________________________________________
subroutine impl_donorcell_adv_diff_delta(n_x,x,Diff,v,g,h,K,L,flim,u,dt,pl,pr,ql,qr,rl,rr,coagulation_method,A,B,C,D)
implicit none

integer,intent(in)             :: n_x,coagulation_method
doubleprecision, intent(in)    :: x(1:n_x),Diff(1:n_x),g(1:n_x),h(1:n_x),K(1:n_x),L(1:n_x),flim(1:n_x)
doubleprecision, intent(in)    :: v(1:n_x) ! array(n) = value @ n-1/2
doubleprecision, intent(inout) :: u(1:n_x)
doubleprecision, intent(in)    :: dt
doubleprecision, intent(out)   :: A(1:n_x),B(1:n_x),C(1:n_x),D(1:n_x)

doubleprecision :: rhs(1:n_x),u2(1:n_x)
doubleprecision :: D05(1:n_x),h05(1:n_x),vol
doubleprecision :: pl,pr,ql,qr,rl,rr
integer :: i

! ----- calculate the arrays at the interfaces
do i = 2,n_x
    D05(i) = flim(i) * 0.5d0 * (Diff(i-1) + Diff(i))
    h05(i) = 0.5d0 * (h(i-1) + h(i))
enddo

! ----- calculate the entries of the tridiagonal matrix
do i = 2,n_x-1
    vol = 0.5d0*(x(i+1)-x(i-1))
    A(i) = -dt/vol *  &
                & ( &
                    & + max(0.d0,v(i))  &
                    & + D05(i) * h05(i) * g(i-1) / (  (x(i)-x(i-1)) * h(i-1)  ) &
                & )
    B(i) = 1.d0 - dt*L(i) + dt/vol * &
                & ( &
                    & + max(0.d0,v(i+1))   &
                    & - min(0.d0,v(i))  &
                    & + D05(i+1) * h05(i+1) * g(i)   / (  (x(i+1)-x(i)) * h(i)    ) &
                    & + D05(i)   * h05(i)   * g(i)   / (  (x(i)-x(i-1)) * h(i)    ) &
                & )
    C(i) = dt/vol *  &
                & ( &
                    & + min(0.d0,v(i+1))  &
                    & - D05(i+1) * h05(i+1)  * g(i+1) / (  (x(i+1)-x(i)) * h(i+1)  ) &
                & )
    D(i) = -dt * K(i)
enddo

! ----- boundary Conditions
A(1)   = 0.d0
B(1)   = ql - pl*g(1) / (h(1)*(x(2)-x(1)))
C(1)   =      pl*g(2) / (h(2)*(x(2)-x(1)))
D(1)   = u(1)-rl

A(n_x) =    - pr*g(n_x-1) / (h(n_x-1)*(x(n_x)-x(n_x-1)))
B(n_x) = qr + pr*g(n_x)  / (h(n_x)*(x(n_x)-x(n_x-1)))
C(n_x) = 0.d0
D(n_x) = u(n_x)-rr
!
! if coagulation_method==2,
!  then we change the arrays and
!  give them back to the calling routine
! otherwise, we solve the equation
!
if (coagulation_method==2) then
    A = A/dt
    B = (B - 1d0)/dt
    C = C/dt
    D = D/dt
else
    !
    ! the old way
    !
    !rhs = u - D

    !
    ! the delta-way
    !
    do i = 2,n_x-1
        rhs(i) = u(i) - D(i) - (A(i)*u(i-1)+B(i)*u(i)+C(i)*u(i+1))
    enddo
    rhs(1)   = rl - (                B(1)*u(1)     + C(1)*u(2))
    rhs(n_x) = rr - (A(n_x)*u(n_x-1)+B(n_x)*u(n_x)            )

    ! solve for u2
    call tridag(A,B,C,rhs,u2,n_x)

    ! update u
    !u = u2   ! old way
    u = u+u2  ! delta way

endif

end subroutine impl_donorcell_adv_diff_delta
! =============================================================================


! _____________________________________________________________________________
! the tridag routine from Numerical Recipes in F77 rewritten to F95
!
! where:    a         =    lower diagonal entries
!            b        =    diagonal entries
!            c        =    upper diagonal entries
!            r        =    right hand side vector
!            u        =    result vector
!            n        =    size of the vectors
! _____________________________________________________________________________
subroutine tridag(a,b,c,r,u,n)
integer        :: n
doubleprecision    :: a(n),b(n),c(n),r(n),u(n)
integer, parameter :: NMAX=10000000
doubleprecision    :: bet,gam(NMAX)
integer :: j


if (b(1).eq.0.) stop 'tridag: rewrite equations'

bet = b(1)

u(1)=r(1)/bet

do j=2,n
    gam(j)    = c(j-1)/bet
    bet    = b(j)-a(j)*gam(j)
    if(bet.eq.0.) stop 'tridag failed'
    u(j)    = (r(j)-a(j)*u(j-1))/bet
enddo

do j=n-1,1,-1
    u(j)=u(j)-gam(j+1)*u(j+1)
enddo
end subroutine tridag
! =============================================================================

! ____________________________________________________________________________
! calculates stokes number
!
! INPUT:    a       = grain size
!           T       = temperature
!           sigma   = surface density
!           R       = radius
!
! RETURNS:  stokes number of paticle
!
! ____________________________________________________________________________
double precision function get_ST(a,T,sigma,R)
use constants,ONLY: pi,mu,m_p,sig_h2,k_b,Grav
use variables,ONLY: m_star,rho_s,stokes_factor
use switches, ONLY: stokes_regime
implicit none
doubleprecision,intent(in) :: a
doubleprecision,intent(in) :: T
doubleprecision,intent(in) :: sigma
doubleprecision,intent(in) :: R

doubleprecision :: lambda
doubleprecision :: cs
doubleprecision :: n
doubleprecision :: omega_k
doubleprecision :: H

cs      = sqrt(k_b*T/mu/m_p)
omega_k = sqrt(Grav*M_star/R**3d0)
H       = cs/omega_k
n       = sigma/(sqrt(2*pi)*H*m_p)
lambda  = 0.5d0/(sig_h2*n)

!
! check if we are in the Eppstein regime
! else, we are in the stokes regime
!
if ((a/lambda<9d0/4d0).or.(stokes_regime==0)) then
    get_ST = a*rho_s/sigma*stokes_factor
else
    !
    ! assume to be in the first stokes regime
    !
    get_ST = 2d0*pi/9d0  *a*a*rho_s/sigma/lambda
endif

end function get_ST
! ============================================================================
