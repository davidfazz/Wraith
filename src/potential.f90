! potential.f90. 
! WRAITH - david gallina - gallina@uni-kassel.de - last update: 06/12/23.
!...............................................................................
! this module acts as an interface between the various optimization routines and 
! the different implemented potentials. 
!...............................................................................

module potential
  use constants 
  use functions
  use dpls1d
  use dpls2d
  use open1d
  use unsurf
  use chains
  implicit none
  save

contains

!...............................................................................
! initialization routines of each implemented potential.

  subroutine potential_init()
    implicit none

    if (dpls1d_log) then 
      call dpls1d_init()
    elseif (open1d_log) then
      call open1d_init()
    elseif (dpls2d_log) then
      call dpls2d_init()
    elseif (unsurf_log) then 
      call unsurf_init()
    elseif (chains_log) then 
      call chains_init()
    endif

    return
  end subroutine potential_init

!...............................................................................
! deinitialization routines of each implemented potential.

  subroutine potential_deinit()
    implicit none

    if (dpls1d_log) then
      call dpls1d_deinit()
    elseif (open1d_log) then
      call open1d_deinit()
    elseif (dpls2d_log) then
      call dpls2d_deinit()
    elseif (unsurf_log) then 
      call unsurf_deinit()
    endif

    return
  end subroutine potential_deinit

!...............................................................................
! routine to calculate the potential energy and the associated gradient.

  subroutine potential_(x,e,g)
    implicit none
    double precision, intent(in) :: x(C_N)
    double precision, intent(out) :: e, g(C_N)

    if (dpls1d_log) then
      call dpls1d_potential(x,e,g)
    elseif (open1d_log) then
      call open1d_potential(x,e,g)
    elseif (dpls2d_log) then
      call dpls2d_potential(x,e,g)
    elseif (unsurf_log) then 
      call unsurf_potential(x,e,g)
    elseif (chains_log) then 
      call chains_potential(x,e,g)
    endif

    return
  end subroutine potential_

!...............................................................................
! routine to calculate the hessian matrix.

  subroutine hessian_(x,h)
    implicit none
    double precision, intent(in) :: x(C_N)
    double precision, intent(out) :: h(C_N,C_N)

    if (dpls1d_log) then 
      call dpls1d_hessian(x,h)
    elseif (open1d_log) then 
      call open1d_hessian(x,h)
    elseif (dpls2d_log) then
      call dpls2d_hessian(x,h)
    elseif (unsurf_log) then 
      call unsurf_hessian(x,h)
    elseif (chains_log) then 
      call chains_hessian(x,h)
    endif

    return
  end subroutine hessian_

!...............................................................................
! function to rescale a configuration back into its defined interval.

  function potential_rescale(x)
    implicit none
    double precision, intent(in) :: x(C_N)
    double precision :: potential_rescale(C_N)
    
    if (dpls1d_log) then 
      potential_rescale = f_rescale_pol(x)
    elseif (open1d_log) then 
      potential_rescale = f_rescale_pol(x)
    elseif (dpls2d_log) then
      potential_rescale = f_rescale_sph(x)
    elseif (unsurf_log) then 
      potential_rescale = f_rescale_sph(x)
    elseif (chains_log) then 
      potential_rescale = f_rescale_pol(x)
    endif

    return
  end function potential_rescale

!...............................................................................
! function to calculate the distance between two configurations.

  function potential_distance(x,y)
    implicit none
    double precision, intent(in) :: x(C_N), y(C_N)
    double precision :: potential_distance

    if (dpls1d_log) then 
      potential_distance = f_polar_distance(x,y)
    elseif (open1d_log) then 
      potential_distance = f_polar_distance(x,y)
    elseif (dpls2d_log) then
      potential_distance = f_spherical_distance(x,y)
    elseif (unsurf_log) then 
      potential_distance = f_spherical_distance(x,y)
    elseif (chains_log) then 
      potential_distance = f_polar_distance(x,y)
    endif

    return
  end function potential_distance

!...............................................................................
! function to create a random configuration.

  function potential_random()
    implicit none
    double precision :: potential_random(C_N)

    if (dpls1d_log) then 
      potential_random = dpls1d_random()
    elseif (open1d_log) then 
      potential_random = open1d_random()
    elseif (dpls2d_log) then
      potential_random = dpls2d_random()
    elseif (unsurf_log) then 
      potential_random = unsurf_random()
    elseif (chains_log) then 
      potential_random = chains_random()
    endif

    return
  end function potential_random

!...............................................................................
! routine to calculate the logarithm of the vibrational entropy.

  subroutine potential_entropy(x,s,nev)
    implicit none
    double precision, intent(in) :: x(C_N)
    double precision, intent(out) :: s
    integer, intent(out) :: nev 

    if (dpls1d_log) then 
      call dpls1d_entropy(x,s,nev)
    elseif (open1d_log) then 
      call open1d_entropy(x,s,nev)
    elseif (dpls2d_log) then
      call dpls2d_entropy(x,s,nev)
    elseif (unsurf_log) then 
      call unsurf_entropy(x,s,nev)
    elseif (chains_log) then 
      call chains_entropy(x,s,nev)
    endif

    return
  end subroutine potential_entropy

!...............................................................................
! function to create the time-inverted configuration of an input configuration.

  function potential_reversal(x)
    implicit none 
    double precision, intent(in) :: x(C_N)
    double precision :: potential_reversal(C_N)

    if (dpls1d_log) then 
      potential_reversal = dpls1d_reversal(x)
    elseif (open1d_log) then 
      potential_reversal = open1d_reversal(x)
    elseif (dpls2d_log) then
      potential_reversal = dpls2d_reversal(x)
    elseif (unsurf_log) then 
      potential_reversal = unsurf_reversal(x)
    elseif (chains_log) then 
      potential_reversal = chains_reversal(x)
    endif

    return 
  end function potential_reversal

!...............................................................................
! function to calculate the magnetic moment of a configuration.

  function potential_moment(x)
    implicit none 
    double precision, intent(in) :: x(C_N)
    double precision :: potential_moment(3)

    if (dpls1d_log) then 
      potential_moment = f_moment_pol(x)
    elseif (open1d_log) then 
      potential_moment = f_moment_pol(x)
    elseif (dpls2d_log) then 
      potential_moment = f_moment_sph(x)
    elseif (unsurf_log) then 
      potential_moment = f_moment_sph(x)
    elseif (chains_log) then 
      potential_moment = f_moment_pol(x)
    endif 

    return 
  end function potential_moment 

!...............................................................................

end module potential
