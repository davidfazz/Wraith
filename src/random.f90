! random.f90 
! WRAITH - david gallina - last update: 04/12/2023
!...............................................................................
! marsaglia and tsang generator for random numbers based on "the ziggurat method 
! for generating random variables" [mar00]. 
!...............................................................................

module random 
  use constants 
  use functions 
  implicit none 

  integer, save, private :: jsr = 123456789

  private :: random_shr3
  public :: random_uniform, random_gauss, random_set
  
contains

!...............................................................................
! initialize ziggurat algorithm.

  subroutine random_set(iseed)
    implicit none 
    integer, intent(in) :: iseed 
    integer :: i 

    jsr = iseed

    return
  end subroutine random_set

!...............................................................................
! main function to find a random integer.

  function random_shr3()
    implicit none 
    integer :: random_shr3, jz

    jz = jsr 
    jsr = ieor(jsr,ishft(jsr,13))
    jsr = ieor(jsr,ishft(jsr,-17))
    jsr = ieor(jsr,ishft(jsr,5))
    random_shr3 = jz + jsr
    
    return
  end function random_shr3

!...............................................................................
! find uniform random number in [0,1).

  function random_uniform()
    implicit none 
    double precision :: random_uniform

    random_uniform = 0.5d0 + 0.2328306d-9*random_shr3()

    return
  end function random_uniform

!...............................................................................
! draw random number from a gaussian distribution.

  function random_gauss(am,sd)
    implicit none 
    double precision, intent(in) :: am, sd 
    double precision :: random_gauss, uni 
    integer :: i

    uni = 0.d0
    do i=1,12
      uni = uni + random_uniform()
    enddo 
    random_gauss = (uni - 6.d0)*sd + am

    return 
  end function random_gauss

!...............................................................................

end module random