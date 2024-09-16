module chains
  use constants
  use functions 
  use random
  implicit none 

! parameters.

  logical, private :: p_periodic 
  integer, private :: p_fixed(50), p_nfixed, p_nn 
  double precision, private :: exc_(6), dmi_(6), mae_

contains 

!...............................................................................

  subroutine chains_init()
    implicit none 

! setup interaction constants.

    C_N = 50
    C_P = 50

    p_periodic = .false.
    p_fixed = 0 
    p_nfixed = 0 
    p_nn = 4

    exc_ = [-85.2733d0,-42.292d0, 13.4593d0,-2.98483d0, 0.d0, 0.d0 ]
    dmi_ = [ 0.537629d0,-0.0456757d0, 0.0442406d0, 0.0161289d0, 0.d0, 0.d0 ]
    mae_ = 0.d0

    return 
  end subroutine chains_init

!...............................................................................

  subroutine chains_potential(x,e,g)
    implicit none 
    double precision :: x(C_N), e, g(C_N), cosx(C_N), sinx(C_N), si(3), sk(3)
    double precision :: di(3)
    integer :: i, j, k
    
!     e = 0.d0 
!     g = 0.d0

!     cosx = cos(x)
!     sinx = sin(x)

! ! loop through all atoms of the chain. then, loop through all positive (+x) and 
! ! negative (-x) neighbors.

!     do i=1,C_N 
!       si = [ cosx(i), 0.d0, sinx(i) ]
!       di = [-sinx(i), 0.d0, cosx(i) ]
  
! ! positive neighbors. 

!       do j=-p_nn,p_nn 
!         if (j .eq. 0) then 
!           cycle 
!         endif 

!         if (p_periodic) then 
!           if ((i+j .lt. 0) .or. (i+j .gt. C_N)) then 
!             cycle 
!           endif 
!           k = i + j 
!         elseif (.not. p_periodic) then
!           if (i+j .lt. 0) then 
            

!         endif

!         if (i+j .gt. C_N) then
!           if (p_periodic) then 
!             k = i + j - C_N
!           else 
!             cycle 
!           endif
!         else 
!           k = i + j
!         endif
!         sk = [ cosx(k), 0.d0, sinx(k) ]

! ! exchange.

!         e = e - exc_(j)*DOT_PRODUCT(si,sk)
!         g(i) = g(i) - exc_(j)*DOT_PRODUCT(di,sk)

! ! dmi.

!         e = e - dmi_(j)*(sk(3)*si(1) - sk(1)*si(3))
!         g(i) = g(i) - dmi_(j)*(sk(3)*di(1) - sk(1)*di(3))
!       enddo
      
! ! negative neighbors.

!       do j=1,nn 
!         if (i-j .lt. 1) then 
!           k = C_N + i - j
!         else 
!           k = i - j
!         endif
!         sk = [ sinx(k), 0.d0, cosx(k) ]

! ! exchange. 

!         e = e - exc_(j)*DOT_PRODUCT(si,sk)
!         g(i) = g(i) - exc_(j)*DOT_PRODUCT(di,sk)

! !dmi.

!         e = e + dmi_(j)*(sk(3)*si(1) - sk(1)*si(3))
!         g(i) = g(i) + dmi_(j)*(sk(3)*di(1) - sk(1)*di(3))
!       enddo
!     enddo
!     g(1) = 0.d0
!     e = e / 2.d0

    return 
  end subroutine chains_potential 

!...............................................................................

  subroutine chains_hessian(x,h)
    implicit none 
    double precision :: x(C_N), h(C_N,C_N), cosx(C_N), sinx(C_N), si(3), sk(3) 
    double precision :: di(3), ti(3), dk(3)
    integer :: i, j, k 

!     h = 0.d0 

!     cosx = cos(x)
!     sinx = sin(x)

! ! loop through all atoms of the chain. then, loop through all positive (+x) and 
! ! negative (-x) neighbors.

!     do i=1,C_N 
!       si = [ sinx(i), 0.d0, cosx(i) ]
!       di = [ cosx(i), 0.d0,-sinx(i) ]
!       ti = [-sinx(i), 0.d0,-cosx(i) ]
  
! ! positive neighbors. 

!       do j=1,nn 
!         if (i+j .gt. C_N) then 
!           k = i + j - C_N
!         else 
!           k = i + j
!         endif
!         sk = [ sinx(k), 0.d0, cosx(k) ]
!         dk = [ cosx(k), 0.d0,-sinx(k) ]

! ! offdiagonal elements.

!         h(i,k) = h(i,k) - exc_(j)*DOT_PRODUCT(di,dk)
!         h(i,k) = h(i,k) - dmi_(j)*(dk(3)*di(1) - dk(1)*di(3))

! ! diagonal elements.

!         h(i,i) = h(i,i) - exc_(j)*DOT_PRODUCT(ti,sk)
!         h(i,i) = h(i,i) - dmi_(j)*(sk(3)*ti(1) - sk(1)*ti(3))
!       enddo
      
! ! negative neighbors.

!       do j=1,nn 
!         if (i-j .lt. 1) then 
!           k = C_N + i - j
!         else 
!           k = i - j
!         endif
!         sk = [ sinx(k), 0.d0, cosx(k) ]
!         dk = [ cosx(k), 0.d0,-sinx(k) ]

! ! offdiagonal elements.

!         h(i,k) = h(i,k) - exc_(j)*DOT_PRODUCT(di,dk)
!         h(i,k) = h(i,k) + dmi_(j)*(dk(3)*di(1) - dk(1)*di(3))

! ! diagonal elements.

!         h(i,i) = h(i,i) - exc_(j)*DOT_PRODUCT(ti,sk)
!         h(i,i) = h(i,i) + dmi_(j)*(sk(3)*ti(1) - sk(1)*ti(3))
!       enddo
!     enddo

!     h(1,:) = 0.d0
!     h(:,1) = 0.d0

    return
  end subroutine chains_hessian

!...............................................................................
! rescaling of the vector x into the interval [0,2pi];

  function chains_rescale(x)
    implicit none
    double precision, intent(in) :: x(C_N)
    double precision :: chains_rescale(C_N)

    chains_rescale = f_rescale_pol(x)

    return 
  end function chains_rescale

!...............................................................................
! create random magnetic configuration:

  function chains_random()
    implicit none
    double precision :: chains_random(C_N) 
    integer :: i

    do i=1,C_N
      chains_random(i) = random_uniform()*C_PI*2.d0
    enddo
    chains_random(1) = 0.d0

    return 
  end function chains_random

!...............................................................................
! calculate distance between vectors.

  function chains_distance(x,y)
    implicit none 
    double precision, intent(in) :: x(C_N), y(C_N)
    double precision :: chains_distance

    chains_distance = f_polar_distance(x,y)

    return 
  end function chains_distance

!...............................................................................
! find state having the opposite magnetization of x.

  function chains_reversal(x) 
    implicit none 
    double precision, intent(in) :: x(C_N)
    double precision :: chains_reversal(C_N)

    chains_reversal = x + C_PI 
    chains_reversal = chains_rescale(chains_reversal)

    return 
  end function chains_reversal

!...............................................................................
! vibrational entropy of magnetic configuration x:

  subroutine chains_entropy(x,s,nev)
    implicit none 
    double precision, intent(in) :: x(C_N) 
    double precision, intent(out) :: s 
    integer, intent(out) :: nev
    double precision :: h(C_N,C_N), evec(C_N,C_N), eval(C_N)
    integer :: i, nzero

    call chains_hessian(x,h)
    call f_dsyev(h,evec,eval)

    nev = 0
    s = 0.d0
    nzero = 0
    do i=1,C_N
      if (ABS(eval(i)) .lt. 1.d-5) then
        nzero = nzero + 1 
      elseif (eval(i) .lt. -1.d-5) then
        nev = nev + 1 
      else 
        s = s + LOG(eval(i))
      endif 
    enddo

    return 
  end subroutine chains_entropy

!...............................................................................

end module chains
