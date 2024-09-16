! open1d.f90 
! WRAITH - david gallina - gallina@uni-kassel.de - last update: 07/12/23.
!...............................................................................
! this potential describes an ensemble of interacting magnetic nanoparticles, 
! where each mnp carries a magnetic moment that is described by a single polar 
! angle. the particles are placed on a two-dimensional surface with different 
! degrees of disorder and different kinds of underlying lattice structures.
! currently, the following magnetic effects are implemented:
!  - dipolar interactions between the mnps.
!  - zeeman interactions between the mnps and external magnetic fields.
! to be implemented:
!  - multipole interactions between the mnps.
!...............................................................................

module open1d
  use constants
  use functions
  use random
  implicit none 
    
! interaction arrays.

  double precision, allocatable, private :: xx(:,:), xy(:,:), yy(:,:)
  
! typa arrays.
  
  type(t_np), allocatable, private :: p_ens(:)

contains

!...............................................................................
! initialization routine:

  subroutine open1d_init()
    implicit none 
    double precision :: rvec(3), rabs, rab5
    integer :: i, j


! output of settings:

    print *,'    > open1d potential initializing..      '
    print *,'                                           '

! read potential parameters.

    call open1d_readparms()

    print *,'      > number of particles:',C_P

! calculate the interactions:

    allocate(xx(C_P,C_P),yy(C_P,C_P),xy(C_P,C_P))
    xx = 0.d0 
    xy = 0.d0 
    yy = 0.d0
    do i=1,C_P 
      do j=1,C_P 
        rvec = p_ens(i)%xyz - p_ens(j)%xyz
        rabs = NORM2(rvec)
        rab5 = rabs**5
        if (i .eq. j) then
          cycle 
        endif
        xx(i,j) = xx(i,j) + rvec(1)*rvec(1) / rab5
        xy(i,j) = xy(i,j) + rvec(1)*rvec(2) / rab5
        yy(i,j) = yy(i,j) + rvec(2)*rvec(2) / rab5 
      enddo
    enddo 

    return
  end subroutine open1d_init

!...............................................................................
! deinitialization.

  subroutine open1d_deinit()
    implicit none

    deallocate(p_ens)
    deallocate(xx)
    deallocate(xy)
    deallocate(yy)

    return 
  end subroutine open1d_deinit

!...............................................................................
! energy and gradient of input configuration x.

  subroutine open1d_potential(x,e,g)
    implicit none
    double precision, intent(in) :: x(C_N) 
    double precision, intent(out) :: g(C_N), e  
    double precision :: cosx(C_N), sinx(C_N), si(3), sj(3), di(3), term(3) 
    integer :: i, j

! initialize variables.

    sinx = SIN(x)
    cosx = COS(x)

! calculate dipolar energy and gradient.

    e = 0.d0
    g = 0.d0
    do i=1,C_P
      si = p_ens(i)%moment * [ cosx(i), sinx(i), 0.d0 ]
      di = p_ens(i)%moment * [-sinx(i), cosx(i), 0.d0 ]
      
      term = 0.d0
      do j=1,C_P
        sj = p_ens(j)%moment * [ cosx(j), sinx(j), 0.d0 ]

        term(1) = term(1) + yy(i,j)*sj(1) - 2.d0*xx(i,j)*sj(1) - 3.d0*xy(i,j)*sj(2)
        term(2) = term(2) + xx(i,j)*sj(2) - 2.d0*yy(i,j)*sj(2) - 3.d0*xy(i,j)*sj(1)
        term(3) = 0.d0
      enddo
      e = e + C_MU0/(8.d0*C_PI) * DOT_PRODUCT(si,term)
      g(i) = C_MU0/(4.d0*C_PI) * DOT_PRODUCT(di,term)

    enddo
! calculate zeeman energy and gradient.

    do i=1,C_P
      si = p_ens(i)%moment * [ cosx(i), sinx(i), 0.d0 ]
      di = p_ens(i)%moment * [-sinx(i), cosx(i), 0.d0 ]

      e = e - DOT_PRODUCT(si,C_FIELD)
      g(i) = g(i) - DOT_PRODUCT(di,C_FIELD)
    enddo

    return 
  end subroutine open1d_potential

!...............................................................................
! calculate hessian matrix of input configuration x.
  
  subroutine open1d_hessian(x,h)
    implicit none 
    double precision, intent(in) :: x(C_N)
    double precision, intent(out) :: h(C_N,C_N)
    double precision :: sinx(C_N), cosx(C_N), sj(3), dj(3), di(3), ti(3) 
    double precision :: oterm(3), dterm(3) 
    integer :: i, j

! initialize variables:

    sinx = SIN(x)
    cosx = COS(x)

! calculate hessian of dipole interaction;

    h = 0.d0
    do i=1,C_P
      di = p_ens(i)%moment * [ -sinx(i),  cosx(i), 0.d0 ]
      ti = p_ens(i)%moment * [ -cosx(i), -sinx(i), 0.d0 ]

      dterm = 0.d0
      do j=1,C_P
        sj = p_ens(j)%moment * [ cosx(j), sinx(j), 0.d0 ]
        dj = p_ens(j)%moment * [-sinx(j), cosx(j), 0.d0 ]

! off-diagonal terms:

        oterm(1) = yy(i,j)*dj(1) - 2.d0*xx(i,j)*dj(1) - 3.d0*xy(i,j)*dj(2)
        oterm(2) = xx(i,j)*dj(2) - 2.d0*yy(i,j)*dj(2) - 3.d0*xy(i,j)*dj(1)
        oterm(3) = xx(i,j)*dj(3) + yy(i,j)*dj(3)
        h(i,j) = h(i,j) + C_MU0/(4.d0*C_PI) * DOT_PRODUCT(di,oterm)

! diagonal terms:

        dterm(1) = dterm(1) + yy(i,j)*sj(1) - 2.d0*xx(i,j)*sj(1) - 3.d0*xy(i,j)*sj(2)
        dterm(2) = dterm(2) + xx(i,j)*sj(2) - 2.d0*yy(i,j)*sj(2) - 3.d0*xy(i,j)*sj(1)
        dterm(3) = 0.d0
      enddo

      h(i,i) = h(i,i) + C_MU0/(4.d0*C_PI) * DOT_PRODUCT(ti,dterm)
    enddo

! calculate hessian of zeeman interaction;

    do i=1,c_n
      ti = p_ens(i)%moment * [ -cosx(i), -sinx(i), 0.d0 ] 
      
      h(i,i) = h(i,i) - DOT_PRODUCT(ti,C_FIELD)
    enddo

    return
  end subroutine open1d_hessian

!...............................................................................
! rescaling of the vector x into the interval [0,2pi];

  function open1d_rescale(x)
    implicit none
    double precision, intent(in) :: x(C_N)
    double precision :: open1d_rescale(C_N)

    open1d_rescale = f_rescale_pol(x)

    return 
  end function open1d_rescale

!...............................................................................
! create random magnetic configuration:

  function open1d_random()
    implicit none
    double precision :: open1d_random(C_P) 
    integer :: i

    do i=1,C_P
      open1d_random(i) = random_uniform()*C_PI*2.d0
    enddo

    return 
  end function open1d_random

!...............................................................................
! calculate distance between vectors.

  function open1d_distance(x,y)
    implicit none 
    double precision, intent(in) :: x(C_P), y(C_P)
    double precision :: open1d_distance

    open1d_distance = f_polar_distance(x,y)

    return 
  end function open1d_distance

!...............................................................................
! find state having the opposite magnetization of x.

  function open1d_reversal(x) 
    implicit none 
    double precision, intent(in) :: x(C_N)
    double precision :: open1d_reversal(C_N)

    open1d_reversal = x + C_PI 
    open1d_reversal = open1d_rescale(open1d_reversal)

    return 
  end function open1d_reversal

!...............................................................................
! vibrational entropy of magnetic configuration x:

  subroutine open1d_entropy(x,s,nev)
    implicit none 
    double precision, intent(in) :: x(C_N) 
    double precision, intent(out) :: s 
    integer, intent(out) :: nev
    double precision :: h(C_N,C_N), evec(C_N,C_N), eval(C_N)
    integer :: i, nzero

    call open1d_hessian(x,h)
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
  end subroutine open1d_entropy

!...............................................................................
! read parameters:

  subroutine open1d_readparms()
    implicit none 
    integer   :: info, i, ierr

! magnetic field.

    call f_read_('magnetic_field','params.dat',C_FIELD,info)

! ensemble. 

    open(unit=11,file='ensemble.data',action='read',status='old')
    
    C_P = 0
    do
      read(11,*,iostat=ierr)
      if (ierr .ne. 0) then 
        exit 
      endif 
      C_P = C_P + 1
    enddo 
    allocate(p_ens(C_P))

    rewind(11)
    do i=1,C_P 
      read(11,*) p_ens(i)%xyz, p_ens(i)%radius, p_ens(i)%moment
    enddo 
    
    close(11)

    C_N = 1*C_P
    C_S = 2*C_P
    C_C = 3*C_P

    return 
  end subroutine open1d_readparms 

!...............................................................................

end module open1d
