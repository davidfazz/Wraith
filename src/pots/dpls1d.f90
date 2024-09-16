! dpls1d.f90 
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

module dpls1d
  use constants
  use functions
  use random
  implicit none 
    
! interaction arrays.

  double precision, allocatable, private :: xx(:,:), xy(:,:), yy(:,:)
  
! typa arrays.
  
  type(t_np), allocatable, private :: p_ens(:)

! lattice properties.

  double precision, private :: p_dx, p_dy, p_angle, p_dsigma, p_hangle

! particle properties.

  double precision, private :: p_radius, p_height, p_rsigma 
  character(len=30), private :: p_material, p_ptype

! unit cell properties.

  integer, private :: p_nx, p_ny, p_npbcx, p_npbcy

! basis properties.

  double precision, private, allocatable :: p_dbasis(:)
  integer, private :: p_nbasis

! debugging.

  logical, private :: p_debug = .false.

contains

!-------------------------------------------------------------------------------
! initialization routine:

  subroutine dpls1d_init()
    implicit none 
    double precision, allocatable :: xgrad(:)
    double precision :: cpos(3), cvec(3), rvec(3), dabs, rabs, rab5, dvec(3) 
    double precision :: rada, radb, t, maxd, dx, dy, ycorrect, rcorrect
    double precision :: volume, moment
    integer :: inp, ibasis, ip, ix, iy, ibx, iby, i, j, k
    logical :: update, overlap

! read potential parameters.

    call dpls1d_readparms()

! allocate the required memory.

    allocate(p_ens(C_P))

! output of settings:

    print *,'    > dpls1d potential initializing..      '
    print *,'                                           '
    print *,'    >>> lattice settings..                 '
    print *,'      > x lattice constant: ',f_trimmed_(p_dx)
    print *,'      > y lattice constant: ',f_trimmed_(p_dy)
    print *,'      > lattice angle:      ',f_trimmed_(p_angle)
    print *,'      > lattice disorder:   ',f_trimmed_(p_dsigma)
    print *,'      > lattice basis:      ',f_trimmed_(p_nbasis)
    print *,'    >>> particle settings..                '
    print *,'      > particle radius:    ',f_trimmed_(p_radius)
    print *,'      > particle height:    ',f_trimmed_(p_height)
    print *,'      > particle disorder:  ',f_trimmed_(p_rsigma)
    print *,'      > particle material:  ',f_trimmed_(p_material)
    print *,'      > particle type:      ',f_trimmed_(p_ptype)
    print *,'    >>> unit cell settings..               '
    print *,'      > x dimensions:       ',f_trimmed_(p_nx)
    print *,'      > y dimensions:       ',f_trimmed_(p_ny)
    print *,'      > x boundary cells:   ',f_trimmed_(p_npbcx)
    print *,'      > y boundary cells:   ',f_trimmed_(p_npbcy)
    print *,'                                           '

! create magnetic nanoparticles.

    dx = DBLE(p_nx)*p_dx 
    dy = DBLE(p_ny)*p_dy!*COS(p_angle)
    do inp=1,C_P
      p_ens(inp)%radius = random_gauss(p_radius,p_rsigma)
      p_ens(inp)%height = random_gauss(p_height,p_rsigma)
      if (trim(p_ptype) .eq. 'sphere') then 
        p_ens(inp)%volume = 4.d0/3.d0*C_PI*p_ens(inp)%radius**3
      elseif (trim(p_ptype) .eq. 'cylinder') then 
        p_ens(inp)%volume = C_PI*p_ens(inp)%radius**2*p_ens(inp)%height
      else 
        p_ens(inp)%volume = 4.d0/3.d0*C_PI*p_ens(inp)%radius**3
      endif
      p_ens(inp)%moment = p_ens(inp)%volume*c_material(p_material,1)*C_MUB
    enddo

! calculate scaling parameter.

    if (trim(p_ptype) .eq. 'sphere') then 
      volume = 4.d0/3.d0*C_PI*p_radius**3 
    elseif (trim(p_ptype) .eq. 'cylinder') then 
      volume = C_PI*p_radius**2*p_height 
    endif
    moment = volume*c_material(p_material,1)*C_MUB
    C_SCALE = C_MU0/(4.d0*C_PI)*moment**2/p_dx**3

! create lattice position for each np.

    inp = 0
    do ix=1,p_nx 
      do iy=1,p_ny 
        rvec(1) = (DBLE(ix) - 0.5d0)*p_dx + (DBLE(iy) - 0.5d0)*p_dy*COS(p_angle) 
        rvec(2) = (DBLE(iy) - 0.5d0)*p_dy*SIN(p_angle)
        rvec(3) = 0.d0 
        do ibasis=1,p_nbasis 
          inp = inp + 1
          rvec = rvec + p_dbasis(3*ibasis-2:3*ibasis)
          p_ens(inp)%xyz = rvec
        enddo
      enddo 
    enddo

! create xy disorder for each np.

    do inp=1,C_P 
      rada = random_gauss(0.d0,p_dx*p_dsigma)
      radb = random_gauss(0.d0,p_dy*p_dsigma)
      t = random_uniform()*2.d0*C_PI
      
      rvec(1) = rada*COS(t)
      rvec(2) = radb*SIN(t)
      rvec(3) = 0.d0 

      p_ens(inp)%xyz = p_ens(inp)%xyz + rvec
    enddo

! relaxe particles such that there is no overlap.

    p_hangle = ATAN(-1.d0/TAN(p_angle))
    if (p_hangle .lt. 0.d0) then 
      p_hangle = p_hangle + C_PI 
    endif

    allocate(xgrad(C_C))
    do i=1,10000
      xgrad = 0.d0
      overlap = .false.

! check for overlap between nps.

      do j=1,C_P 
        do k=1,C_P 
          if (j .eq. k) then 
            cycle 
          endif 
          dvec = p_ens(k)%xyz - p_ens(j)%xyz 
          dabs = NORM2(dvec)
          if (dabs .gt. p_ens(j)%radius + p_ens(k)%radius) then 
            cycle 
          else 
            overlap = .true.
            xgrad(3*j-2:3*j) = - dvec 
          endif
        enddo 
      enddo

! check for overlap with boundaries.

      do j=1,C_P 

! x boundaries.

        ycorrect = p_ens(j)%radius*ABS(SIN(p_hangle))
        rcorrect = ABS((COS(p_hangle) / 1.d0))
        if (p_ens(j)%xyz(1) - p_ens(j)%radius*rcorrect .lt. (p_ens(j)%xyz(2)+ycorrect)/TAN(p_angle)) then 
          xgrad(3*j-2) = 1.d0
          overlap = .true.
        elseif (p_ens(j)%xyz(1) + p_ens(j)%radius*rcorrect .gt. dx + (p_ens(j)%xyz(2)+ycorrect)/TAN(p_angle)) then
          xgrad(3*j-2) =-1.d0
          overlap = .true.
        endif

! y boundaries.

        if (p_ens(j)%xyz(2) - p_ens(j)%radius .lt. 0.d0) then 
          xgrad(3*j-1) = 1.d0
          overlap = .true.
        elseif (p_ens(j)%xyz(2) + p_ens(j)%radius .gt. dy) then 
          xgrad(3*j-1) =-1.d0
          overlap = .true.
        endif
      enddo

      if (overlap) then 
        do j=1,C_P 
          p_ens(j)%xyz = p_ens(j)%xyz + 1.d-2*xgrad(3*j-2:3*j)
        enddo
      else 
        exit 
      endif
    enddo
    deallocate(xgrad)

! we calculate the interaction matrices..

    print *,''
    print *,'    > periodic boundaries..'
    print *,''

! calculate the center of the unit cell:

    cpos(1) = dx + dy*COS(p_angle)
    cpos(2) = dy*SIN(p_angle)
    cpos(3) = 0.d0
    cpos = cpos / 2.d0
    maxd = 1.d8
    
! calculate the interactions:

    allocate(xx(C_P,C_P),yy(C_P,C_P),xy(C_P,C_P))
    xx = 0.d0 
    xy = 0.d0 
    yy = 0.d0
    do ibx=-p_npbcx,p_npbcx
      do iby=-p_npbcy,p_npbcy
        cvec(1) = dble(ibx)*dx + dble(iby)*dy*COS(p_angle)
        cvec(2) = dble(iby)*dy*SIN(p_angle)
        cvec(3) = 0.d0
        do i=1,C_P 
          do j=1,C_P 
            if (NORM2(p_ens(j)%xyz + cvec - cpos) .gt. maxd) then
              cycle 
            endif
            rvec = p_ens(i)%xyz - p_ens(j)%xyz - cvec
            rabs = NORM2(rvec)
            rab5 = rabs**5

            if ((ibx .eq. 0) .and. (iby .eq. 0)) then
              if (i .eq. j) then
                cycle 
              endif
              xx(i,j) = xx(i,j) + rvec(1)*rvec(1) / rab5
              xy(i,j) = xy(i,j) + rvec(1)*rvec(2) / rab5
              yy(i,j) = yy(i,j) + rvec(2)*rvec(2) / rab5 
            else 
              xx(i,j) = xx(i,j) + rvec(1)*rvec(1) / rab5
              xy(i,j) = xy(i,j) + rvec(1)*rvec(2) / rab5
              yy(i,j) = yy(i,j) + rvec(2)*rvec(2) / rab5 
            endif 
          enddo
        enddo
      enddo
    enddo 

! output:
    
    open(unit=11,file='ensemble.data',action='write',status='replace')
    do i=1,c_p 
      write(11,*) p_ens(i)%xyz,p_ens(i)%radius,p_ens(i)%moment
    enddo 
    close(11)

    return
  end subroutine dpls1d_init

!...............................................................................
! deinitialization.

  subroutine dpls1d_deinit()
    implicit none

    deallocate(p_ens)
    deallocate(xx)
    deallocate(xy)
    deallocate(yy)

    return 
  end subroutine dpls1d_deinit

!...............................................................................
! energy and gradient of input configuration x.

  subroutine dpls1d_potential(x,e,g)
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
  end subroutine dpls1d_potential

!...............................................................................
! calculate hessian matrix of input configuration x.
  
  subroutine dpls1d_hessian(x,h)
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
  end subroutine dpls1d_hessian

!...............................................................................
! rescaling of the vector x into the interval [0,2pi];

  function dpls1d_rescale(x)
    implicit none
    double precision, intent(in) :: x(C_N)
    double precision :: dpls1d_rescale(C_N)

    dpls1d_rescale = f_rescale_pol(x)

    return 
  end function dpls1d_rescale

!...............................................................................
! create random magnetic configuration:

  function dpls1d_random()
    implicit none
    double precision :: dpls1d_random(C_P) 
    integer :: i

    do i=1,C_P
      dpls1d_random(i) = random_uniform()*C_PI*2.d0
    enddo

    return 
  end function dpls1d_random

!...............................................................................
! calculate distance between vectors.

  function dpls1d_distance(x,y)
    implicit none 
    double precision, intent(in) :: x(C_P), y(C_P)
    double precision :: dpls1d_distance

    dpls1d_distance = f_polar_distance(x,y)

    return 
  end function dpls1d_distance

!...............................................................................
! find state having the opposite magnetization of x.

  function dpls1d_reversal(x) 
    implicit none 
    double precision, intent(in) :: x(C_N)
    double precision :: dpls1d_reversal(C_N)

    dpls1d_reversal = x + C_PI 
    dpls1d_reversal = dpls1d_rescale(dpls1d_reversal)

    return 
  end function dpls1d_reversal

!...............................................................................
! vibrational entropy of magnetic configuration x:

  subroutine dpls1d_entropy(x,s,nev)
    implicit none 
    double precision, intent(in) :: x(C_N) 
    double precision, intent(out) :: s 
    integer, intent(out) :: nev
    double precision :: h(C_N,C_N), evec(C_N,C_N), eval(C_N)
    integer :: i, nzero

    call dpls1d_hessian(x,h)
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
  end subroutine dpls1d_entropy

!...............................................................................
! read parameters:

  subroutine dpls1d_readparms()
    implicit none 
    integer   :: info

! unit cell parameters.

    call f_read_('unitcell_nx','params.dat',p_nx,info)
    call f_read_('unitcell_ny','params.dat',p_ny,info)

    call f_read_('unitcell_npbcx','params.dat',p_npbcx,info)
    call f_read_('unitcell_npbcy','params.dat',p_npbcy,info)

! lattice parameters.

    call f_read_('lattice_dx','params.dat',p_dx,info)
    call f_read_('lattice_dy','params.dat',p_dy,info)
    call f_read_('lattice_angle','params.dat',p_angle,info)
    p_angle = p_angle / 180.d0 * C_PI
    call f_read_('lattice_sigma','params.dat',p_dsigma,info)

! basis parameters.

    call f_read_('basis_nvecs','params.dat',p_nbasis,info)
    allocate(p_dbasis(3*p_nbasis))
    call f_read_('basis_dvecs','params.dat',p_dbasis,info)

! particle parameters.

    call f_read_('particle_radius','params.dat',p_radius,info)
    call f_read_('particle_height','params.dat',p_height,info)
    call f_read_('particle_sigma','params.dat',p_rsigma,info)
    call f_read_('particle_material','params.dat',p_material,info)
    call f_read_('particle_type','params.dat',p_ptype,info)

! magnetic field parameters.

    C_FIELD = 0.d0
    call f_read_('magnetic_field','params.dat',C_FIELD,info)

! define some parameters.

    C_P = p_nx*p_ny*p_nbasis
    C_N = C_P 
    C_S = 2*C_P 
    C_C = 3*C_P

    return 
  end subroutine dpls1d_readparms 

!...............................................................................

end module dpls1d
