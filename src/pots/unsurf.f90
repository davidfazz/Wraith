! unsurf.f90 
! WRAITH - david gallina - gallina@uni-kassel.de - last update: 08/12/23.
!...............................................................................

module unsurf
  use constants
  use functions
  use random
  implicit none 
    
! interaction arrays.

  double precision, allocatable, private :: xx(:,:), xy(:,:), yy(:,:)
  double precision, allocatable, private :: xz(:,:), yz(:,:), zz(:,:)
  
! typa arrays.
  
  type(t_np), allocatable, private :: p_ens(:)

! lattice properties.

  double precision, private :: p_dx, p_dy, p_angle, p_dsigma, p_hangle, p_zvary

! particle properties.

  double precision, private :: p_radius, p_height, p_rsigma 
  character(len=30), private :: p_material, p_ptype

! unit cell properties.

  integer, private :: p_nx, p_ny, p_npbcx, p_npbcy

! basis properties.

  double precision, private, allocatable :: p_dbasis(:)
  integer, private :: p_nbasis

! debugging.

  logical, private :: p_debug = .false., p_uneven = .false.
  
contains
  
!...............................................................................
! initialization routine:

  subroutine unsurf_init()
    implicit none 
    double precision, allocatable :: xgrad(:)
    double precision :: cpos(3), cvec(3), rvec(3), dabs, rabs, rab5, dvec(3) 
    double precision :: rada, radb, t, maxd, dx, dy, ycorrect, rcorrect
    double precision :: volume, moment
    integer :: inp, ibasis, ip, ix, iy, ibx, iby, i, j, k
    logical :: update, overlap

! read potential parameters.

    call unsurf_readparms()

! allocate the required memory.

    allocate(p_ens(C_P))

! output of settings:

    print *,'    > unsurf potential initializing..      '
    print *,'                                           '
    print *,'    >>> lattice settings..                 '
    print *,'      > x lattice constant: ',f_trimmed_(p_dx)
    print *,'      > y lattice constant: ',f_trimmed_(p_dy)
    print *,'      > lattice angle:      ',f_trimmed_(p_angle)
    print *,'      > lattice disorder:   ',f_trimmed_(p_dsigma)
    print *,'      > lattice basis:      ',f_trimmed_(p_nbasis)
    print *,'      > lattice z variance: ',f_trimmed_(p_zvary)
    print *,'    >>> particle settings..                '
    print *,'      > particle radius:    ',f_trimmed_(p_radius)
    print *,'      > particle height:    ',f_trimmed_(p_height)
    print *,'      > particle disorder:  ',f_trimmed_(p_rsigma)
    print *,'      > particle material:  ',f_trimmed_(p_material)
    print *,'      > particle type:      ',f_trimmed_(p_ptype)
    print *,'      > uneven surface:     ',p_uneven
    print *,'    >>> unit cell settings..               '
    print *,'      > x dimensions:       ',f_trimmed_(p_nx)
    print *,'      > y dimensions:       ',f_trimmed_(p_ny)
    print *,'      > x boundary cells:   ',f_trimmed_(p_npbcx)
    print *,'      > y boundary cells:   ',f_trimmed_(p_npbcy)
    print *,'                                           '

! create magnetic nanoparticles.

    dx = DBLE(p_nx)*p_dx 
    dy = DBLE(p_ny)*p_dy
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
    moment = volume*c_material(p_material,1)
    C_SCALE = C_MU0/(8.d0*C_PI)*moment*2/p_dx**3

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

! create xyz disorder for each np.

    do inp=1,C_P 
      rada = random_gauss(0.d0,p_dx*p_dsigma)
      radb = random_gauss(0.d0,p_dy*p_dsigma)
      t = random_uniform()*2.d0*C_PI
      
      rvec(1) = rada*COS(t)
      rvec(2) = radb*SIN(t)
      rvec(3) = random_gauss(0.d0,p_zvary)

      if (p_uneven) then 
        rvec(3) = rvec(3) + p_ens(inp)%radius 
      endif
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

    print *,'    >>> calculating periodic boundaries..'
    print *,''

! calculate the center of the unit cell:

    cpos(1) = dx + dy*COS(p_angle)
    cpos(2) = dy*SIN(p_angle)
    cpos(3) = 0.d0
    cpos = cpos / 2.d0
    maxd = 1.d8
    
! calculate the interactions:
    
    allocate(xx(C_P,C_P),yy(C_P,C_P),xy(C_P,C_P),xz(C_P,C_P),yz(C_P,C_P))
    allocate(zz(C_P,C_P))

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
              xz(i,j) = xz(i,j) + rvec(1)*rvec(3) / rab5
              yy(i,j) = yy(i,j) + rvec(2)*rvec(2) / rab5 
              yz(i,j) = yz(i,j) + rvec(2)*rvec(3) / rab5 
              zz(i,j) = zz(i,j) + rvec(3)*rvec(3) / rab5 
            else 
              xx(i,j) = xx(i,j) + rvec(1)*rvec(1) / rab5
              xy(i,j) = xy(i,j) + rvec(1)*rvec(2) / rab5
              xz(i,j) = xz(i,j) + rvec(1)*rvec(3) / rab5
              yy(i,j) = yy(i,j) + rvec(2)*rvec(2) / rab5 
              yz(i,j) = yz(i,j) + rvec(2)*rvec(3) / rab5 
              zz(i,j) = zz(i,j) + rvec(3)*rvec(3) / rab5 
            endif 
          enddo
        enddo
      enddo
    enddo 

! output:
    
    open(unit=11,file='ensemble.data',action='write',status='replace')
    do i=1,c_p 
      write(11,'(5(F15.5))') p_ens(i)%xyz,p_ens(i)%radius,p_ens(i)%moment
    enddo 
    close(11)

    print *,'    > unsurf potential initialized.'
    print *,''

    return
  end subroutine unsurf_init

!...............................................................................
! deinitialization:

  subroutine unsurf_deinit()
    implicit none

    deallocate(p_ens)
    deallocate(xx,xy,xz,yy,yz,zz)

    return 
  end subroutine unsurf_deinit

!...............................................................................
! energy and gradient of configuration x;

  subroutine unsurf_potential(x,e,g)
    implicit none
    double precision, intent(in) :: x(C_N)
    double precision, intent(out) :: e, g(C_N)
    double precision :: cosx(C_N), sinx(C_N), si(3), sj(3), dti(3), dpi(3) 
    double precision :: term(3) 
    integer :: i, j, i1, i2, j1, j2

! initialize variables:

    sinx = SIN(x)
    cosx = COS(x)

! calculate dipolar energy and gradient - the derivatives with respect to phi 
! are divided by sin(theta);
  
    e = 0.d0
    g = 0.d0
    do i=1,C_P 
      i1 = 2*i-1
      i2 = 2*i 
      si = p_ens(i)%moment * [ sinx(i1)*cosx(i2), sinx(i1)*sinx(i2), cosx(i1) ]
      dti = p_ens(i)%moment * [ cosx(i1)*cosx(i2), cosx(i1)*sinx(i2),-sinx(i1) ]
      dpi = p_ens(i)%moment * [-sinx(i2), cosx(i2), 0.d0 ]

      term = 0.d0 
      do j=1,C_P 
        j1 = 2*j-1 
        j2 = 2*j
        sj = p_ens(j)%moment * [ sinx(j1)*cosx(j2), sinx(j1)*sinx(j2), cosx(j1) ]
        term(1) = term(1) + (yy(i,j) - 2.d0*xx(i,j) + zz(i,j))*sj(1) &
                          - 3.d0*xy(i,j)*sj(2) - 3.d0*xz(i,j)*sj(3)
        term(2) = term(2) + (xx(i,j) - 2.d0*yy(i,j) + zz(i,j))*sj(2) &
                          - 3.d0*xy(i,j)*sj(1) - 3.d0*yz(i,j)*sj(3) 
        term(3) = term(3) + (xx(i,j) + yy(i,j) - 2.d0*zz(i,j))*sj(3) &
                          - 3.d0*xz(i,j)*sj(1) - 3.d0*yz(i,j)*sj(2) 
      enddo 
      e = e + C_MU0 / (8.d0*C_PI) * DOT_PRODUCT(si,term)
      g(i1) = g(i1) + C_MU0 / (4.d0*C_PI) * DOT_PRODUCT(dti,term)
      g(i2) = g(i2) + C_MU0 / (4.d0*C_PI) * DOT_PRODUCT(dpi,term)

      e = e - DOT_PRODUCT(C_FIELD,si)
      g(i1) = g(i1) - DOT_PRODUCT(C_FIELD,dti)
      g(i2) = g(i2) - DOT_PRODUCT(C_FIELD,dpi)
    enddo

    return 
  end subroutine unsurf_potential

!...............................................................................
! hessian of configuration x.

  subroutine unsurf_hessian(x,h)
    implicit none 
    double precision, intent(in) :: x(C_N)
    double precision, intent(inout) :: h(C_N,C_N)
    double precision :: sinx(C_N), cosx(C_N), g(C_N), si(3), dti(3), dpi(3)
    double precision :: sj(3), dtj(3), dpj(3), dt2(3), dp2(3), dtp(3), oterm(3)
    double precision :: dterm(3), e 
    integer :: i, j, i1, i2, j1, j2

! initialize variables:

    sinx = SIN(x)
    cosx = COS(x)

! calculate hessian:

    h = 0.d0
    do i=1,C_P
      i1 = 2*i-1 
      i2 = 2*i 
      dti = p_ens(i)%moment * [ cosx(i1)*cosx(i2), cosx(i1)*sinx(i2),-sinx(i1) ]
      dt2 = p_ens(i)%moment * [-sinx(i1)*cosx(i2),-sinx(i1)*sinx(i2),-cosx(i1) ]
      dpi = p_ens(i)%moment * [-sinx(i2), cosx(i2), 0.d0 ]
      dp2 = p_ens(i)%moment * [-sinx(i1)*cosx(i2),-sinx(i1)*sinx(i2), 0.d0 ]
      dtp = p_ens(i)%moment * [-cosx(i1)*sinx(i2), cosx(i1)*cosx(i2), 0.d0 ]

      dterm = 0.d0
      do j=1,C_P 
        j1 = 2*j-1
        j2 = 2*j
        sj = p_ens(j)%moment * [ sinx(j1)*cosx(j2), sinx(j1)*sinx(j2), cosx(j1) ]
        dtj = p_ens(j)%moment * [-cosx(j1)*sinx(j2), cosx(j1)*cosx(j2),-sinx(j1) ]
        dpj = p_ens(j)%moment * [-sinx(j2), cosx(j2), 0.d0 ]

! off-diagonal terms.
        
        oterm(1) = (yy(i,j) - 2.d0*xx(i,j) + zz(i,j))*dtj(1) &
                  - 3.d0*xy(i,j)*dtj(2) - 3.d0*xz(i,j)*dtj(3)
        oterm(2) = (xx(i,j) - 2.d0*yy(i,j) + zz(i,j))*dtj(2) &
                  - 3.d0*xy(i,j)*dtj(1) - 3.d0*yz(i,j)*dtj(3) 
        oterm(3) = (xx(i,j) + yy(i,j) - 2.d0*zz(i,j))*dtj(3) &
                  - 3.d0*xz(i,j)*dtj(1) - 3.d0*yz(i,j)*dtj(2)

        h(i1,j1) = h(i1,j1) + C_MU0/(4.d0*c_pi)*DOT_PRODUCT(dti,oterm)
        h(i2,j1) = h(i2,j1) + C_MU0/(4.d0*c_pi)*DOT_PRODUCT(dpi,oterm)

        oterm(1) = (yy(i,j) - 2.d0*xx(i,j) + zz(i,j))*dpj(1) &
                  - 3.d0*xy(i,j)*dtj(2) - 3.d0*xz(i,j)*dpj(3)
        oterm(2) = (xx(i,j) - 2.d0*yy(i,j) + zz(i,j))*dpj(2) &
                  - 3.d0*xy(i,j)*dtj(1) - 3.d0*yz(i,j)*dpj(3) 
        oterm(3) = (xx(i,j) + yy(i,j) - 2.d0*zz(i,j))*dpj(3) &
                  - 3.d0*xz(i,j)*dtj(1) - 3.d0*yz(i,j)*dpj(2)

        h(i1,j2) = h(i1,j2) + C_MU0/(4.d0*c_pi)*DOT_PRODUCT(dti,oterm)
        h(i2,j2) = h(i2,j2) + C_MU0/(4.d0*c_pi)*DOT_PRODUCT(dpi,oterm)

! diagonal terms.

        dterm(1) = dterm(1) + (yy(i,j) - 2.d0*xx(i,j) + zz(i,j))*sj(1) &
                            - 3.d0*xy(i,j)*sj(2) - 3.d0*xz(i,j)*sj(3) 
        dterm(2) = dterm(2) + (xx(i,j) - 2.d0*yy(i,j) + zz(i,j))*sj(2) &
                            - 3.d0*xy(i,j)*sj(1) - 3.d0*yz(i,j)*sj(3)
        dterm(3) = dterm(3) + (xx(i,j) + yy(i,j) - 2.d0*zz(i,j))*sj(3) &
                            - 3.d0*xz(i,j)*sj(1) - 3.d0*yz(i,j)*sj(2)
      enddo

      h(i1,i1) = h(i1,i1) + C_MU0/(4.d0*c_pi) * DOT_PRODUCT(dt2,dterm)
      h(i1,i2) = h(i1,i2) + C_MU0/(4.d0*c_pi) * DOT_PRODUCT(dtp,dterm)
      h(i2,i1) = h(i2,i1) + C_MU0/(4.d0*c_pi) * DOT_PRODUCT(dtp,dterm)
      h(i2,i2) = h(i2,i2) + C_MU0/(4.d0*c_pi) * DOT_PRODUCT(dp2,dterm)

      h(i1,i1) = h(i1,i1) - dot_product(dt2,C_FIELD)
      h(i1,i2) = h(i1,i2) - dot_product(dtp,C_FIELD)
      h(i2,i1) = h(i2,i1) - dot_product(dtp,C_FIELD)
      h(i2,i2) = h(i2,i2) - dot_product(dp2,C_FIELD)
    enddo

! correct hessian according to [Gra11]:

    call unsurf_potential(x,e,g)
    do i=1,C_P 
      i1 = 2*i-1
      i2 = 2*i
      h(i1,i1) = h(i1,i1)
      h(i1,i2) = h(i1,i2)/sinx(i1) - cosx(i1)/sinx(i1)*g(i2)
      h(i2,i1) = h(i2,i1)/sinx(i1) - cosx(i1)/sinx(i1)*g(i2)
      h(i2,i2) = h(i2,i2)/sinx(i1)/sinx(i1) + cosx(i1)/sinx(i1)*g(i1)
    enddo

    return
  end subroutine unsurf_hessian

!...............................................................................
! rescaling of each element of vector x into ([0,pi],[0,2pi]):

  subroutine unsurf_rescale(x)
    implicit none
    double precision, dimension(c_n), intent(inout) :: x
    integer                                         :: i

    x = f_rescale_sph(x)

    return 
  end subroutine unsurf_rescale

!...............................................................................
! create random configuration in spherical coordinates.

  function unsurf_random()
    implicit none
    double precision :: unsurf_random(C_N)
    integer :: i 

    do i=1,C_N,2 
      unsurf_random(i+0) = random_uniform()*C_PI 
      unsurf_random(i+1) = random_uniform()*C_PI*2.d0
    enddo
    
    return 
  end function unsurf_random

!...............................................................................
! find the time inverted state of x.

  function unsurf_reversal(x)
    implicit none 
    double precision, intent(in) :: x(C_N)
    double precision :: unsurf_reversal(C_N)
    integer :: i 

    do i=1,C_N,2 
      unsurf_reversal(i+0) = x(i+0) + C_PI 
      unsurf_reversal(i+1) = x(i+1)
    enddo 
    unsurf_reversal = f_rescale_sph(unsurf_reversal)

    return 
  end function unsurf_reversal

!...............................................................................
! vibrational entropy:

  subroutine unsurf_entropy(x,s,nev)
    implicit none 
    double precision, dimension(c_n), intent(in)  :: x 
    double precision, intent(out)                 :: s 
    integer, intent(out)                          :: nev
    double precision, dimension(c_n,c_n)          :: h 
    double precision, dimension(c_n,c_n)          :: evec
    double precision, dimension(c_n)              :: eval 
    integer                                       :: i
    integer                                       :: nzero

    call unsurf_hessian(x,h)
    call f_dsyev(h,evec,eval)

    nev = 0
    s = 0.d0
    nzero = 0
    do i=1,c_n
      if (abs(eval(i)) .lt. 1.d-5) then
        nzero = nzero + 1 
      elseif (eval(i) .lt. -1.d-5) then
        nev = nev + 1 
      else 
        s = s + log(eval(i))
      endif 
    enddo

    return 
  end subroutine unsurf_entropy

!...............................................................................
! read parameters:

  subroutine unsurf_readparms()
    implicit none 
    integer   :: info

! unit cell parameters.

    p_nx = 10
    call f_read_('unitcell_nx','params.dat',p_nx,info)
    
    p_ny = 10
    call f_read_('unitcell_ny','params.dat',p_ny,info)

    p_npbcx = 0
    call f_read_('unitcell_npbcx','params.dat',p_npbcx,info)

    p_npbcy = 0
    call f_read_('unitcell_npbcy','params.dat',p_npbcy,info)

! lattice parameters.

    p_dx = 10.d0
    call f_read_('lattice_dx','params.dat',p_dx,info)

    p_dy = 10.d0
    call f_read_('lattice_dy','params.dat',p_dy,info)

    p_angle = 90.d0
    call f_read_('lattice_angle','params.dat',p_angle,info)
    p_angle = p_angle / 180.d0 * C_PI

    p_dsigma = 0.d0
    call f_read_('lattice_sigma','params.dat',p_dsigma,info)

! basis parameters.

    p_nbasis = 1
    call f_read_('basis_nvecs','params.dat',p_nbasis,info)
    allocate(p_dbasis(3*p_nbasis))

    p_dbasis = 0.d0
    call f_read_('basis_dvecs','params.dat',p_dbasis,info)

! particle parameters.

    p_radius = 1.d0
    call f_read_('particle_radius','params.dat',p_radius,info)
    
    p_height = 0.d0
    call f_read_('particle_height','params.dat',p_height,info)
    
    p_rsigma = 0.d0
    call f_read_('particle_sigma','params.dat',p_rsigma,info)
    
    p_material = 'iron'
    call f_read_('particle_material','params.dat',p_material,info)

    p_ptype = 'sphere'
    call f_read_('particle_type','params.dat',p_ptype,info)
    
    p_uneven = .false.
    call f_read_('particle_uneven','params.dat',p_uneven,info)

    p_zvary = 0.d0
    call f_read_('particle_zvariance','params.dat',p_zvary,info)

! magnetic field parameters.

    C_FIELD = 0.d0
    call f_read_('magnetic_field','params.dat',C_FIELD,info)

! define some parameters.

    C_P = p_nx*p_ny*p_nbasis
    C_N = 2*C_P 
    C_S = 2*C_P 
    C_C = 3*C_P

    return 
  end subroutine unsurf_readparms 

!-------------------------------------------------------------------------------

end module unsurf
