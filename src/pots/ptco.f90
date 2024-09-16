! ptco.f90 
! WRAITH - david gallina - gallina@uni-kassel.de - last update: 07/12/23.
!...............................................................................
! this potential describes an ensemble of interacting co nanoislands located on 
! a two-dimensional pt surface. the islands are cylindrically shaped with 
! varying heights. each nanoisland is described by a single magnetic moment.
!...............................................................................

module ptco
  use constants
  use functions
  implicit none

! module parameters:

  double precision, private :: p_radius, p_height, p_psigma, p_lambda, p_lsigma 
  double precision, private :: p_nx, p_ny, p_bnx, p_bny

! material parameters:
! bulk and interface anisotropies and magnetization of cobalt.

  double precision, private, parameter :: co_mae_interface = 10.d0
  double precision, private, parameter :: co_mae_bulk = 2.684d0
  double precision, private, parameter :: co_magnet_bulk = 153.484d0

! interaction arrays:

  double precision, private, allocatable :: xx(:,:), xy(:,:), yy(:,:)

contains 

!...............................................................................

  subroutine ptco_init()
    implicit none 
    double precision :: rvec(3), cvec(3), angle, length, dx, dy, rabs, rab5
    integer :: i, j, imnp, ibx, iby

! read parameters.

    call ptco_readparms()

! calculate for each particle its height and width.

    do i=1,C_P 
      p_mnp(i)%radius = random_gauss(p_radius,p_radius*p_psigma) 
      p_mnp(i)%height = random_gauss(p_height,p_height*p_psigma)
      p_mnp(i)%area = p_mnp(i)%radius*p_mnp(i)%radius
      p_mnp(i)%volume = p_mnp(i)%area*p_mnp(i)%height

      p_mnp(i)%imae = p_mnp(i)%area*p_mae_interface
      p_mnp(i)%vmae = p_mnp(i)%volume*p_mae*volume
      p_mnp(i)%moment = p_mnp(i)%volume*p_magnet_bulk
    enddo

! calculate for each particle its position.

    imnp = 0
    do i=1,p_nx 
      do j=1,p_ny
        rvec = [ DBLE(i)*p_lambda, DBLE(j)*p_lambda, 0.d0 ]
        rvec = rpos - [ p_lambda, p_lambda, 0.d0]/2.d0

        imnp = imnp + 1
        p_mnp(imnp)%xyz = rvec
      enddo 
    enddo
    
! create small amounts of structural disorder.

    do i=1,C_P 
      angle = random_number()*2.d0*C_PI
      length = random_gauss(p_lambda,p_lambda*p_lsigma)
      rvec = length * [ COS(angle), SIN(angle), 0.d0 ]
      p_mnp(i)%xyz = p_mnp(i)%xyz + rvec
    enddo

! calculate interaction constants.

    allocate(xx(C_P,C_P),xy(C_P,C_P),yy(C_P,C_P))
    xx = 0.d0 
    xy = 0.d0 
    yy = 0.d0

    dx = DBLE(p_nx)*p_lambda
    dy = DBLE(p_ny)*p_lambda
    do ibx=-p_bnx,p_bnx
      do iby=-p_bny,p_bny
        cvec(1) = dble(ibx)*dx 
        cvec(2) = dble(iby)*dy        
        cvec(3) = 0.d0
        
        do i=1,C_P 
          do j=1,C_P
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

! create output.

    open(unit=11,file='ensemble.data',action='write',status='replace')
    do i=1,C_P 
      write(11,'(5(F12.5))') p_mnp(i)%xyz,p_mnp(i)%radius,p_mnp(i)%moment
    enddo 
    close(11)

    return 
  end subroutine ptco_init 

!...............................................................................

  subroutine ptco_deinit()
    implicit none 

    return 
  end subroutine ptco_deinit 

!...............................................................................
! cartesian coordinates.
!...............................................................................

  subroutine ptco_energy(x,e)
    implicit none 
    double precision :: x(C_C), e, edip, eint, evol, ezee, si(3), sj(3)
    integer :: i, j, i1, i3, j1, j3

! initialize total energy.

    e = 0.d0
    g = 0.d0 

! calculate individual energy contributions.

    edip = 0.d0
    eint = 0.d0
    evol = 0.d0
    ezee = 0.d0
    do i=1,C_P 
      i1 = 3*i - 2
      i3 = 3*i 
      si = x(i1:i3) 

! dipole-dipole interaction.

      term = 0.d0
      do j=1,C_P 
        j1 = 3*j - 2 
        j3 = 3*j 
        sj = x(j1:j3)
        term(1) = term(1) + yy(i,j)*sj(1) - 2.d0*xx(i,j)*sj(1) - 3.d0*xy(i,j)*sj(2)
        term(2) = term(2) + xx(i,j)*sj(2) - 2.d0*yy(i,j)*sj(2) - 3.d0*xy(i,j)*sj(1)
        term(3) = term(3) + xx(i,j)*sj(3) + yy(i,j)*sj(3)
      enddo 
      edip = edip + C_MU0/(8.d0*C_PI)*DOT_PRODUCT(si,sj)

! zeeman interaction.

      ezee = ezee - DOT_PRODUCT(si,C_FIELD)

! interface anisotropy.

      eint = eint - p_ens(i)%imae*si(3)*si(3)

! volume anisotropy.

      evol = evol - p_ens(i)%vmae*si(3)*si(3)
    enddo 

    return 
  end subroutine ptco_energy 

!...............................................................................

  subroutine ptco_gradient(x,g)
    implicit none 
    double precision :: x(C_C), g(C_C), si(3), sj(3) 
    integer :: i, j, i1, i3, j1, j3

! calculate gradient.

    g = 0.d0
    do i=1,C_P 
      i1 = 3*i - 2
      i2 = 3*i - 1
      i3 = 3*i 
      si = x(i1:i3)

! dipole-dipole interaction.

      term = 0.d0
      do j=1,C_P 
        j1 = 3*j - 2
        j3 = 3*j 
        sj = x(j1:j3)
        term(1) = term(1) + yy(i,j)*sj(1) - 2.d0*xx(i,j)*sj(1) - 3.d0*xy(i,j)*sj(2)
        term(2) = term(2) + xx(i,j)*sj(2) - 2.d0*yy(i,j)*sj(2) - 3.d0*xy(i,j)*sj(1)
        term(3) = term(3) + xx(i,j)*sj(3) + yy(i,j)*sj(3)
        g(i1) = g(i1) + C_MU0/(4.d0*C_PI)*term(1) 
        g(i2) = g(i2) + C_MU0/(4.d0*C_PI)*term(2) 
        g(i3) = g(i3) + C_MU0/(4.d0*C_PI)*term(3)
      enddo 

! zeeman interaction.

      g(i1) = g(i1) - C_FIELD(1) 
      g(i2) = g(i2) - C_FIELD(2)
      g(i3) = g(i3) - C_FIELD(3)

! interface anisotropy.

      g(i1) = g(i1) 
      g(i2) = g(i2)
      g(i3) = g(i3) - 2.d0*p_ens(i)%imae*si(3)

! volume anisotropy.

      g(i1) = g(i1) 
      g(i2) = g(i2) 
      g(i3) = g(i3) - 2.d0*p_ens(i)%vmae*si(3)
    enddo 

    return 
  end subroutine ptco_gradient

!...............................................................................
! potential routine that calculates the energy as well as the gradient.

  subroutine ptco_potential(x,e,g)
    implicit none 
    double precision :: x(C_C), e, g(C_C), edip, eint, evol, ezee, si(3), sj(3)
    integer :: i, j, i1, i2, i3, j1, j3

! calculate energy and gradient.

    call ptco_energy(x,e)
    call ptco_gradient(x,g)

! project gradient on sphere.

    g = f_project_sph(x,g)

    return 
  end subroutine ptco_potential 

!...............................................................................
! routine to calculate the euclidean hessian.

  subroutine ptco_hessian_eucl(x,h)
    implicit none 
    double precision :: x(C_C), h(C_C,C_C)
    integer :: i, j, i1, i2, i3, j1, j2, j3 

! calculate hessian matrix.

    h = 0.d0
    do i=1,C_P 
      i1 = 3*i - 2 
      i2 = 3*i - 1
      i3 = 3*i 

! dipole-dipole interaction.

      do j=1,C_P 
        j1 = 3*j - 2
        j2 = 3*j - 1
        j3 = 3*j

        if (j .eq. i) then 
          cycle 
        endif

        h(i1,j1) = h(i1,j1) + C_MU0/(4.d0*C_PI)*(yy(i,j) - 2.d0*xx(i,j))
        h(i1,j2) = h(i1,j2) 
        h(i1,j3) = h(i1,j3) 
        h(i2,j1) = h(i2,j1)  
        h(i2,j2) = h(i2,j2) + C_MU0/(4.d0*C_PI)*(xx(i,j) - 2.d0*yy(i,j)) 
        h(i2,j3) = h(i2,j3)  
        h(i3,j1) = h(i3,j1)  
        h(i3,j2) = h(i3,j2) 
        h(i3,j3) = h(i3,j3) + C_MU0/(4.d0*C_PI)*(xx(i,j) + yy(i,j))
      enddo 

! interface anisotropy.

      h(i3,i3) = h(i3,i3) - 2.d0*p_ens(i)%imae 

! volume anisotropy.

      h(i3,i3) = h(i3,i3) - 2.d0*p_ens(i)%vmae 

    enddo 

    return 
  end subroutine ptco_hessian_eucl 

!...............................................................................
! routine to calculate the riemannian hessian.

  subroutine ptco_hessian(x,h)
    implicit none 
    double precision :: x(C_C), h(C_C,C_C), g(C_C), heuc(C_C,C_C), geuc(C_C)
    double precision :: unity(3,3), xxt(3,3), hpar(3,3)
    integer :: i1, i3, j1, j3

! unit 3x3 matrix.

    unity = 0.d0
    unity(1,1) = 1.d0 
    unity(2,2) = 1.d0 
    unity(3,3) = 1.d0

! calculate euclidean hessian.

    call ptco_hessian_eucl(x,heuc)

! calculate euclidean gradient.

    call ptco_gradient(x,geuc)

! calculate hessian projections.

    do i=1,C_P 
      i1 = 3*i - 2
      i3 = 3*i

      xxt = f_projector(x(i1:i3))
      do j=1,C_P
        j1 = 3*j - 2
        j3 = 3*j 

        hpar = heuc(i1:i3,j1:j3)
        if (i .eq. j) then 
          hpar = hpar - DOT_PRODUCT(geuc(i1:i3),x(i1:i3))*unity
          h(i1:i3,i1:i3) = MATMUL(unity - xxt,hpar)
        elseif (i .ne. j) then 
          h(i1:i3,j1:j3) = MATMUL(xxt,hpar)
        endif
      enddo 
    enddo

    return 
  end subroutine ptco_hessian 

!...............................................................................
! calculate the finite differences hessian along the direction u;
!   > see eq.(10.41) in [Bou22].

  subroutine ptco_hessu(x,u,h)
    implicit none 
    double precision :: x(C_C), u(C_C), h(C_C), gu(C_C), g(C_C), t, e

! set step length - step length ||t*u|| should be 2^(-14);

    t = 1.d0 / (f_norm(u) * 16384.d0)

! calculation of the gradient at x and at Rx(t*u);

    call ptco_potential(f_retract(x,t,u),e,gu)
    call ptco_potential(x,e,g)

! transport gradient back to TxM;

    gu = f_transport(x,gu)

! calculate finite differences hessian;

    h = (gu - g) / t

    return 
  end subroutine ptco_hessu

!...............................................................................
! spherical coordinates.
!...............................................................................


!...............................................................................

  subroutine ptco_readparms()
    implicit none 
    integer   :: info

! unit cell parameters.

    call f_read_('unitcell_nx','params.dat',p_nx,info)
    call f_read_('unitcell_ny','params.dat',p_ny,info)

    call f_read_('unitcell_npbcx','params.dat',p_bnx,info)
    call f_read_('unitcell_npbcy','params.dat',p_bny,info)

! lattice parameters.

    call f_read_('lattice_lambda','params.dat',p_lambda,info)
    call f_read_('lattice_sigma','params.dat',p_lsigma,info)

! particle parameters.

    call f_read_('particle_radius','params.dat',p_radius,info)
    call f_read_('particle_height','params.dat',p_height,info)
    call f_read_('particle_sigma','params.dat',p_psigma,info)

! define some parameters.

    C_P = p_nx*p_ny
    C_N = C_P 
    C_S = 2*C_P 
    C_C = 3*C_P

    return 
  end subroutine ptco_readparms 

!...............................................................................

end module ptco