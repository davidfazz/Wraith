module follow
  use constants
  use functions
  use potential 
  use dneb
  use lbfgs
  implicit none 
  double precision, parameter, private :: bdir(3) = [ 1.d0, 0.d0, 0.d0 ] !/ SQRT(2.d0)
  double precision, parameter, private :: bint = 0.d0
  double precision, parameter, private :: bend = 50.d0
  double precision, parameter, private :: binc = 1.d-1

  integer, parameter, private :: istate = 5
  character(len=10), parameter, private :: input = 'min.pts'

contains 

!...............................................................................

  subroutine follow_()
    implicit none
    double precision :: xflw(C_N), xtwo(C_N), g(C_N), e, m(3), mabs, xsp(C_N)
    double precision :: eone, etwo, xone(C_N), esp, x(C_N), psi, mzero(3), mplus(3)
    double precision :: mone, mtwo
    integer :: i, ierr, ipos, ix, iy, j
    logical :: exist, conv
    character(len=20) :: file, newfile


    C_FIELD = 6.d0 * [ 1.d0, 0.d0, 0.d0] 
    C_FIELD = 6.d0 * [ 1.d0, 1.d0, 0.d0] / SQRT(2.d0)

! get initial configuration to follow.
    
    xflw = 0.d0
    call follow_optim(xflw)
    call potential_(xflw,eone,g)      
    m = 0.d0
    do j=1,C_N 
      m(1) = m(1) + COS(xflw(j))
      m(2) = m(2) + SIN(xflw(j))
    enddo 
    m = m / DBLE(C_N)
    print *,SQRT(DOT_PRODUCT(g,g)),eone,eone/49.0653
    
    xflw = C_PI/4.d0
    call follow_optim(xflw)
    call potential_(xflw,eone,g)      
    m = 0.d0
    do j=1,C_N 
      m(1) = m(1) + COS(xflw(j))
      m(2) = m(2) + SIN(xflw(j))
    enddo 
    m = m / DBLE(C_N)
    print *,SQRT(DOT_PRODUCT(g,g)),eone,eone/49.0653

    xflw = C_PI/2.d0
    call follow_optim(xflw)
    call potential_(xflw,eone,g)      
    m = 0.d0
    do j=1,C_N 
      m(1) = m(1) + COS(xflw(j))
      m(2) = m(2) + SIN(xflw(j))
    enddo 
    m = m / DBLE(C_N)
    print *,SQRT(DOT_PRODUCT(g,g)),eone,eone/49.0653

    xflw =7.d0*C_PI/4.d0
    call follow_optim(xflw)
    call potential_(xflw,eone,g)      
    m = 0.d0
    do j=1,C_N 
      m(1) = m(1) + COS(xflw(j))
      m(2) = m(2) + SIN(xflw(j))
    enddo 
    m = m / DBLE(C_N)
    print *,SQRT(DOT_PRODUCT(g,g)),eone,eone/49.0653

    stop



    do i=100,-100,-1
      C_FIELD = DBLE(i)*bdir*binc
      call follow_optim(xflw)
      call potential_(xflw,eone,g)
      
      m = 0.d0
      do j=1,C_N 
        m(1) = m(1) + COS(xflw(j))
        m(2) = m(2) + SIN(xflw(j))
      enddo 
      m = m / DBLE(C_N)

      print *,dble(i)*binc*bdir(1),DOT_PRODUCT(m,bdir),eone
    enddo
    stop       



    open(unit=11,file=trim(input),action='read',status='old',iostat=ierr)
    if (ierr .ne. 0) then 
      print *,'    follow > input file is missing or corrupted.'
      print *,''
      stop 
    endif 

    do i=1,4!20
      read(11,*) xflw 
      call potential_(xflw,e,g)
      print *,i,e

      if (i .eq. 2) then 
        xone = xflw 
      endif 

      if (i .eq. 4) then 
        xtwo = xflw 
      endif
    enddo
    close(11)

    ! open(unit=11,file='energies.data',action='write',status='replace')
    ! open(unit=12,file='configs.data',action='write',status='replace')
    
    ! xflw = 0.d0
    ! do i=1,2500
    !   xflw = potential_random()
    !   xflw(1) = 0.d0
    !   call lbfgs_(C_N,5,xflw,1.d-6,conv,e,1000,j,.FALSE.)
    !   call potential_(xflw,e,g)
    !   print *,i,e,dot_product(g,g)
    !   write(11,*) i,e 
    !   write(12,*) xflw
    ! enddo 
    ! close(11)
    ! close(12)
    ! stop

    ! open(unit=15,file='suscep.data',action='write',status='replace')
    ! do i=0,200
    !   psi = dble(i)*2.d0*C_PI / 200.d0

    !   ipos = 0
    !   do ix=1,16
    !     do iy=1,16
    !       ipos = ipos + 1
    !       if ((mod(ix,2) .eq. 0) .and. (mod(iy,2) .eq. 0)) then 
    !         x(ipos) =-psi 
    !       elseif ((mod(ix,2) .eq. 1) .and. (mod(iy,2) .eq. 0)) then 
    !         x(ipos) = psi
    !       elseif ((mod(ix,2) .eq. 0) .and. (mod(iy,2) .eq. 1)) then 
    !         x(ipos) = psi - C_PI
    !       elseif ((mod(ix,2) .eq. 1) .and. (mod(iy,2) .eq. 1)) then 
    !         x(ipos) = C_PI - psi
    !       endif
    !     enddo 
    !   enddo 

    !   !print *,x
    !   C_FIELD = 0.d0 * bdir
    !   mzero = potential_moment(x)
    !   mone = DOT_PRODUCT(bdir,mzero)
    !   call potential_(x,e,g)
    !   !print *,'A:',psi,e,mzero(1)

    !   C_FIELD = 1.d-2 * bdir
    !   call follow_optim(x)
    !   mplus = potential_moment(x)
    !   mtwo = DOT_PRODUCT(bdir,mplus)
    !   call potential_(x,e,g)
    !   !print *,'B:',psi,e,mplus(1),(mplus(1) - mzero(1))/ 1.d-2

    !   write(15,'(2(F15.5))') psi,(mtwo-mone)/ 1.d-2
    ! enddo 
    ! close(15)
    ! stop

    do i=0,70
      C_FIELD = DBLE(i)*binc*bdir
      call follow_optim(xone)
      call follow_optim(xtwo)

      if (mod(i,5) .eq. 0) then
        !call dneb_(xone,xtwo,xsp,15,spring=1.d0,printpath='full')
      
        call dneb_(xone,xtwo,xsp,15,spring=1.d0)
        call potential_(xone,eone,g)
        call potential_(xtwo,etwo,g)
        call potential_(xsp,esp,g)
        !call potential_(xtwo,eout,g)
        print *,dble(i)*binc*bdir(1),eone,esp,etwo
      endif
    enddo
    stop    



! follow state as magnetic field is increased.

    open(unit=11,file='follow.data',action='write',status='replace')
    C_FIELD = bdir*bint 
    do
      call follow_optim(xflw)
      call potential_(xflw,e,g)
      m = potential_moment(xflw)
      mabs = DOT_PRODUCT(C_FIELD,bdir)

      write(11,*) mabs,e,DOT_PRODUCT(m,bdir)
      print *,mabs,e,DOT_PRODUCT(m,bdir)
      
      C_FIELD = C_FIELD + binc*bdir
      if (mabs .gt. bend) then 
        exit 
      elseif (mabs .lt. 0.d0) then 
        exit
      endif
    enddo

    do
      call follow_optim(xflw)
      call potential_(xflw,e,g)
      m = potential_moment(xflw)
      mabs = DOT_PRODUCT(C_FIELD,bdir)

      write(11,*) mabs,e,DOT_PRODUCT(m,bdir)
      print *,mabs,e,DOT_PRODUCT(m,bdir)
      C_FIELD = C_FIELD - binc*bdir
      if (mabs .gt. bend) then 
        exit 
      elseif (mabs .lt. 0.d0) then 
        exit
      endif
    enddo  
    close(11)

    return
  end subroutine follow_

!...............................................................................

  subroutine follow_optim(x)
    implicit none 
    double precision :: x(C_N), g(C_N), e
    integer :: i 
    
    do i=1,10000
      call potential_(x,e,g)
      !if ((i .gt. 25) .and. (f_rms(g) .lt. 1.d-6)) then 
      !  exit 
      !endif
      x = x - 1.d-3*g 
    enddo 

    if (i .eq. 1000) then 
      print *,'NOT CONVERGED'
    endif

    return 
  end subroutine follow_optim

!...............................................................................

end module follow

