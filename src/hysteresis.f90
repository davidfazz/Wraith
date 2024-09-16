module hysteresis
  use constants
  use functions 
  use potential
  implicit none 
  double precision, private, parameter :: bsat = 50.d0, db = 1.d0, binit = 0.d0
  double precision, private, parameter :: bdir(3) = [ 1.d0, 0.d0, 0.d0 ]
  double precision, private, parameter :: p_stp = 1.d-3, p_deps = 1.d-6
  integer, private, parameter :: istate = 6057, imaxopt = 5000
  character(len=20), private, parameter :: ifile = 'min.pts'

contains 

!...............................................................................

  subroutine hysteresis_()
    implicit none
    double precision :: xinit(C_N), x(C_N), g(C_N), m(3), babs, e
    integer :: ierr, i

! dg: load initial configuration from file. this is necessary as we want to 
!     start from a specific state.

    open(unit=11,file=trim(ifile),action='read',status='old',iostat=ierr)
    if (ierr .ne. 0) then 
      print *,''
      print *,'    hysteresis > configuration file is missing!'
      print *,''
      stop 
    endif 

    do i=1,istate 
      read(11,*,iostat=ierr) xinit
      if (ierr .ne. 0) then  
        print *,''
        print *,'    hysteresis > configuration file is missing!'
        print *,''
        stop 
      endif 
    enddo 
    close(11)

! dg: perform virgin curve. starting from magnetic field 'binit' and increase 
!     field incrementally by 'db' until 'bsat' is reached. the current magnetic 
!     field amplitude is stored in 'babs' and the field direction in 'bdir'. the 
!     current magnetic configuration is stored in 'x' and is optimized at each 
!     step by a simple gradient following routine.
    
    open(unit=19,file='vrgn.data',action='write',status='replace')

    babs = binit
    x = xinit
    do 
      call hysteresis_optimize(x)
      call potential_(x,e,g)
      m = potential_moment(x)

      write (19,'(5(F15.5))') babs, e, m

      babs = babs + db
      C_FIELD = babs*bdir 

      if (babs .gt. bsat) then 
        exit 
      endif
    enddo

    close(19)


! dg: perform hysteresis curve. starting from the configuration obtained in the 
!     virgin curve, the magnetic field runs from 'bsat' to '-bsat' and back to 
!     'bsat' in steps of 'db'. at each step, the configuration 'x' is optimized 
!     by gradient following.

    open(unit=19,file='hysteresis.data',action='write',status='replace')

! dg: go from 'bsat' to '-bsat'.

    babs = bsat
    do 
      call hysteresis_optimize(x) 
      call potential_(x,e,g)
      m = potential_moment(x) 

      write (19,'(5(F15.5))') babs, e, m

      babs = babs - db 
      C_FIELD = babs*bdir 

      if (babs .lt. -bsat) then 
        exit 
      endif 
    enddo 

! dg: go from '-bsat' to 'bsat'.

    do
      call hysteresis_optimize(x) 
      call potential_(x,e,g)
      m = potential_moment(x) 

      write (19,'(5(F15.5))') babs, e, m

      babs = babs + db 
      C_FIELD = babs*bdir 

      if (babs .gt. bsat) then 
        exit 
      endif 
    enddo 

    close(19)

    return 
  end subroutine hysteresis_

!...............................................................................
! simple routine to optimize a given configuration.

  subroutine hysteresis_optimize(x)
    implicit none 
    double precision :: x(C_N), g(C_N), e
    integer :: i


    do i=1,imaxopt
      call potential_(x,e,g)
      if ((f_rms(g) .lt. p_deps) .and. (i .gt. 20)) then 
        exit 
      endif 

      x = x - p_stp*g 
    enddo

    return 
  end subroutine hysteresis_optimize

!...............................................................................
  
end module hysteresis
