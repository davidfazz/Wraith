! functions.f90 
! WRATIH - david gallina - last update 04/12/23 
!...............................................................................
! module contains a large number of auxiliary functions.
! mathematical routines.

! basic routines.
!  1) f_outer: outer product between two 3-vectors.
!  2) f_rms: root mean square of a n-vector.
!  3) f_euclidean_distance: euclidean distance between two n-vectors.
!  4) f_riemann_distance: riemannian distance between two n-vectors.
!  5) f_spherical_distance: riemannian distance in spherical coordinates.
!  6) f_polar_distance: riemannian distance in polar coordinates.
!  7) f_convert_s2c: coordinate transformation from spherical to cartesian.
!  8) f_convert_c2s: coordinate transformation from cartesian to spherical.
!  9) f_convert_p2c: coordinate transformation from polar to cartesian.
! 10) f_convert_c2p: coordinate transformation from cartesian to polar.
! 11) f_angle: calculates the angle between two 3-vectors.
! 12) f_dysev: calculates the eigenvectors/eigenvalues of a symmetric matrix.
! 13)* f_dgeev: calculates the eigenvectors/eigenvalues of a nonsymmetric matrix.
! 14) f_mean: calculates the mean of a set of values.
! 15) f_variance: calculates the variance of a set of values.
! 16) f_stdev: calculates the standard deviation of a set of values.
! 17)* f_skew: calculates the skew of a set of values.
! 
!...............................................................................

module functions
  use constants 
  implicit none 

!...............................................................................
! interfaces.

! interface for f_read_ function.

  interface f_read_
    procedure f_read_int
    procedure f_read_dbl
    procedure f_read_chr
    procedure f_read_log
    procedure f_read_vec
  end interface f_read_

! interface for f_read_store_ function.

  interface f_read_store_
    procedure f_read_store_int
    procedure f_read_store_dbl
    procedure f_read_store_chr
    procedure f_read_store_log
    procedure f_read_store_vec
  end interface f_read_store_  

! interface for f_trimmed_ function.

  interface f_trimmed_
    procedure f_trimmed_int
    procedure f_trimmed_dbl
    procedure f_trimmed_chr
    procedure f_trimmed_log
    procedure f_trimmed_vec
  end interface f_trimmed_ 

contains 

!...............................................................................
! mathematical routines
!...............................................................................
! currently imple

  pure function f_outer(x,y)
    double precision, intent(in) :: x(3), y(3)
    double precision :: f_outer(3)

    f_outer(1) = x(2)*y(3) - x(3)*y(2)
    f_outer(2) = x(3)*y(1) - x(1)*y(3)
    f_outer(3) = x(1)*y(2) - x(2)*y(1)

    return 
  end function f_outer 

!...............................................................................
! root mean square of an input vector x.

  pure function f_rms(x)
    implicit none
    double precision, intent(in) :: x(:)
    double precision :: f_rms 

    f_rms = DOT_PRODUCT(x,x)
    f_rms = SQRT(f_rms/DBLE(SIZE(x)))

    return
  end function f_rms

!...............................................................................

  pure function f_euclidean_distance(x,y)
    implicit none 
    double precision, intent(in) :: x(:), y(:) 
    double precision :: f_euclidean_distance 

    f_euclidean_distance = SQRT(DOT_PRODUCT(x-y,x-y))

    return 
  end function f_euclidean_distance

!...............................................................................
! riemannian distance.

  pure function f_riemann_distance(x,y)
    implicit none 
    double precision, intent(in) :: x(:), y(:)
    double precision :: f_riemann_distance, angle 
    integer :: i 

    f_riemann_distance = 0.d0
    do i=1,SIZE(x),3
      angle = DOT_PRODUCT(x(3*i-2:3*i),y(3*i-2:3*i))
      angle = MIN(angle,1.d0)
      angle = MAX(angle,-1.d0)
      angle = ACOS(angle)
      f_riemann_distance = f_riemann_distance + angle*angle 
    enddo 
    f_riemann_distance = SQRT(f_riemann_distance)

    return 
  end function f_riemann_distance

!...............................................................................
! riemannian distance in spherical coordinates.

  pure function f_spherical_distance(x,y)
    implicit none 
    double precision, intent(in) :: x(:), y(:)
    double precision :: f_spherical_distance, angle
    integer :: i, i1, i2

    f_spherical_distance = 0.d0
    do i=1,SIZE(x),2 
      i1 = i+0
      i2 = i+1
      angle = SIN(x(i1))*COS(x(i2))*SIN(y(i1))*COS(y(i2)) &
            + SIN(x(i1))*SIN(x(i2))*SIN(y(i1))*SIN(y(i2)) &
            + COS(x(i1))*COS(y(i1))
      angle = ACOS(MAX(MIN(angle,1.d0),-1.d0))
      f_spherical_distance = f_spherical_distance + angle*angle
    enddo 
    f_spherical_distance = SQRT(f_spherical_distance)

    return 
  end function f_spherical_distance

!...............................................................................
! riemannian distance in polar coordinates.

  pure function f_polar_distance(x,y)
    implicit none 
    double precision, intent(in) :: x(:), y(:)
    double precision :: f_polar_distance, angle 
    integer :: i 

    f_polar_distance = 0.d0 
    do i=1,SIZE(x)
      angle = MIN(ABS(x(i) - y(i)), 2.d0*C_PI - ABS(x(i) - y(i)))
      f_polar_distance = f_polar_distance + angle*angle
    enddo 
    f_polar_distance = SQRT(f_polar_distance)

    return 
  end function f_polar_distance

!...............................................................................
! transforms unit vector from spherical to euclidean coordinates.

  pure function f_convert_s2c(x)
    implicit none
    double precision, intent(in) :: x(2) 
    double precision :: f_convert_s2c(3)

    f_convert_s2c(1) = SIN(x(1))*COS(x(2)) 
    f_convert_s2c(2) = SIN(x(1))*SIN(x(2))
    f_convert_s2c(3) = COS(x(1))

    return 
  end function f_convert_s2c 

!...............................................................................
! convert unit vector from cartesian coordinates to spherical coordinates.

  pure function f_convert_c2s(x)
    implicit none
    double precision, intent(in) :: x(3) 
    double precision :: f_convert_c2s(2)
    
    f_convert_c2s(1) = ACOS(x(3))
    f_convert_c2s(2) = ATAN2(x(2),x(1))

    return 
  end function f_convert_c2s

!...............................................................................
! convert unit vector from polar coordinates to cartesian coordinates.

  pure function f_convert_p2c(x)
    implicit none
    double precision, intent(in) :: x
    double precision :: f_convert_p2c(2) 

    f_convert_p2c(1) = COS(x)
    f_convert_p2c(2) = SIN(x)

    return 
  end function f_convert_p2c 

!...............................................................................
! convert cartesian coordinates to polar coordinates:

  pure function f_convert_c2p(x)
    implicit none
    double precision, intent(in) :: x(2)
    double precision :: f_convert_c2p 

    f_convert_c2p = atan2(x(2),x(1))

    return 
  end function f_convert_c2p

!...............................................................................
! angle between two 3-cartesian vectors.

  pure function f_angle(x,y)
    implicit none
    double precision, intent(in) :: x(3), y(3)
    double precision :: f_angle
    
    f_angle = ACOS(DOT_PRODUCT(x/NORM2(x),y/NORM2(y)))

    return 
  end function f_angle

!...............................................................................
! function to project vectors on sphere.

  function f_project_sph(x,g)
    implicit none 
    double precision :: x(C_C), g(C_C), f_project_sph(C_C)
    integer :: i, i1, i3 

    do i=1,C_P
      i1 = 3*i - 2
      i3 = 3*i 

      f_project_sph(i1:i3) = g(i1:i3) - DOT_PRODUCT(g(i1:i3),x(i1:i3))*x(i1:i3)
    enddo 

    return 
  end function f_project_sph

!...............................................................................
! create projector matrix.

  function f_projector(x)
    implicit none 
    double precision :: x(3), f_projector(3,3)

    f_projector(1,1) = x(1)*x(1) 
    f_projector(1,2) = x(1)*x(2)
    f_projector(1,3) = x(1)*x(3)
    f_projector(2,1) = x(2)*x(1)
    f_projector(2,2) = x(2)*x(2)
    f_projector(2,3) = x(2)*x(3)
    f_projector(3,1) = x(3)*x(1)
    f_projector(3,2) = x(3)*x(2)
    f_projector(3,3) = x(3)*x(3)

    return 
  end function f_projector 

!...............................................................................
! norm on sphere.

  function f_norm_sph(x)
    implicit none 
    double precision :: x(C_C), f_norm_sph(C_C)
    integer :: i, i1, i3 

    do i=1,C_P 
      i1 = 3*i - 2
      i3 = 3*i 

      f_norm_sph(i1:i3) = x(i1:i3) / NORM2(x(i1:i3))
    enddo 
    return 
  end function f_norm_sph

!...............................................................................
! simplified retraction function.

  function f_retract(x,t,u)
    implicit none 
    double precision :: x(C_C), t, u(C_C), f_retract(C_C)
    integer :: i, i1, i3 

    f_retract = x + t*u 
    f_retract = f_norm_sph(f_retract)

    return 
  end function f_retract

!...............................................................................
! simplified transportation function.

  function f_transport(x,v)
    implicit none 
    double precision :: x(C_C), v(C_C), f_transport(C_C)
    integer :: i, i1, i3 

    do i=1,C_P 
      i1 = 3*i - 2
      i3 = 3*i 

      f_transport(i1:i3) = v(i1:i3) - DOT_PRODUCT(x(i1:i3),v(i1:i3))*x(i1:i3)
    enddo 

    return 
  end function f_transport

!...............................................................................
! eigenvalues and eigenvectors of a symmetric matrix hsym:

  subroutine f_dsyev(hsym,evec,eval)
    implicit none
    double precision, intent(in) :: hsym(:,:)
    double precision, intent(out) :: evec(:,:), eval(:)
    double precision, allocatable :: copy(:,:), work(:)
    integer :: lwork, info, n, i

! define work arrays.

    n = SIZE(hsym,1)
    lwork = 33*n

! allocate memory.

    allocate(copy(n,n))
    allocate(work(lwork))
    copy = hsym
    work = 0.d0

  ! call dsyev routine.

    call dsyev('V','U',n,copy,n,eval,work,lwork,info)
    if (info .ne. 0) then
      print *,''
      print *,'    f_dsyev > diagonalisation failed!'
      print *,''
      return
    endif

  ! store eigenvalues and eigenvectors:  

    do i=1,n 
      evec(i,:) = copy(:,i)
    enddo

    deallocate(copy)
    deallocate(work)

    return
  end subroutine f_dsyev  

!...............................................................................
! average of a set of values.

  pure function f_mean(x)
    implicit none 
    double precision, dimension(:), intent(in) :: x
    double precision :: f_mean 
    integer :: i 

    f_mean = 0.d0 
    do i=1,SIZE(x)
      f_mean = f_mean + x(i)
    enddo
    f_mean = f_mean / dble(size(x))

    return 
  end function f_mean

!...............................................................................
! variance of a set of values.

  pure function f_variance(x,mean)
    implicit none 
    double precision, intent(in) :: x(:)
    double precision, intent(in), optional :: mean
    double precision :: f_variance, dmean 
    integer :: i

    if (present(mean)) then 
      dmean = mean 
    else 
      dmean = f_mean(x)
    endif 

    f_variance = 0.d0 
    do i=1,SIZE(x)
      f_variance = f_variance + (x(i) - dmean)*(x(i) - dmean)
    enddo 
    f_variance = f_variance / DBLE(SIZE(x))

    return 
  end function f_variance

!...............................................................................
! standard deviation of a set of values.

  pure function f_stdev(x,mean)
    implicit none 
    double precision, intent(in) :: x(:)
    double precision, intent(in), optional :: mean 
    double precision :: f_stdev, dmean

    if (present(mean)) then 
      dmean = mean 
    else 
      dmean = f_mean(x) 
    endif 

    f_stdev = SQRT(f_variance(x,dmean))

    return 
  end function f_stdev

!...............................................................................
! skew of a set of values.

  pure function f_skew(x)
    implicit none 
    double precision, intent(in) :: x(:) 
    double precision :: f_skew 

    f_skew = x(1)

    return 
  end function f_skew

!...............................................................................
! rescaling of a configuration in cartesian coordinates.

  pure function f_rescale_car(x)
    implicit none 
    double precision, intent(in) :: x(C_C)
    double precision :: f_rescale_car(C_C)
    integer :: i, i1, i3

    f_rescale_car = x
    do i=1,C_P
      i1 = 3*i-2 
      i3 = 3*i 
      f_rescale_car(i1:i3) = f_rescale_car(i1:i3) / NORM2(f_rescale_car(i1:i3)) 
    enddo

    return 
  end function f_rescale_car 

!...............................................................................
! rescaling of a configuration in spherical coordinates.

  pure function f_rescale_sph(x)
    implicit none 
    double precision, intent(in) :: x(C_N)
    double precision :: f_rescale_sph(C_N)
    integer :: i 

    f_rescale_sph = x
    do i=1,C_N,2
      do while (f_rescale_sph(i) .lt. 0.d0)
        f_rescale_sph(i) = f_rescale_sph(i) + 2.d0*C_PI 
      enddo 
      do while (f_rescale_sph(i) .gt. 2.d0*C_PI)
        f_rescale_sph(i) = f_rescale_sph(i) - 2.d0*C_PI 
      enddo 

      if (f_rescale_sph(i) .gt. C_PI) then 
        f_rescale_sph(i) = 2.d0*C_PI - f_rescale_sph(i)
        f_rescale_sph(i+1) = f_rescale_sph(i+1) + C_PI
      endif 

      do while (f_rescale_sph(i+1) .gt. 2.d0*C_PI)
        f_rescale_sph(i+1) = f_rescale_sph(i+1) - 2.d0*C_PI 
      enddo 
      do while (f_rescale_sph(i+1) .lt. 0.d0)
        f_rescale_sph(i+1) = f_rescale_sph(i+1) + 2.d0*C_PI
      enddo
    enddo

    return 
  end function f_rescale_sph

!...............................................................................

  pure function f_rescale_pol(x)
    implicit none 
    double precision, intent(in) :: x(C_P)
    double precision :: f_rescale_pol(C_P)
    integer :: i

    f_rescale_pol = x
    do i=1,C_P 
      do while (f_rescale_pol(i) .lt. 0.d0)
        f_rescale_pol(i) = f_rescale_pol(i) + 2.d0*C_PI 
      enddo 
      do while (f_rescale_pol(i) .gt. 2.d0*C_PI)
        f_rescale_pol(i) = f_rescale_pol(i) - 2.d0*C_PI 
      enddo
    enddo

    return 
  end function f_rescale_pol

!...............................................................................
! unit cell moment in cartesian coordinates.

  pure function f_moment_car(x)
    implicit none 
    double precision, intent(in) :: x(C_C)
    double precision :: f_moment_car(3)
    integer :: i 

    f_moment_car = 0.d0
    do i=1,C_C,3 
      f_moment_car = f_moment_car + x(i:i+2)
    enddo
    f_moment_car = f_moment_car / DBLE(C_P)

    return 
  end function f_moment_car 

!...............................................................................
! unit cell moment in spherical coordinates.

  pure function f_moment_sph(x)
    implicit none 
    double precision, intent(in) :: x(C_N)
    double precision :: f_moment_sph(3)
    integer :: i 

    f_moment_sph = 0.d0 
    do i=1,C_N,2 
      f_moment_sph(1) = f_moment_sph(1) + SIN(x(i))*COS(x(i+1))
      f_moment_sph(2) = f_moment_sph(2) + SIN(x(i))*SIN(x(i+1))
      f_moment_sph(3) = f_moment_sph(3) + COS(x(i))
    enddo 
    f_moment_sph = f_moment_sph / DBLE(C_P)

    return 
  end function f_moment_sph

!...............................................................................
! unit cell moment in polar coordinates.

  pure function f_moment_pol(x)
    implicit none 
    double precision, intent(in) :: x(C_P)
    double precision :: f_moment_pol(3)
    integer :: i 

    f_moment_pol = 0.d0 
    do i=1,C_P 
      f_moment_pol(1) = f_moment_pol(1) + COS(x(i))
      f_moment_pol(2) = f_moment_pol(2) + SIN(x(i))
      f_moment_pol(3) = f_moment_pol(3) + 0.d0
    enddo
    f_moment_pol = f_moment_pol / DBLE(C_P)

    return 
  end function f_moment_pol

!-------------------------------------------------------------------------------
! STRING routines
!-------------------------------------------------------------------------------
! removes all white spaces from a string:

  subroutine f_stripspaces(chr)
    implicit none
    character(len=*), intent(inout) :: chr
    integer                         :: stringlen
    integer                         :: last
    integer                         :: actual

    stringlen = len(chr)
    last = 1
    actual = 1
    do while (actual .lt. stringlen)
      if (chr(last:last) .eq. ' ') then
        actual = actual + 1
        chr(last:last) = chr(actual:actual)
        chr(actual:actual) = ' '
      else
        last = last + 1
        if (actual .lt. last) then
          actual = last
        endif
      endif
    enddo
    
    return
  end subroutine f_stripspaces

!-------------------------------------------------------------------------------
! read integers:

  subroutine f_read_int(chartag,charfile,value,info)
    implicit none
    character(len=*), intent(in)  :: chartag
    character(len=*), intent(in)  :: charfile
    integer, intent(out)          :: value
    integer, intent(out)          :: info
    character(len=25)             :: tag
    integer                       :: ioerr

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_int > ',trim(charfile),' could not be openend. no value is'
      print *,'            returned.'
      return
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif
      
      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,value
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)

    return
  end subroutine f_read_int

!-------------------------------------------------------------------------------
! read doubles:

  subroutine f_read_dbl(chartag,charfile,value,info)
    implicit none
    character(len=*), intent(in)                :: chartag
    character(len=*), intent(in)                :: charfile
    double precision                            :: value
    integer, intent(out)                        :: info
    character(len=25)                           :: tag
    integer                                     :: ioerr

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_dp > ',trim(charfile),' could not be openend. no value is '
      print *,'           returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,value
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)

    return
  end subroutine f_read_dbl

!-------------------------------------------------------------------------------
! read characters:

  subroutine f_read_chr(chartag,charfile,value,info)
    implicit none
    character(len=*), intent(in)  :: chartag
    character(len=*), intent(in)  :: charfile
    character(len=*), intent(out) :: value
    integer, intent(out)          :: info
    character(len=25)             :: tag
    integer                       :: ioerr

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_char > ',trim(charfile),' could not be openend. no value'
      print *,'             is returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,value
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    value = trim(value)
    close(12)

    return
  end subroutine f_read_chr

!-------------------------------------------------------------------------------
! read vectors:

  subroutine f_read_vec(chartag,charfile,value,info)
    implicit none
    character(len=*), intent(in)                :: chartag
    character(len=*), intent(in)                :: charfile
    double precision, dimension(:), intent(out) :: value
    integer, intent(out)                        :: info
    character(len=25)                           :: tag
    integer                                     :: ioerr

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_vec > ',trim(charfile),' could not be openend. no value is '
      print *,'             returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,value
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)

    return
  end subroutine f_read_vec

!-------------------------------------------------------------------------------
! read logicals:

  subroutine f_read_log(chartag,charfile,value,info)
    implicit none
      character(len=*), intent(in)    :: chartag
      character(len=*), intent(in)    :: charfile
      logical, intent(out)            :: value
      integer, intent(out)            :: info
      character(len=30)               :: tag
      integer                         :: ioerr

    value = .false.

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_char > ',trim(charfile),' could not be openend. no value'
      print *,'             is returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        info = 0
        value = .true.
        exit
      endif
    enddo
    close(12)

    return
  end subroutine f_read_log

!-------------------------------------------------------------------------------
! read integers:

  subroutine f_read_store_int(chartag,charfile,values,default,info)
    implicit none
    character(len=*), intent(in)  :: chartag
    character(len=*), intent(in)  :: charfile
    integer, intent(out)          :: values
    integer, intent(in)           :: default
    integer, intent(out)          :: info
    character(len=25)             :: tag
    integer                       :: ioerr

    values = 0
    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_int > ',trim(charfile),' could not be openend. no value is'
      print *,'            returned.'
      return
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif
      
      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,values
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)

    if (info .ne. 0) then
      values = default 
    endif

    return
  end subroutine f_read_store_int

!-------------------------------------------------------------------------------
! read doubles:

  subroutine f_read_store_dbl(chartag,charfile,values,default,info)
    implicit none
      character(len=*), intent(in)                :: chartag
      character(len=*), intent(in)                :: charfile
      double precision, intent(out)               :: values
      double precision, intent(in)                :: default
      integer, intent(out)                        :: info
      character(len=25)                           :: tag
      integer                                     :: ioerr

    values = 0.d0

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_dp > ',trim(charfile),' could not be openend. no value is '
      print *,'           returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,values
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)

    if (info .ne. 0) then 
      values = default 
    endif

    return
  end subroutine f_read_store_dbl

!-------------------------------------------------------------------------------
! read characters:

  subroutine f_read_store_chr(chartag,charfile,values,default,info)
    implicit none
      character(len=*), intent(in)  :: chartag
      character(len=*), intent(in)  :: charfile
      character(len=*), intent(out) :: values
      character(len=*), intent(in)  :: default
      integer, intent(out)          :: info
      character(len=25)             :: tag
      integer                       :: ioerr

    values = ''

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_char > ',trim(charfile),' could not be openend. no value'
      print *,'             is returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,values
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)
    
    if (info .eq. 0) then
      values = trim(values)
    else 
      values = trim(default)
    endif

    return
  end subroutine f_read_store_chr

!-------------------------------------------------------------------------------
! read vectors:

  subroutine f_read_store_vec(chartag,charfile,values,default,info)
    implicit none
    character(len=*), intent(in)                :: chartag
    character(len=*), intent(in)                :: charfile
    double precision, dimension(:), intent(out) :: values
    double precision, dimension(:), intent(in)  :: default
    integer, intent(out)                        :: info
    character(len=25)                           :: tag
    integer                                     :: ioerr

    values = 0.d0

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_vec > ',trim(charfile),' could not be openend. no value is '
      print *,'             returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,values
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)

    if (info .ne. 0) then
      values = default 
    endif

    return
  end subroutine f_read_store_vec

!-------------------------------------------------------------------------------
! read logicals:

  subroutine f_read_store_log(chartag,charfile,values,default,info)
    implicit none
      character(len=*), intent(in)    :: chartag
      character(len=*), intent(in)    :: charfile
      logical, intent(out)            :: values
      logical, intent(in)             :: default
      integer, intent(out)            :: info
      character(len=30)               :: tag
      integer                         :: ioerr

    values = .false.

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_char > ',trim(charfile),' could not be openend. no value'
      print *,'             is returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        info = 0
        values = .true.
        exit
      endif
    enddo
    close(12)
    if (info .ne. 0) then
      values = default 
    endif

    return
  end subroutine f_read_store_log

!-------------------------------------------------------------------------------
! trim integer string:
  
  function f_trimmed_int(int)
    implicit none
    character(len=15) :: f_trimmed_int
    integer           :: int
    character(len=30) :: dchr

    write(dchr,'(I15)') int
    call f_stripspaces(dchr)
    f_trimmed_int = trim(dchr)

    return
  end function f_trimmed_int

!-------------------------------------------------------------------------------
! trim character string:

  function f_trimmed_chr(chr)
    implicit none
    character(len=15) :: f_trimmed_chr
    character(len=30) :: chr
    character(len=30) :: dchr

    dchr = trim(chr)
    call f_stripspaces(dchr)
    f_trimmed_chr = trim(dchr)

    return
  end function f_trimmed_chr

!-------------------------------------------------------------------------------
! trim logical string:

  function f_trimmed_log(log)
    implicit none
    character(len=15) :: f_trimmed_log
    logical           :: log

    if (log) then
      f_trimmed_log = 'true'
    else
      f_trimmed_log = 'false'
    endif

    return
  end function f_trimmed_log

!-------------------------------------------------------------------------------
! trim double string:

  function f_trimmed_dbl(dbl)
    implicit none
    character(len=15) :: f_trimmed_dbl
    double precision  :: dbl
    character(len=30) :: dchr

    write(dchr,'(F15.5)') dbl
    call f_stripspaces(dchr)
    f_trimmed_dbl = trim(dchr)

    return
  end function f_trimmed_dbl
 
!-------------------------------------------------------------------------------
! trim vector:

  function f_trimmed_vec(dbl)
    implicit none
    character(len=15)               :: f_trimmed_vec
    double precision, dimension(3)  :: dbl
    character(len=15)               :: dchr

    write(dchr,100) dbl
100 format(F4.1,' ',F4.1,' ',F4.1)
    f_trimmed_vec = dchr 

    return
  end function f_trimmed_vec

!-------------------------------------------------------------------------------
! errors

  subroutine f_critical_(chr,line)
    implicit none 
    character(len=*), intent(in) :: chr 
    integer, intent(in)           :: line 

    print *,'    > critical error in line ',f_trimmed_(line)
    print *,'      of file ',trim(chr)
    stop 

    return 
  end subroutine f_critical_

!-------------------------------------------------------------------------------
! create a file name based on the initial file name provided by the user.

  subroutine f_getfilename(chrin,chrout)
    implicit none 
    character(len=*) :: chrin, chrout 
    character(len=4) :: chrnum
    integer :: i 
    logical :: res

    inquire(file=trim(chrin), exist=res)
    if (.not. res) then 
      chrout = chrin 
      return
    endif

    do i=1,5000
      write(chrnum,'(I4)') i 
      chrout = trim(chrnum) // trim(chrin)
      call f_stripspaces(chrout)

      inquire(file=trim(chrout), exist=res)
      if (.not. res) then 
        exit 
      endif 
    enddo
  end subroutine f_getfilename

!...............................................................................

end module functions