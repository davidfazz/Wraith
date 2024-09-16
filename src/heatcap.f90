module heatcap
  use constants
  use functions 
  implicit none 

contains 

  subroutine heatcap_
    implicit none 
    type(t_lm), allocatable :: mvec(:)
    double precision, allocatable :: pvec(:)
    double precision :: mine, maxs, zcan, kt, deltat, heat, mean, variance, square
    integer :: i, nm, it, ierr

! get minima.
    
    nm = 0
    open(unit=11,file='min3d.data',action='read',status='old')
    do
      read(11,*,iostat=ierr)
      if (ierr .eq. 0) then 
        nm = nm + 1
      else 
        exit 
      endif 
    enddo 

    rewind(11)
    allocate(mvec(nm))
    do i=1,nm 
      read(11,*) mvec(i)%e, mvec(i)%s 
    enddo 
    close(11)

! get minimum energy and maximum entropy.

    mine = 1.d8 
    maxs = 0.d0 
    do i=1,nm 
      if (mvec(i)%e .lt. mine) then 
        mine = mvec(i)%e 
      endif 
      if (mvec(i)%s .gt. maxs) then 
        maxs = mvec(i)%s 
      endif 
    enddo 

! calculate heat capacity for different values of T.

    open(unit=11,file='heatcapacity.dat',action='write',status='replace')

    allocate(pvec(nm))

    deltat = 5.d0
    do it=1,80 
      kt = C_KB * deltat * DBLE(it)

! calculate canonical partition function.

      zcan = 0.d0 
      do i=1,nm 
        zcan = zcan + SQRT(EXP(maxs - mvec(i)%s))*EXP(-(mvec(i)%e - mine)/kt)
      enddo 

! calculate probabilities.

      pvec = 0.d0
      do i=1,nm 
        pvec(i) = SQRT(EXP(maxs - mvec(i)%s))*EXP(-(mvec(i)%e - mine)/kt)/zcan
      enddo 

     ! f = -kt*144.d0*LOG(2.d0*C_PI*kt)

      print *,pvec(1:3)

! calculate variance.

      mean = 0.d0 
      do i=1,nm 
        mean = mean + pvec(i)*mvec(i)%e 
      enddo 

      square = 0.d0 
      do i=1,nm 
        square = square + pvec(i)*mvec(i)%e*mvec(i)%e
      enddo 

      variance = square - mean*mean 

      heat = variance / (kt * deltat * DBLE(it)) + 144.d0/2.d0*C_KB 
      mean = 144.d0/2.d0 * kt + mean

      write(11,'(F12.5,F12.5,F15.3,F15.3)') deltat*DBLE(it), heat, mean, variance
    enddo
    close(11)

    deallocate(mvec)
    deallocate(pvec)
    return
  end subroutine heatcap_

end module heatcap