module rates
  use potential
  use constants
  use functions 
  use omp_lib
  implicit none 
  double precision, parameter :: C_EDD = 49.0653d0
  type(t_lm), allocatable :: mvec(:) 
  type(t_ts), allocatable :: svec(:)
  integer :: nm, ns

contains 

  subroutine rates_()
    implicit none 
    double precision :: e, g(C_N)
    integer :: i, j, ierr

! get local minima.

    open(unit=11,file='min.data',action='read',status='old')
    open(unit=12,file='min3d.pts',action='read',status='old')
    !open(unit=12,file='min.pts',action='read',status='old')

    nm = 0
    do
      read(11,*,iostat=ierr)
      if (ierr .eq. 0) then 
        nm = nm + 1
      else 
        exit 
      endif 
    enddo 

    allocate(mvec(nm))
    do i=1,nm 
      allocate(mvec(i)%x(C_N))
    enddo 
    rewind(11)
    do i=1,nm 
      read(11,*) mvec(i)%e, mvec(i)%s 
      read(12,*) mvec(i)%x 
    enddo 
    close(11)
    close(12)

! update entropies.

    do i=1,nm 
      call rates_entropy(mvec(i)%x,mvec(i)%s)
    enddo

! get transition states:

    open(unit=11,file='ts.data',action='read',status='old')
    open(unit=12,file='ts3d.pts',action='read',status='old') 

    ns = 0
    do
      read(11,*,iostat=ierr)
      if (ierr .eq. 0) then 
        ns = ns + 1 
      else 
        exit 
      endif 
    enddo 

    allocate(svec(ns))
    do i=1,ns 
      allocate(svec(i)%x(C_N))
    enddo 
    rewind(11) 
    do i=1,ns 
      read(11,*) svec(i)%e, svec(i)%s, j, svec(i)%iloc(1), svec(i)%iloc(2)
      read(12,*) svec(i)%x 
    enddo 
    close(11)
    close(12)
    print *,ns

! calculate rates and saddle entropy for each transition.

!$omp parallel do

    do i=1,ns 
      call rates_entropy(svec(i)%x,svec(i)%s)
      call rates_htst(svec(i)%iloc(1),i,svec(i)%k(1))
      call rates_htst(svec(i)%iloc(2),i,svec(i)%k(2))
      if (mod(i,1000) .eq. 0) then 
        print *,i,ns 
      endif
    enddo 

!$omp end parallel do

    open(unit=11,file='ts.data',action='write',status='replace')
    do i=1,ns
      write(11,'(F15.7,F15.7,I6,I6,I6,F15.3,F15.3,I6)') svec(i)%e*144.d0, svec(i)%s, &
      0, svec(i)%iloc(1), svec(i)%iloc(2), svec(i)%k(1), svec(i)%k(2), 0
    enddo
    close(11)

    open(unit=11,file='min.data',action='write',status='replace')
    do i=1,nm 
      write(11,'(F15.5,F15.5,I6,I6,I6,I6)') mvec(i)%e*144.d0, mvec(i)%s, 0, 0, 0, 0
    enddo
    close(11)

    deallocate(mvec)
    deallocate(svec)


    return 
  end subroutine rates_

!...............................................................................

  subroutine rates_htst(ilm,its,rate)
    implicit none 
    double precision :: xlm(C_N), xts(C_N), hlm(C_N,C_N), hts(C_N,C_N)
    double precision :: eveclm(C_N,C_N), evects(C_N,C_N), evallm(C_N) 
    double precision :: evalts(C_N), alm(C_N), rate, tsum, psum, ci, di, fi, gi 
    double precision :: root, slm
    integer :: ilm, its, i, j 

! find configurations.

    xlm = mvec(ilm)%x 
    xts = svec(its)%x 

! find eigenvalues and eigenvectors.

    call hessian_(xlm,hlm)
    call hessian_(xts,hts)

    call f_dsyev(hlm,eveclm,evallm)
    call f_dsyev(hts,evects,evalts)

! calculate the a's from the bessarab paper.

    alm = 0.d0 
    do j=2,C_N 
      tsum = 0.d0 
      psum = 0.d0 
      do i=1,C_N,2 
        ci = C_GAMMA / (C_M*(1.d0 + C_ALPHA*C_ALPHA))
        di = ci 
        fi = ci*C_ALPHA 
        gi = fi 

        tsum = tsum + evects(i,1)*(-fi*evects(i,j) - di*evects(i+1,j))
        psum = psum + evects(i+1,1)*(di*evects(i,j)  - gi*evects(i+1,j))
      enddo
      alm(j) = evalts(j)*(tsum + psum)
    enddo

! calculate the root.

    root = 0.d0 
    do i=2,C_N 
      root = root + alm(i)*alm(i)/evalts(i)
    enddo 
    root = SQRT(root)

! calculate rate. 

    rate = root*SQRT(evallm(1))
    do i=2,C_N 
      rate = rate*SQRT(evallm(i)/evalts(i))
    enddo

    return 
  end subroutine rates_htst

!...............................................................................

  subroutine rates_entropy(x,s)
    implicit none 
    double precision :: x(C_N), h(C_N,C_N), eval(C_N), evec(C_N,C_N), s
    integer :: i 

    call hessian_(x,h)
    call f_dsyev(h,evec,eval)
    s = 0.d0
    do i=1,C_N 
      if (eval(i) .le. 1.d-10) then 
        cycle 
      endif 
      s = s + log(eval(i))
    enddo 

    return 
  end subroutine rates_entropy

!...............................................................................

end module rates