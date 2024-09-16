! explore.f90 - elapid - david gallina - dg@physik.uni-kassel.de - 2022
!-------------------------------------------------------------------------------
! algorithm to create a connected network of local minima and adjacent first-
! order saddle points for a given potential energy surface. more information on
! the routine can be found in the faq.
!-------------------------------------------------------------------------------
! double parameters:
! p_perturb - perturbation in random ef searches.
! p_offset - offset at the saddle point along the single negative eigenvalue.
! p_deps - accuracy of ef and lbfgs searches - rmsgrad < p_deps.
! p_dvec - criterion to check if two states are different - based on configs.
! p_dnrg - criterion to check if two states are different - based on energy.
!
! integer parameters:
! p_nmodes - number of eigenmodes to follow in ef searches.
! p_nrands - number of random directions to follow in ef searches.
! p_npool - pool of local minima to start the algorithm.
! p_iter - number of iterations.
! p_nm - maximum number of storable local minima.
! p_ns - maximum number of storable saddle points.
!
! boolean parameters:
! p_timerev_add - add states found through time-reversal symmetry.
! p_print_rates - print transition rates obtained through TST in 'rates.data'.
! p_print_order - print order parameters of all local minima in 'order.data'
! p_print_cfgs - print configurations of all states in '*.pts' files.
!-------------------------------------------------------------------------------

module explore
  use constants 
  use functions
  use potential
  use lbfgs
  use efol
  use dneb
  use omp_lib
  implicit none

! explore parameters.

  double precision, private :: p_deps, p_dvec, p_dnrg, p_offset, p_delta
  double precision, private :: p_scale
  integer, private, allocatable :: itoj(:)
  integer, private :: p_modes, p_bands, p_pool, p_nm, p_ns, p_out
  logical, private :: p_write_rates, p_write_cfgs, p_write_scale, p_timerev

! storage arrays.

  type(t_lm), allocatable :: mvec(:)
  type(t_ts), allocatable :: svec(:)

! control variables.

  integer :: nm, ns, storednm, storedns

contains

!-------------------------------------------------------------------------------

  subroutine explore_()
    implicit none
    double precision :: xopt(C_N), energy, g(C_N)
    type(t_lm) :: lmdummy 
    type(t_ts) :: tsdummy
    integer :: pos, i, j, iter
    logical :: conv, sorted
    
! read parameters:

    call explore_readparms()

! output.

    print *,''
    print *,'    >>> explore settings...'
    print *,'      > max. number of lm:      ',trim(f_trimmed_(p_nm))
    print *,'      > max. number of ts:      ',trim(f_trimmed_(p_ns))
    print *,'      > number of pool lm:      ',trim(f_trimmed_(p_pool))
    print *,'      > number of ef searches:  ',trim(f_trimmed_(p_modes))
    print *,'      > number of neb searches: ',trim(f_trimmed_(p_bands))
    print *,''

! allocate minima and saddle point arrays.

    allocate(mvec(p_nm))
    do i=1,p_nm
      allocate(mvec(i)%x(C_N))
    enddo
    nm = 0

    allocate(svec(p_ns))
    do i=1,p_ns
      allocate(svec(i)%x(C_N))
    enddo
    ns = 0

    if (CONTINUE_LOG) then 

! get minima and saddle points from file and set iteration to iter.

      call explore_initfromfile(iter)
      print *,'' 
      print *,'    >>> continuing exploration'
      print *,'      > current interation:  ',f_trimmed_(iter)
      print *,'      > current number of lm:',f_trimmed_(storednm)
      print *,'      > current number of ts:',f_trimmed_(storedns)
      print *,''
    else

! create initial pool of local minima.

      do i=1,p_pool
        xopt = potential_random()
        call lbfgs_(C_N,5,xopt,1.d-6,conv,energy,1000,j,.false.)
        if (.not. conv) then
          cycle
        endif
        call explore_add_lm(xopt,pos)
        call potential_(xopt,energy,g)

        if (p_timerev) then
          xopt = potential_reversal(xopt)
          call lbfgs_(C_N,5,xopt,1.d-6,conv,energy,1000,j,.false.)
          if (.not. conv) then
            cycle
          endif
          call explore_add_lm(xopt,pos)      
        endif 
      enddo

! create output files.

      open(unit=11,file='min.data',action='write',status='replace')
      close(11)
      open(unit=11,file='min.pts',action='write',status='replace')
      close(11)
      open(unit=11,file='ts.data',action='write',status='replace')
      close(11)
      open(unit=11,file='ts.pts',action='write',status='replace')
      close(11)

! set iteration to 1.

      iter = 1
      storednm = 0
      storedns = 0
    endif

! start of the main algorithm.

    do while (iter .le. nm)

! every p_out iterations, the located stationary configurations are stored.

      if (mod(iter-1,p_out) .eq. 0) then
        call explore_writetofile(iter)

        print *,''
        print *,'    explore > explored lm:   ',trim(f_trimmed_(iter))
        print *,'            > total lm:      ',trim(f_trimmed_(nm))
      endif

      call explore_search_ef(iter)

      call explore_search_neb(iter)
    
      iter = iter + 1
    enddo

! sort minima and saddle points.

    sorted = .false.
    do while (.not. sorted)
      sorted = .true.
      do i=1,nm-1
        if (mvec(i)%e .gt. mvec(i+1)%e) then 
          lmdummy = mvec(i+1)
          mvec(i+1) = mvec(i)
          mvec(i) = lmdummy
          sorted = .false.
        endif
      enddo 
    enddo

    sorted = .false.
    do while (.not. sorted)
      sorted = .true.
      do i=1,ns-1
        if (svec(i)%e .gt. svec(i+1)%e) then 
          tsdummy = svec(i+1)
          svec(i+1) = svec(i)
          svec(i) = tsdummy
          sorted = .false.
        endif 
      enddo 
    enddo 

    allocate(itoj(nm))
    do i=1,nm 
      itoj(mvec(i)%pos) = i 
    enddo 

    do i=1,ns 
      svec(i)%iloc(1) = itoj(svec(i)%iloc(1))
      svec(i)%iloc(2) = itoj(svec(i)%iloc(2))
    enddo

! store minima and saddle points.

    open(unit=15,file='min.data',action='write',status='replace')
    do i=1,nm
      write(15,'(F15.5,F15.5,I8,F12.5,F12.5,F12.5)') &
            mvec(i)%e/p_scale,mvec(i)%s,i,mvec(i)%m
    enddo
    close(15)

    open(unit=15,file='ts.data',action='write',status='replace')
    do i=1,ns
      write(15,'(F15.5,F15.5,I3,I8,I8,I3,I3,I3)') &
            svec(i)%e/p_scale,svec(i)%s,0,svec(i)%iloc,0,0,0
    enddo
    close(15)

! store the associated configurations.

    if (p_write_cfgs) then 
      open(unit=15,file='min.pts',action='write',status='replace')
      do i=1,nm 
        write(15,*) mvec(i)%x 
      enddo 
      close(15)

      open(unit=15,file='ts.pts',action='write',status='replace')
      do i=1,ns 
        write(15,*) svec(i)%x 
      enddo 
      close(15)
    endif 

! store scaled energies.

    if (p_write_scale) then 
      open(unit=15,file='sclmin.data',action='write',status='replace')
      do i=1,nm
        write(15,'(F15.5,F15.5,I8,F12.5,F12.5,F12.5)') &
              mvec(i)%e/C_SCALE,mvec(i)%s,i,mvec(i)%m
      enddo
      close(15)
  
      open(unit=15,file='sclts.data',action='write',status='replace')
      do i=1,ns
        write(15,'(F15.5,F15.5,I3,I8,I8,I3,I3,I3)') &
              svec(i)%e/C_SCALE,svec(i)%s,0,svec(i)%iloc,0,0,0
      enddo
      close(15)
    endif

! deallocate:

    deallocate(mvec)
    deallocate(svec)

    return
  end subroutine explore_

!...............................................................................
! routine to add a local minimum to the current database of local minima.

  subroutine explore_add_lm(x,pos)
    implicit none
    double precision, dimension(C_N), intent(in) :: x
    integer, intent(out) :: pos
    double precision :: g(C_N), s, e, m(3)
    integer :: i, nev

! at first, we calculate the energy and vibrational entropy.

    call potential_(x,e,g)
    call potential_entropy(x,s,nev)
    m = potential_moment(x) 
    if (nev .gt. 0) then
      pos = -1
      return
    endif

! see if lm is new.

    pos = 0
    if (nm .gt. 0) then
      do i=1,nm 
        if (ABS(mvec(i)%e - e) .lt. p_dnrg) then
          if (potential_distance(mvec(i)%x,x) .lt. p_dvec) then   
            pos = i
            exit 
          endif
        endif
      enddo
    endif
    
    if (pos .eq. 0) then
      if (nm + 1 .le. p_nm) then
        nm = nm + 1
        mvec(nm)%x = x 
        mvec(nm)%e = e
        mvec(nm)%s = s  
        mvec(nm)%m = m
        mvec(nm)%pos = nm
      else 
        pos =-1
      endif
    endif

    return
  end subroutine explore_add_lm

!...............................................................................
! routine to add a first-order saddle point to the current database:

  subroutine explore_add_ts(x,evec)
    implicit none
    double precision, intent(in) :: x(C_N), evec(C_N)
    double precision :: g(C_N), xp(C_N), xm(C_N), m(3), s, e, ep, em
    integer :: mpos, ppos, nev, i, j
    logical :: conv

! at first, we calculate the energy and moment.

    call potential_(x,e,g)
    if (ns .gt. 0) then
      do i=1,ns 
        if (abs(svec(i)%e - e) .lt. p_dnrg) then
          if (potential_distance(svec(i)%x,x) .lt. p_dvec) then
            return 
          endif 
        endif 
      enddo
    endif

! check if array is full.

    if (ns+1 .gt. p_ns) then
      return
    endif

! if the state is new, we calculate its vibrational entropy.

    call potential_entropy(x,s,nev)

! we calculate the adjacent minima and its time-reversed configuration.

    xp = x + p_offset * evec
    call lbfgs_(C_N,5,xp,p_deps,conv,ep,5000,j,.false.)
    if (.not. conv) then
      return
    endif
    call explore_add_lm(xp,ppos)

    xm = x - p_offset * evec
    call lbfgs_(C_N,5,xm,p_deps,conv,em,5000,j,.false.)
    if (.not. conv) then
      return
    endif
    call explore_add_lm(xm,mpos)

    if ((mpos .le. 0) .or. (ppos .le. 0)) then
      return
    endif

    if (mpos .eq. ppos) then 
      return 
    endif

! we update the transition state array:

    ns = ns + 1
    svec(ns)%x = x
    svec(ns)%e = e
    svec(ns)%s = s
    svec(ns)%iloc = (/ mpos, ppos /)

    if (p_timerev) then 
      xp = potential_reversal(xp)
      call lbfgs_(C_N,5,xp,p_deps,conv,ep,5000,j,.false.)
      if (.not. conv) then
        return
      endif
      call explore_add_lm(xp,ppos)
  
      xm = potential_reversal(xm)
      call lbfgs_(C_N,5,xm,p_deps,conv,em,5000,j,.false.)
      if (.not. conv) then
        return
      endif
      call explore_add_lm(xm,mpos)
    endif

    return
  end subroutine explore_add_ts

!...............................................................................
! routine to find first-order saddle points from a given local minimum ipos by 
! employing a number of eigenvector following searches.

  subroutine explore_search_ef(ipos)
    implicit none
    integer, intent(in) :: ipos
    double precision :: xopt(C_N), evec(C_N), prod 
    integer :: i, mode 
    logical :: conv 

! searching along the modes of the hessian of the local minimum ipos.

!$omp parallel do private(xopt,evec,prod,conv)

    do mode=-p_modes,p_modes 
      if (mode .eq. 0) then 
        cycle 
      endif 

      xopt = mvec(ipos)%x
      call efol_(xopt,C_N,mode,evec,prod,conv)
      if (.not. conv) then
        cycle 
      endif
      !print *,xopt(1),xopt(2)
!$omp critical (addTS) 

!$omp flush(ns)

      call explore_add_ts(xopt,evec)

!$omp end critical (addTS)

      if (p_timerev) then 
        xopt = potential_reversal(xopt)
        call efol_(xopt,C_N,1,evec,prod,conv)
        if (.not. conv) then 
          cycle 
        endif 

!$omp critical (addTS)
  
!$omp flush(ns)

        call explore_add_ts(xopt,evec)

!$omp end critical (addTS)
      
      endif
    enddo

!$omp end parallel do 

    return
  end subroutine explore_search_ef

!...............................................................................
! routine to find first-order saddle points from a given local minimum by 
! performing a number of neb searches.

  subroutine explore_search_neb(ipos)
    implicit none 
    integer, intent(in) :: ipos 
    double precision :: distances(p_bands), dist, delta, xts(C_N), evec(C_N)
    double precision :: prod
    integer :: ineighs(p_bands), i, j
    logical :: conv

! initialize neighbors and distances.

    ineighs = 0 
    distances = 1.d8
    
! find closest neighbors whose energy differs by p_delta or less.
    
    do i=1,nm 
      if (i .eq. ipos) then 
        cycle 
      endif
      
      dist = potential_distance(mvec(ipos)%x,mvec(i)%x)
      delta = ABS(mvec(ipos)%e - mvec(i)%e)
      do j=1,p_bands
        if ((dist .lt. distances(j)) .and. (delta .lt. p_delta)) then 
          distances(j+1:p_bands) = distances(j:p_bands-1)
          distances(j) = dist 
          ineighs(j+1:p_bands) = ineighs(j:p_bands-1)
          ineighs(j) = i 
        endif 
      enddo
    enddo

!$omp parallel do private(xts,evec,prod,conv)

    do i=1,p_bands
      if (ineighs(i) .eq. 0) then 
        cycle 
      endif 

      call dneb_(mvec(ipos)%x,mvec(ineighs(i))%x, xts, 11, spring=1.d0, & 
                 tangent='imptau',forces='dneb',optimizer='verlet', & 
                 dnebfactor=1.d-1)
      call efol_(xts,C_N,1,evec,prod,conv)
      if (.not. conv) then 
        cycle 
      endif 

!$omp critical (addTS)

!$omp flush(ns)

      call explore_add_ts(xts,evec)
    
!$omp end critical (addTS)

      if (p_timerev) then 
        xts = potential_reversal(xts)
        call efol_(xts,C_N,1,evec,prod,conv)
        if (.not. conv) then 
          cycle 
        endif 
    
!$omp critical (addTS)

!$omp flush(ns)

        call explore_add_ts(xts,evec)

!$omp end critical (addTS)

      endif
    enddo

!$omp end parallel do 

    return 
  end subroutine explore_search_neb

!...............................................................................
! init from file. 

  subroutine explore_initfromfile(iteration)
    implicit none 
    integer :: iteration, i, j, ierr

! local minima.

    open(unit=11,file='min.data',status='old',action='read',iostat=ierr)
    if (ierr .ne. 0) then 
      print *,'    > files are missing for continuation job. aborted.'
      stop 
    endif 

    open(unit=12,file='min.pts',status='old',action='read',iostat=ierr)
    if (ierr .ne. 0) then 
      print *,'    > files are missing for continuation job. aborted.'
      stop 
    endif 

    do 
      read(11,*,iostat=ierr) 
      if (ierr .ne. 0) then 
        exit 
      else 
        nm = nm + 1
      endif 
    enddo 

    rewind(11)
    do i=1,nm 
      print *,i,nm
      read(11,*) mvec(i)%e,mvec(i)%s,j,mvec(i)%m 
      read(12,*) mvec(i)%x 
    enddo 

    close(11)
    close(12)

! transition states.

    open(unit=11,file='ts.data',status='old',action='read',iostat=ierr)
    if (ierr .ne. 0) then 
      print *,'    > files are missing for continuation job. aborted.'
      stop 
    endif 

    open(unit=12,file='ts.pts',status='old',action='read',iostat=ierr)
    if (ierr .ne. 0) then 
      print *,'    > files are missing for continuation job. aborted.'
      stop 
    endif 

    do 
      read(11,*,iostat=ierr) 
      if (ierr .ne. 0) then 
        exit 
      else 
        ns = ns + 1
      endif 
    enddo 

    rewind(11)
    do i=1,ns 
      read(11,*) svec(i)%e,svec(i)%s,j,svec(i)%iloc 
      read(12,*) svec(i)%x 
    enddo 

    close(11)
    close(12)

    open(unit=11,file='STATUS',action='read',status='old',iostat=ierr)
    if (ierr .ne. 0) then 
      print *,'    > files are missing for continuation job. aborted.'
      stop 
    endif 

    read(11,*) iteration, storednm, storedns
    close(11)

    nm = storednm 
    ns = storedns

    return 
  end subroutine explore_initfromfile

!...............................................................................

  subroutine explore_writetofile(iteration)
    implicit none 
    integer :: i, iteration

    open(unit=11,file='min.data',status='old',action='write',position='append')
    open(unit=12,file='min.pts',status='old',action='write',position='append')
    do i=storednm+1,nm 
      write(11,*) mvec(i)%e,mvec(i)%s,i,mvec(i)%m 
      write(12,*) mvec(i)%x
    enddo 
    close(11)
    close(12)

    open(unit=11,file='ts.data',status='old',action='write',position='append')
    open(unit=12,file='ts.pts',status='old',action='write',position='append')
    do i=storedns+1,ns 
      write(11,*) svec(i)%e,svec(i)%s,0,svec(i)%iloc
      write(12,*) svec(i)%x
    enddo 
    close(11)
    close(12)

    open(unit=11,file='STATUS',action='write',status='replace')
    write(11,*) iteration,nm,ns 
    close(11)
    storednm = nm 
    storedns = ns
    return 
  end subroutine explore_writetofile

!...............................................................................
! parameters of the explore module.

  subroutine explore_readparms()
    implicit none
    integer           :: info

! convergence criteria.

    p_deps = 1.d-6 
    call f_read_('explore_deps','params.dat',p_deps,info)

    p_dvec = 1.d-2 
    call f_read_('explore_dvec','params.dat',p_dvec,info)

    p_dnrg = 1.d-4
    call f_read_('explore_dnrg','params.dat',p_dnrg,info)

! explore criteria.

    p_pool = 100
    call f_read_('explore_pool','params.dat',p_pool,info)

    p_modes = 5 
    call f_read_('explore_modes','params.dat',p_modes,info)

    p_bands = 5
    call f_read_('explore_bands','params.dat',p_bands,info)

    p_nm = 10000
    call f_read_('explore_nm','params.dat',p_nm,info)

    p_ns = 50000
    call f_read_('explore_ns','params.dat',p_ns,info)

    p_out = 100
    call f_read_('explore_out','params.dat',p_out,info)

! neb criteria.

    p_delta = 80.d0
    call f_read_('explore_delta','params.dat',p_delta,info)

! connection criteria.

    p_offset = 1.d-2
    call f_read_('explore_offset','params.dat',p_offset,info)

! storage criteria.

    p_write_cfgs = .true.
    call f_read_('explore_write_cfgs','params.dat',p_write_cfgs,info)

    p_write_rates = .false.
    call f_read_('explore_write_rates','params.dat',p_write_rates,info)

    p_write_scale = .false.
    call f_read_('explore_write_scale','params.dat',p_write_scale,info)

    p_timerev = .true.
    call f_read_('explore_timerev','params.dat',p_timerev,info)

    p_scale = 1.d0
    call f_read_('explore_scale','params.dat',p_scale,info)

    return
  end subroutine explore_readparms

!...............................................................................

end module explore





