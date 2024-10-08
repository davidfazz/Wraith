! efol.f90.
! WRAITH - david gallina - gallina@uni-kassel.de - last update: 07/12/23.
!...............................................................................
! implementation of the eigenvector-following algorithm [wal98].
!...............................................................................

module efol
  use constants 
  use functions
  use potential
  implicit none

!...............................................................................
! parameters that are necessary for the algorithm.

  double precision, parameter, private    :: efol_pushoff = 2.d-1
  double precision, parameter, private    :: efol_maxmax = 5.d-1
  double precision, parameter, private    :: efol_minmax = 1.d-2
  double precision, parameter, private    :: efol_trad = 2.d-1
  double precision, parameter, private    :: efol_pushcut = 1.d-3
  double precision, parameter, private    :: efol_zero = 1.d-6
  double precision, parameter, private    :: efol_overlap = 0.8d0
  double precision, parameter, private    :: efol_maxstep = 0.2d0
  integer, parameter, private             :: efol_maxit = 1000
  logical, parameter, private             :: efol_debug = .false.

contains

!...............................................................................
! eigenvector-following routine:

  subroutine efol_(xopt,n,mode,evec,prod,conv)
    implicit none
    double precision, intent(inout) :: xopt(n)
    double precision, intent(out) :: evec(n), prod 
    integer, intent(in) :: n, mode 
    logical, intent(out) :: conv 
    double precision :: xhess(n,n), evecs(n,n), xmin(n), gnew(n), step(n) 
    double precision :: cstep(n), pstep(n), stpmax(n), evals(n), fob(n), pfob(n)
    double precision :: rat(n), vec(n), tempa(n), stpmag, sstpmag, enew, olap 
    double precision :: maxolap, drms, rat1, rat2, r1, r2, sum, lp, lp1, lp2 
    double precision :: scale 
    integer :: iter, i, j, nnev, imode, omode, im, icount 
    logical :: zt(n), pzt(n), away 

! initialising variables.
!   > stpmax is the maximum step length for each direction.
!   > vec is the current eigendirection.

    stpmax = efol_maxstep
    vec = 0.d0
    if (ABS(mode) .gt. n) then
      print *,'    efol > defined mode is larger than number of dimensions.'
      print *,'           mode =',mode
      conv = .false.
      return
    endif

! some additional variables.
!   > xmin contains the initial minimum.

    xmin = xopt

! beginning of the optimisation cycle.
!   > iter is the number of the current iteration.

    iter = 1
    do

! we calculate the current energy and gradient as well as the root mean square 
! of the gradient to see if we are converged to a critical point.
!   > xopt contains the current configuration.
!   > gnew contains the current gradient.
!   > enew contains the current energy.
!   > drms contains the current rms gradient.

      call potential_(xopt,enew,gnew)
      drms = f_rms(gnew)

! we calculate the hessian matrix and diagonalise it.
!   > evecs contains the eigenvectors.
!   > evals contains the eigenvalues.

      call hessian_(xopt,xhess)
      call f_dsyev(xhess,evecs,evals)
      evec = evecs(1,:)

! we remove zero eigenvalues to avoid some problems in calculating step sizes. 
! in addition, we calculate the number of negative eigenvalues
!   > zt is true, if the eigenvalue is not zero, and false if it is zero
!   > nnev is the number of negative eigenvalues

      zt = .true.
      do i=1,n
        if (ABS(evals(i)) .lt. efol_zero) then
          zt(i) = .false.
        endif
      enddo

      nnev = 0
      do i=1,n
        if ((zt(i)) .and. (evals(i) .lt. 0.d0)) then
          nnev = nnev + 1
        endif
      enddo

! we find the mode we want to follow uphill. if we are in the first iteration, 
! we take the input mode as mode. in the following iterations, we take the mode 
! with the largest overlap as uphill direction. if the overlap is very small,
! there are different cases. if the new mode is equal or softer than the old 
! mode, we take the new mode. otherwise, we take the old mode.
!   > imode is the new mode
!   > omode is the old mode
!   > olap is the overlap between the eigenvectors and the old mode
!   > maxolap is the maximum overlap
!   > im is the mode with the maximum overlap

      imode = 0
      if (mode .ne. 0) then
        if (zt(ABS(mode))) then
          imode = ABS(mode)
        else
          do i=1,n
            if (zt(i)) then
              imode = i
              exit
            endif
          enddo
        endif
      elseif (mode .eq. 0) then
        imode = 1
      endif

      if ((iter .gt. 1) .and. (mode .ne. 0)) then
        maxolap = 0.d0
        do i=1,n
          olap = DOT_PRODUCT(vec,evecs(i,:))
          if ((ABS(olap) .gt. maxolap) .and. (zt(i))) then
            maxolap = ABS(olap)
            im = i
          endif
        enddo

        if (maxolap .lt. efol_overlap) then
          if (im .lt. omode) then
            if (zt(im)) then
              im = im
            endif
          else 
            if (zt(omode)) then 
              im = omode 
            endif
          endif
        endif

        if (im .gt. omode+8) then
          im = omode
        endif

        if (zt(im)) then
          imode = im
        elseif (zt(omode)) then
          imode = omode
        endif
        omode = imode
      else
        omode = imode
      endif

      if ((iter .gt. 10) .and. (nnev .ge. 1)) then
        imode = 1
      endif

! we calculate the projection of the gradient on the different eigenvectors of 
! the hessian
!   > fob contains the projection for each eigenvector

      fob = 0.d0
      do i=1,n
        fob(i) = dot_product(gnew,evecs(i,:))
      enddo

      if (efol_debug) then
        write(*,*) '    efol > eigenvector projection of the gradient:'
        write(*,'(6(F15.5))') fob
        write(*,*) '    efol > number of negative eigenvalues:',nnev
        write(*,*) '    efol > this mode will be searched uphill:',imode,evals(imode)
      endif

! we calculate the maximum step size for each direction based on a trust radius 
! scheme that was described in "theoretical study of the water pentamer" by 
! d.j. wales and t.r. walsh. at each step we compare the calculated eigenvalue 
! with an estimated eigenvalue from the gradients. if the difference is large, 
! we decrease the step size, if the difference is small, we increase the step 
! size.
!   > rat contains the calculated radius for each direction
!   > tempa is a temporary array
!   > pzt is the zt array from the previous run
!   > pfob is the fob array from the previous run
!   > r1/r2 are the ratio between the estimated and the calculated eigenvalue
!   > rat1/rat2 are the radii for r1/r2

      if (iter .eq. 1) then
        rat = 0.d0
      endif

      j = 0
      if (iter .gt. 1) then
        tempa = stpmax
        do i=1,n
          rat(i) = 0.d0
          if (zt(i)) then
            do while (j .lt. n)
              j = j + 1
              if ((.not. pzt(j)) .and. (j .lt. n)) then
                cycle
              endif

              if (abs(pstep(j)) .gt. 1.d-40) then
                r1 = (fob(i) - pfob(j)) / (evals(i)*pstep(j))
                r2 = (-fob(i) - pfob(j)) / (evals(i)*pstep(j))

                rat1 = ABS(ABS(r1) - 1.d0)
                rat2 = ABS(ABS(r2) - 1.d0)
                rat(i) = MIN(rat1,rat2)

                if (efol_debug) then
                  write(*,'(I4,7(F15.8))') j,fob(i),pfob(i),evals(i),r1,r2,rat(i),pstep(j)
                endif

                if (rat(i) .gt. efol_trad) then
                  stpmax(i) = MAX(tempa(j)/1.11d0,efol_minmax)
                else
                  stpmax(i) = MIN(MAX(tempa(j)*1.09d0,efol_minmax),efol_maxmax)
                endif
              endif
              exit
            enddo
          endif
        enddo
      endif

! we calculate the mean of the eigenvalues and the logarithm of the product of 
! the eigenvalues
!   > sum is the sum of eigenvalues
!   > prod is the logarithm of the product of the eigenvalues
!   > vec is the current eigenmode

      sum = 0.d0
      prod = 0.d0
      do i=1,n
        if (zt(i)) then
          if (evals(i) .gt. 0.d0) then
            sum = sum + ABS(evals(i))
            prod = prod + LOG(evals(i))
          endif
        endif
      enddo
      sum = sum/DBLE(n-nnev)
      vec = evecs(imode,:)

      if (efol_debug) then
        write(*,*) '    efol > mean of positive eigenvalues:',sum
        write(*,*) '    efol > log product of positive eigenvalues:',prod
      endif

! if we converge to a critical point, we step off if it has the wrong number of 
! negative eigenvalues

      away = .false.
      if (drms .lt. efol_pushcut) then
        if (nnev .ne. 1) then
          if (mod(iter-1,4) .eq. 0) then
            away = .true.
          endif
        endif
      endif

! we calculate the lagrange multiplier for each direction as is described in the 
! paper "rearrangements of 55-atom lennard-jones and (C60)55 clusters" by 
! d.j. wales.
!   > step is the calculated step length for each direction

      icount = 0
      do i=1,n
        step(i) = 0.d0
        if (zt(i)) then
          lp1 = ABS(evals(i))/2.d0
          lp2 = 1.d0 + 4.d0*(fob(i)/evals(i))**2
          lp = lp1*(1.d0 + SQRT(lp2))

          if (i .eq. imode) then
            lp =-lp
            icount = icount + 1
          elseif (icount .lt. 1) then
            lp =-lp
            icount = icount + 1
          endif

          step(i) = -fob(i)/lp
        endif
      enddo

      if (efol_debug) then
        do i=1,n
          write(*,*) '    efol > unscaled step for mode',i,step(i)
        enddo
      endif

! if the algorithm converged to a critical point with wrong index, we step off.

      if (away) then
        if (nnev .eq. 0) then
          if (efol_pushoff .ne. 0.d0) then
            step(imode) = efol_pushoff
          else
            step(imode) = stpmax(imode)/10.d0
          endif

          if (mode .lt. 0) then
            step(imode) = -step(imode)
          endif

        elseif (nnev .gt. 1) then
          if (mode .eq. 0) then
            do i=1,n
              if ((zt(i)) .and. (evals(i) .lt. 0.d0)) then
                if (efol_pushoff .ne. 0.d0) then
                  step(i) = efol_pushoff
                else
                  step(i) = stpmax(i) / 10.d0
                endif
              endif
            enddo
          else
            if (efol_pushoff .ne. 0.d0) then
              step(imode) = efol_pushoff
            else
              step(imode) = stpmax(imode) / 10.d0
            endif
            if (mode .lt. 0) then
              step(imode) = -step(imode)
            endif
          endif
        endif
      endif

! we calculate the magnitude of the steps and the scaling factor
!   > stpmag is the maximum unscaled step length
!   > sstpmag is the maximum scaled step length

      stpmag = 1.d-10
      do i=1,n
        if (ABS(step(i)) .gt. stpmag) then
          stpmag = ABS(step(i))
        endif
      enddo

      do i=1,n
        scale = MIN(stpmax(i)/MAX(ABS(step(i)),1.d-10),1.d0)
        if ((i .ne. imode) .or. (.not. away)) then
          step(i) = scale*step(i)
        endif

        if (ABS(evals(i)) .lt. efol_zero) then
          step(i) = 0.d0
        endif
      enddo

      sstpmag = 1.d-10
      do i=1,n
        if (ABS(step(i)) .gt. sstpmag) then
          sstpmag = ABS(step(i))
        endif
      enddo
      scale = 1.d0

      if (efol_debug) then
        write(*,*) '    efol > Vector Gradient Secder Step Max step Trus ratio'
        do i=1,n
          write(*,300) i,fob(i),evals(i),step(i),stpmax(i),rat(i)
300       format(I4,1X,E15.6,1X,E15.6,1X,E13.6,1X,E13.6,1X,E15.6)
        enddo
      endif

      if (efol_debug) then
        write(*,*) '    efol > maximum scaled/unscaled step is:',sstpmag,stpmag
      endif

! we check if the algorith has converged. if not, we check if the maximum number 
! of iterations has been reached. if not, the configuration is updated by adding 
! cstep. the different arrays are stored and the loop is repeated.
!   > cstep is the current calculated eigenvector following step

      if ((nnev .eq. 1) .and. (drms .lt. 1.d-3)) then
        if (iter .gt. 1) then
          conv = .true.
          exit
        endif
      endif

      if (iter .gt. efol_maxit) then
        conv = .false.
        exit
      endif

      cstep = 0.d0
      do i=1,n
        cstep = cstep + step(i)*evecs(i,:)
      enddo

      xopt = xopt + cstep
      xopt = potential_rescale(xopt)

      iter = iter + 1

      pstep = step
      pfob = fob
      pzt = zt
    enddo

    return
  end subroutine efol_

!...............................................................................
    
end module efol
