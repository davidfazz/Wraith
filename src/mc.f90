module mc
  use constants 
  use functions 
  use potential 
  use random
  implicit none 
  double precision, parameter :: tlow = 10.d0, thigh = 210.d0
  integer, parameter :: ntmp = 20, nouter = 100000, maxsteps = 10000
  double precision, parameter :: mc_perturb = 5.d-2

contains 

  subroutine mc_() 
    implicit none 
    double precision :: xcfg(C_N), xnew(C_N), enew, eold, grad(C_N), phi, theta 
    double precision :: rval, moment(3), kt
    integer :: outer, i

    kt = C_KB * 100.d0 * 1000.d0
    xcfg = potential_random()
    call potential_(xcfg,eold,grad)

    open(unit=11,file='mc.out',action='write',status='replace')

! outer step.
    
    print *,eold,C_P
    do outer=1,1000000

! inner step - perturb xcfg.

      do i=1,C_P 
        phi = 2.d0*C_PI*(0.5d0 - random_uniform())*mc_perturb 
        theta = 1.d0*C_PI*(0.5d0 - random_uniform())*mc_perturb
        
        !theta = ACOS(1.d0 - 2.d0*random_uniform())*mc_perturb 
        xnew(2*i-1) = xcfg(2*i-1) + theta 
        xnew(2*i) = xcfg(2*i) + phi 
      enddo 
      xnew = potential_rescale(xnew)
      call potential_(xnew,enew,grad)

      if (enew .lt. eold) then 
        xcfg = xnew 
        eold = enew
      else 
        !print *,exp(-(eold-enew)/kt),eold-enew,kt
        rval = random_uniform() 
        if (rval .lt. exp(-(enew-eold)/kt)) then 
          xcfg = xnew 
          eold = enew 
        endif 
      endif 

      if (mod(outer,100) .eq. 0) then 
       ! stop
        if (outer .gt. 1000) then 
          moment = potential_moment(xcfg)
          write(11,'(F12.5,F12.5,F12.5,F12.5)') eold,moment
        endif 
      endif
    enddo

    close(11)

    return
  end subroutine mc_

end module mc 