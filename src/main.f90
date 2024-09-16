! main.f90 - elapid- david gallina - dg@physik.uni-kassel.de - 2022
!..............................................................................

program elapid
  use constants
  use functions
  use potential
  use explore
  use random
  use follow
  use hysteresis
  use rates
  use heatcap
  use mc
  use omp_lib
  implicit none 
  double precision    :: stime
  double precision    :: etime
  double precision    :: dtime

  print *,'                                                            '
  print *,'   _    _______  ___  _____ _____ _   _      __   _____     '
  print *,'  | |  | | ___ \/ _ \|_   _|_   _| | | |    /  | |  _  |    '
  print *,'  | |  | | |_/ / /_\ \ | |   | | | |_| |    `| | | | | |    '
  print *,'  | |/\| |    /|  _  | | |   | | |  _  |     | | | | | |    '
  print *,'  \  /\  / |\ \| | | |_| |_  | | | | | |    _| |_\ |_/ /    '
  print *,'   \/  \/\_| \_\_| |_/\___/  \_/ \_| |_/    \___(_)___/     '
  print *,'                                                            '
  print *,'                                                            '
                                                             
  stime = omp_get_wtime()

! initialize:

  call readparms_()

! intialise potential:  

  call potential_init()

! explore:

  !call checkpot_()
  call explore_()
  !call follow_()
 !jh call heatcap_()
 ! call mc_()
  !call rates_()
  !call hysteresis_()

  etime = omp_get_wtime()

  dtime = etime - stime

  call potential_deinit()

  print *,''
  print *,'    > wraith has finished.'

contains 

!..............................................................................
! read main parameters:

  subroutine readparms_()
    implicit none 
    character(len=30) :: string
    integer :: info, number
    logical :: true 
    
! potential parameters.

    dpls1d_log = .false.
    dpls2d_log = .false.
    open1d_log = .false. 
    chains_log = .false.
    unsurf_log = .false.

    call f_read_('potential','params.dat',string,info)
    if (trim(string) .eq. 'dpls1d') then 
      dpls1d_log = .true.
    elseif (trim(string) .eq. 'dpls2d') then 
      dpls2d_log = .true.
    elseif (trim(string) .eq. 'unsurf') then 
      unsurf_log = .true.
    elseif (trim(string) .eq. 'chains') then 
      chains_log = .true.
    elseif (trim(string) .eq. 'open1d') then 
      open1d_log = .true.
    endif 

! continue code.

    call f_read_('continue','params.dat',true,info)
    if (true) then 
      CONTINUE_LOG = .true.
    else 
      CONTINUE_LOG = .false.
    endif 

! rng.

    call f_read_('random_seed','params.dat',number,info)
    call random_set(number)

    return 
  end subroutine readparms_

!...............................................................................

end program elapid

