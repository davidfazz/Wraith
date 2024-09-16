module dneb 
  use constants
  use functions
  use potential
  implicit none

! different methods can be used in the dneb algorithm.
! > kind_tau chooses the method to calculate the tangent at each image. it has 
!   the options 
!   > "imptau": tangents are calculated according to the paper [hen00].
! > kind_frc chooses the method to calculate the forces at each image. it has 
!   the options
!   > "stdneb": forces are calculated according to the standard nudged elastic 
!               band method as found in [hen00].
!   >   "dneb": forces are calculated according to the doubly-nudged elastic 
!               band method as found in [try04].
! > kind_opt chooses the optimizer to optimize the band of images.
!   > "verlet": velocity verlet algorithm.   

  double precision :: kneb, dchi
  character(len=6) :: kind_tau, kind_frc, kind_opt 
  logical :: output = .true.

contains 

!...............................................................................
! nudged elastic band routine.

  subroutine dneb_(xst,xnd,xsp,images,spring,tangent,forces,optimizer,&
                   dnebfactor)
    implicit none 
    double precision, intent(in) :: xst(C_N), xnd(C_N)
    double precision, intent(out) :: xsp(C_N)
    double precision, optional :: spring, dnebfactor
    character(len=*), optional :: tangent, forces, optimizer
    integer :: images, ipos, i 
    double precision :: xpath(images*C_N), xgrad(images*C_N), epath(images)
    double precision :: e, g(C_N), m(3), xone(C_N), xtwo(C_N), distance

! standard parameters.

    if (.not. present(spring)) then 
      kneb = 1.d0 
    else 
      kneb = spring 
    endif

    if (.not. present(tangent)) then 
      kind_tau = 'imptau'
    else 
      kind_tau = trim(tangent)
    endif 

    if (.not. present(forces)) then 
      kind_frc = 'stdneb' 
    else 
      kind_frc = trim(forces)
    endif 

    if (.not. present(optimizer)) then 
      kind_opt = 'verlet'
    else 
      kind_opt = trim(optimizer)
    endif

    if (.not. present(dnebfactor)) then 
      dchi = 1.d-2
    else 
      dchi = dnebfactor
    endif


! create path of images.

    call dneb_path(xst,xnd,xpath,images)

! optimize path of images.

    call dneb_optimize(xpath,epath,images)

! find local maximum along path.

    ipos = MAXLOC(epath,1)
    xsp = xpath((ipos-1)*C_N+1:(ipos-1)*C_N+C_N)

    distance = 0.d0 
    do i=2,images 
      xone = xpath((i-1)*C_N+1:(i-1)*C_N+C_N)
      xtwo = xpath((i-2)*C_N+1:(i-2)*C_N+C_N)

      distance = distance + potential_distance(xone,xtwo)
    enddo

    return
  end subroutine 

!...............................................................................
! calculate the forces acting on each image of the chain of images.

  subroutine dneb_grad(xpath,xgrad,epath,images)
    implicit none 
    double precision :: xpath(images*C_N), xgrad(images*C_N), epath(images)
    double precision :: gtrue(images*C_N), xtau(images*C_N), gspring(images*C_N)
    integer :: images, i, i1, in

! calculate gradient and energies for each image.

    do i=1,images
      i1 = (i-1)*C_N + 1
      in = (i-1)*C_N + C_N 
      call potential_(xpath(i1:in),epath(i),gtrue(i1:in))   
    enddo

! calculate tangent vectors at each image.

    call dneb_tau(xpath,xtau,epath,images)

! calculate spring forces for each image.

    call dneb_spring(xpath,gspring,images)

! calculate total force acting on each image.

    xgrad = 0.d0 
    do i=2,images-1 
      i1 = (i-1)*C_N + 1
      in = (i-1)*C_N + C_N 

! perpendicular part of true gradient.

      xgrad(i1:in) = gtrue(i1:in) - DOT_PRODUCT(gtrue(i1:in),xtau(i1:in))*xtau(i1:in)

      if (trim(kind_frc) .eq. 'stdneb') then

! spring forces only along the path tangent.
        
        xgrad(i1:in) = xgrad(i1:in) + DOT_PRODUCT(gspring(i1:in),xtau(i1:in))*xtau(i1:in)
      elseif (trim(kind_frc) .eq. 'dneb') then         

! spring forces along the path tangent.

        xgrad(i1:in) = xgrad(i1:in) + DOT_PRODUCT(gspring(i1:in),xtau(i1:in))*xtau(i1:in)

! and having a small part perpendicular to the path tangent.

        xgrad(i1:in) = xgrad(i1:in) + dchi*(gspring(i1:in) - &
                       DOT_PRODUCT(gspring(i1:in),xtau(i1:in))*xtau(i1:in))
      endif
    enddo

    return 
  end subroutine dneb_grad

!...............................................................................
! calculate the spring gradient of each image.

  subroutine dneb_spring(xpath,gspring,images)
    implicit none 
    double precision :: xpath(images*C_N), gspring(images*C_N)
    integer :: images, i, i1, in, ip1, ipn, im1, imn 

    gspring = 0.d0
    do i=2,images-1
      ip1 = i*C_N + 1 
      ipn = i*C_N + C_N
      i1 = (i-1)*C_N + 1
      in = (i-1)*C_N + C_N
      im1 = (i-2)*C_N + 1
      imn = (i-2)*C_N + C_N 
      
      gspring(i1:in) = kneb*dneb_polarTau(xpath(i1:in),xpath(ip1:ipn)) &
                     - kneb*dneb_polarTau(xpath(im1:imn),xpath(i1:in))
    enddo 
    
    return 
  end subroutine dneb_spring

!...............................................................................

  subroutine dneb_optimize(xpath,epath,images)
    implicit none 
    double precision :: xpath(images*C_N), xgrad(images*C_N), epath(images)
    double precision :: xgradold(images*C_N), v(images*C_N)
    integer :: images, i
    double precision, parameter :: dt=1.d-2, friction=0.9d0

    !print *,trim(kind_opt)
    if (trim(kind_opt) .eq. 'verlet') then 
      v = 0.d0 
      xgrad = 0.d0
      do i=1,500 
        xpath = xpath + v*dt - xgrad*dt*dt/2.d0
        call dneb_grad(xpath,xgrad,epath,images)
        v = v + (-xgrad - xgradold)*dt/2.d0
        v = friction*v
        xgradold = xgrad 
        !print *,i,epath
      enddo
    endif

    return 
  end subroutine dneb_optimize

!...............................................................................
! creating a path of iamges between the start and end point of the chain.

  subroutine dneb_path(xst,xnd,xpath,images)
    implicit none 
    double precision :: xst(C_N), xnd(C_N), xpath(images*C_N), dvec(C_N)
    integer :: images, i, i1,in

    dvec = dneb_polarTau(xst,xnd) / DBLE(images-1)
    xpath = 0.d0
    do i=1,images
      i1 = (i-1)*C_N + 1
      in = (i-1)*C_N + C_N 
      xpath(i1:in) = xst + DBLE(i-1)*dvec
    enddo

    return 
  end subroutine dneb_path

!...............................................................................
! calculate the tangent at each image of the path.

  subroutine dneb_tau(xpath,xtau,eimg,images)
    implicit none 
    double precision :: xpath(images*C_N), xtau(images*C_N), eimg(images)
    double precision :: taup(C_N), taum(C_N), tau(C_N), deltaMax, deltaMin
    integer :: images, ip1, ipn, i1, in, im1, imn, i

! we store the start and end points for the images required to calculate the 
! tangent vectors.

    xtau = 0.d0
    do i=2,images-1 
      ip1 = i*C_N + 1 
      ipn = i*C_N + C_N
      i1 = (i-1)*C_N + 1
      in = (i-1)*C_N + C_N
      im1 = (i-2)*C_N + 1
      imn = (i-2)*C_N + C_N 

! we calculate taup = r(i+1) - r(i).

      taup = dneb_polarTau(xpath(ip1:ipn),xpath(i1:in))

! and we calculate taum = r(i) - r(i-1).

      taum = dneb_polarTau(xpath(i1:in),xpath(im1:imn))

! calculate tau according to reference [hen00]. if the image is neither a 
! minimum nor a maximum in energy, tau is chosen such that it points towards the 
! neighboring image having the higher energy.
! >>> tau = taup [E(i+1) > E(i) > E(i-1)]
! >>> tau = taum [E(i+1) < E(i) < E(i-1)]
! if the image is a minimum or a maximum along the path, a weighted average is 
! used to calculate the tangent vector.

      if ((eimg(i+1) .ge. eimg(i)) .and. (eimg(i) .ge. eimg(i-1))) then
        tau = taup
      elseif ((eimg(i+1) .le. eimg(i)) .and. (eimg(i) .le. eimg(i-1))) then
        tau = taum 
      else 
        deltaMax = MAX(ABS(eimg(i+1)-eimg(i)),ABS(eimg(i-1)-eimg(i)))
        deltaMin = MIN(ABS(eimg(i+1)-eimg(i)),ABS(eimg(i-1)-eimg(i)))
        if (eimg(i+1) .gt. eimg(i-1)) then 
          tau = taup*deltaMax + taum*deltaMin
        else 
          tau = taup*deltaMin + taum*deltaMax
        endif
      endif
      
! store tau in array.

      xtau(i1:in) = tau / NORM2(tau)
    enddo

    return 
  end subroutine dneb_tau

!...............................................................................
! helper function to calculate the angle between to polar vectors.

  function dneb_polarTau(x,y)
    implicit none 
    double precision :: dneb_polarTau(C_N), x(C_N), y(C_N), dangle
    integer :: i 

    do i=1,C_N 
      dangle = y(i) - x(i)
      if (dangle .lt. -C_PI) then 
        dangle = dangle + 2.d0*C_PI 
      elseif (dangle .gt. C_PI) then 
        dangle = dangle - 2.d0*C_PI 
      endif 
      dneb_polarTau(i) = dangle 
    enddo

  end function dneb_polarTau

!...............................................................................
! helper function to calculate the distance between two polar vectors.

  function dneb_distance(x,y)
    implicit none 
    double precision :: dneb_distance, x(C_N), y(C_N), dangle
    integer :: i 

    dneb_distance = 0.d0
    do i=1,C_N 
      dangle = ABS(x(i) - y(i))
      dangle = MIN(dangle, 2.d0*C_PI - dangle)

      dneb_distance = dneb_distance + dangle*dangle
    enddo 
    dneb_distance = SQRT(dneb_distance)

    return
  end function dneb_distance

!...............................................................................

  subroutine dneb_output(epath,xpath,images)
    implicit none
    double precision :: epath(images), xpath(images*C_N), xone(C_N), xtwo(C_N)
    double precision :: distance 
    integer :: i, images
    character(len=25) :: filename ='bums.data'
    
    !call f_getfilename('path.data',filename)
    open(unit=19,file=trim(filename),action='write',status='replace')

    write(19,'(2(F15.5))') 0.d0,epath(1)
    distance = 0.d0
    do i=2,images  
      xone = xpath((i-1)*C_N+1:(i-1)*C_N+C_N)
      xtwo = xpath((i-2)*C_N+1:(i-2)*C_N+C_N)

      distance = distance + potential_distance(xone,xtwo)
      write(19,'(2(F15.5))') distance,epath(i)
    enddo
    close(19)

    return
  end subroutine dneb_output

!...............................................................................

end module dneb
