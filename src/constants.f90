! constants.f90 
! WRAITH - david gallina - last update: 04/12/23
!...............................................................................
! module contains all global constants and variables that are used in the WRAITH
! code. moreover, a number of type structures are defined.
!...............................................................................

module constants
  implicit none 

! mathematical constants.

  double precision, parameter :: C_PI = 3.1415926535897932d0    ! pi 
  double precision, parameter :: C_EXP = 2.7182818284590452d0   ! euler 

! physical constants.

  double precision, parameter :: C_MU0 = 201335.4520789104d0    ! [nm3*mT/meV]
  double precision, parameter :: C_MUB = 5.7883818012d-5        ! [meV/mT]
  double precision, parameter :: C_KB = 8.617333262d-2          ! [eV/K]
  double precision, parameter :: C_GAMMA = 1.760859644d-1       ! [rad/mT*µs]
  double precision, parameter :: C_M = 2.89884d0                ! [meV/mT]

! global variables.

  double precision :: C_FIELD(3)                                ! [mT]
  double precision :: C_TEMP                                    ! [K]
  double precision :: C_ALPHA                                   ! 
  double precision :: C_TIME                                    ! [µs]
  double precision :: C_SCALE
  integer :: C_P                                                ! # particles
  integer :: C_N                                                ! # dimensions
  integer :: C_C                                                ! # cartesians
  integer :: C_S                                                ! # spherics
  integer :: C_SEED                                             ! rng seed

! potential variables.

  character(len=15) :: C_POT                                    ! potential
  logical :: DPLS1D_LOG                                         ! selector
  logical :: DPLS2D_LOG                                         ! selector
  logical :: OPEN1D_LOG                                         ! selector
  logical :: UNSURF_LOG                                         ! selector
  logical :: CHAINS_LOG                                         ! selector
  logical :: CONTINUE_LOG                                       ! continue 

! type definitions.

! nanoparticle type.

  type t_np
    double precision, dimension(3) :: xyz                       ! position [nm]
    double precision :: radius                                  ! radius [nm]
    double precision :: height                                  ! height [nm]
    double precision :: volume                                  ! volume [nm^3]
    double precision :: moment                                  ! moment 
    character(len=10) :: material
  end type t_np

! local minima type.

  type t_lm
    double precision, allocatable :: x(:)                       ! configuration 
    double precision :: m(3)                                    ! moment
    double precision :: v(3)                                    ! vorticity
    double precision :: e                                       ! energy
    double precision :: s                                       ! vib. entropy
    integer :: pos
  end type t_lm 

! transition state type.

  type t_ts  
    double precision, allocatable :: x(:)                       ! configuration
    double precision :: k(2)                                    ! trans. rates
    double precision :: e                                       ! energy
    double precision :: s                                       ! vib. entropy
    integer :: iloc(2)                                          ! adj. minima
  end type t_ts

! material type.

  type t_material 
    character(len=10) :: id 
    double precision :: magnetization, anisotropy, mass 
  end type t_material 

! material storage dummies.

  type(t_material), allocatable, private, save :: c_material_vec(:)
  integer, private, save :: num_materials = 0
  logical, private, save :: log_materials = .false.

contains 

!...............................................................................
! material function.

  function c_material(material,property)
    implicit none 
    character(len=20), intent(in) :: material 
    double precision :: c_material
    integer :: ierr, i, property

! if the routine is called for the first time, we read the parameters from file.

    if (.not. log_materials) then 
      open(unit=13,file='materials.data',action='read',status='old',iostat=ierr)
      if (ierr .ne. 0) then 
        print *,'                                '
        print *,'    > materials file is missing.'
        print *,'                                '
        stop 
      endif 

! get number of materials.

      do  
        read(13,*,iostat=ierr) 
        if (ierr .ne. 0) then 
          exit 
        else 
          num_materials = num_materials + 1
        endif 
      enddo 

! store materials.

      allocate(c_material_vec(num_materials))
      rewind(13)
      do i=1,num_materials 
        read(13,*) c_material_vec(i)%id,c_material_vec(i)%magnetization,&
                   c_material_vec(i)%anisotropy,c_material_vec(i)%mass 
      enddo
      close(13)
      log_materials = .true. 
    endif 

! get the material parameters from the dataset.

    do i=1,num_materials 
      if (trim(c_material_vec(i)%id) .eq. trim(material)) then 
        if (property .eq. 1) then
          c_material = c_material_vec(i)%magnetization
        elseif (property .eq. 2) then
          c_material = c_material_vec(i)%anisotropy
        elseif (property .eq. 3) then 
          c_material = c_material_vec(i)%mass
        endif 
      endif 
    enddo

    return 
  end function c_material

!...............................................................................

end module constants                                                            