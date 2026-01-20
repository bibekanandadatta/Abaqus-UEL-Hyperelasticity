      include '../UEL-NLMech/src/uel_nlmech_pk1.for'

! **********************************************************************
! **********************************************************************

!     using ifort compiler compile the code from intel oneAPI terminal:
!     ifort /Qmkl -o nlmech_pk1 test_nlmech_pk1.for
!     this command will generate an executable call nlmech_pk1
!     execute the pegel on the terminal using: .\nlmech_pk1

! **********************************************************************
! **********************************************************************

      PROGRAM NLMECHANICAL

      use global_parameters, only: wp, zero, one, two, half

      implicit none

      integer, parameter    :: case     = 1
      integer, parameter    :: fileUnit = 15
      character(len=256)    :: fileName

      real(wp), allocatable :: RHS(:,:), AMATRX(:,:), SVARS(:)
      real(wp), allocatable :: PROPS(:), COORDS(:,:), DUall(:,:)
      real(wp), allocatable :: Uall(:), Vel(:), Accn(:)
      real(wp), allocatable :: ADLMAG(:,:), PREDEF(:,:,:), DDLMAG(:,:)

      integer, allocatable  :: JPROPS(:)

      real(wp)  :: ENERGY(8), PNEWDT, TIME(2), DTIME, PARAMS(3), PERIOD

      integer   :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE
      integer   :: KSTEP, KINC, JELEM, NDLOAD, JDLTYP(1,1), NPREDF
      integer   :: MLVARX, MDLOAD, NJPROPS, LFLAGS(5), nlSdv, ngSdv

      integer   :: nDim, nStress, nOrder, nInt, matFlag, nPostVars
      integer   :: fbarFlag
      real(wp)  :: Gshear, kappa, lam_L

      real(wp)  :: elemLen, lam_1_App, lam_2_App, u_1_App, u_2_App
      integer   :: i, j, k, l

      ! additional vectors and matrices for numerical derivatives
      real(wp), allocatable :: Uall_1(:), Uall_2(:)
      real(wp), allocatable :: RHS_1(:,:), AMATRX_1(:,:)
      real(wp), allocatable :: RHS_2(:,:), AMATRX_2(:,:)
      real(wp), allocatable :: AMATRX_h(:,:)
      real(wp)              :: delta_h



      ! default boundary conditions
      lam_1_App   = one

      if (case.eq.1) then
        write(*,*) 'testing case-I'
        fileName = 'NLmechDebug_1_PK1.dat'
      else if (case.eq.2) THEN
        write(*,*) 'testing case-II'
        fileName = 'NLmechDebug_2_PK1.dat'
        lam_1_App = two
      else
        write(*,*) 'wrong case'
        stop
      end if

      elemLen   = 1.0e-3_wp
      lam_2_App = one/lam_1_App
      u_1_App   = (lam_1_App - one)*elemLen
      u_2_App   = (lam_2_App - one)*elemLen




      ! element type, element no, dimension, etc.
      JELEM   = 1        ! no of elements in patch test
      nDim    = 2        ! spatial dimension of the problem
      nOrder  = 1        ! polynomial order of the lagrangian element
      JTYPE   = 7        ! type of element (U3: HEX8, U7: QUAD4)

      if (nDim .eq. 3) then
        nStress = 6
      else if (nDim .eq. 2) then
        nStress = 3
      else
        write(*,*) 'model dimension is wrong'
        stop
      end if

      ! nInt : no of volume integration point, nPostVars: no of post-processed variables
      if ((nDim .eq. 3) .and. (nOrder .eq. 1)) then
        nNode     = 8
        nInt      = 8
        nPostVars = 12
      else if ((nDim .eq. 2) .and. (nOrder .eq. 1)) then
        nNode     = 4
        nInt      = 4
        nPostVars = 6
      else
        write(*,*)  'model dimension or element order is wrong'
        stop
      end if


      ! Abaqus variables
      NDOFEL  = nNode*nDim      ! no of total degrees of freedom
      NRHS    = 1               ! right hand side residual column dimension
      NSVARS  = 1               ! no of state variables
      NPROPS  = 3               ! total no of material properties
      MCRD    = nDim            ! no of dimension (ABAQUS does it differently)
      NJPROPS = 4               ! no of integer properties
      NPREDF  = 1               ! no of predefined field
      MLVARX  = NDOFEL          ! maximum variables (assuming same as NDOFEL)


      allocate(RHS(MLVARX,NRHS), AMATRX(NDOFEL,NDOFEL), SVARS(NSVARS),
     &      PROPS(NPROPS), COORDS(MCRD,NNODE), DUALL(NDOFEL,1),
     &      UALL(NDOFEL), Vel(NDOFEL), Accn(NDOFEL), ADLMAG(1,1),
     &      PREDEF(2,NPREDF,1), DDLMAG(NNODE,1), JPROPS(NJPROPS))

      ! initializing the output of UEL subroutine
      RHS     = zero            ! residual vector
      AMATRX  = zero            ! stiffness matrix
      SVARS   = zero            ! array containing state variables
      ENERGY  = zero            ! array containing energy
      PNEWDT  = zero            ! time stepping flag

      ! initialize some other not-so-necessary (for statics) input to UEL
      LFLAGS  = [1, 1, 0, 0, 0]   ! Abaqus step flag
      VEL     = zero              ! velocity
      ACCN    = zero              ! acceleration
      PARAMS  = zero              ! time parameters (irrelevant now)
      PREDEF  = zero              ! no predefined field

      ! no distributed load
      NDLOAD = zero
      MDLOAD = zero

      KSTEP   = 1
      KINC    = 1
      PERIOD  = 1.0e-2_wp
      TIME    = 1.0e-2_wp
      DTIME   = 1.0e-2_wp


      ! material properties
      Gshear      = 100.0e3_wp        ! Shear modulus
      kappa       = 100.0_wp*Gshear   ! Bulk modulus
      lam_L       = 0.0_wp            ! Locking stretch (only for AB model)

      fbarFlag    = 0
      matFlag     = 1

      ! define material and element properties
      PROPS(1:NPROPS)   = [ Gshear, kappa, lam_L ]
      JPROPS(1:NJPROPS) = [ nInt, fbarFlag, matFlag, nPostVars ]



      ! define the element coordinates and "simulated" nodal solution
      if (nDim.eq.2) then
        COORDS(1,:) = [0, 1, 1, 0]
        COORDS(2,:) = [0, 0, 1, 1]
        Uall(1:nDOFEL) = [zero,     zero,
     &                    u_1_App,  zero,
     &                    u_1_App,  u_2_App,
     &                    zero,     u_2_App]

      else if (nDim.eq.3) then
        COORDS(1,:) = [0, 1, 1, 0, 0, 1, 1, 0]
        COORDS(2,:) = [0, 0, 1, 1, 0, 0, 1, 1]
        COORDS(3,:) = [0, 0, 0, 0, 1, 1, 1, 1]
        Uall(1:nDOFEL) = [zero, zero, zero,
     &                    zero, zero, zero,
     &                    zero, zero, zero,
     &                    zero, zero, zero,
     &                    zero, zero, zero,
     &                    zero, zero, zero,
     &                    zero, zero, zero,
     &                    zero, zero, zero]
      else
        write(*,*)  'model dimension is wrong', nDim
        stop
      end if

      COORDS = COORDS*elemLen



      ! open the debugging file
      open(unit=fileUnit, file=fileName, status='unknown')


      ! call the UEL subroutine
      call UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)


      !! numerical tangent matrix calculation
      allocate( UALL_1(NDOFEL), RHS_1(MLVARX,NRHS),
     &          UALL_2(NDOFEL), RHS_2(MLVARX,NRHS),
     &          AMATRX_1(NDOFEL,NDOFEL), AMATRX_2(NDOFEL,NDOFEL),
     &          AMATRX_h(NDOFEL,NDOFEL) )


      !! perturb the degrees of freedom and calculate the tangent matrix
      delta_h   = sqrt(epsilon(one))

      do i = 1, ndofel
        do j = 1, ndofel

          ! perturb in +ve direction
          RHS_1     = zero
          AMATRX_1  = zero
          SVARS     = zero
          Uall_1    = Uall
          Uall_1(j) = Uall_1(j) + delta_h

          call UEL(RHS_1,AMATRX_1,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &      PROPS,NPROPS,COORDS,MCRD,NNODE,Uall_1,DUall,Vel,Accn,JTYPE,
     &      TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     &      PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     &      NJPROPS,PERIOD)

          ! perturb in -ve direction
          RHS_2     = zero
          AMATRX_2  = zero
          SVARS     = zero
          Uall_2    = Uall
          Uall_2(j) = Uall_2(j) - delta_h

          call UEL(RHS_2,AMATRX_2,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &      PROPS,NPROPS,COORDS,MCRD,NNODE,Uall_2,DUall,Vel,Accn,JTYPE,
     &      TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     &      PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     &      NJPROPS,PERIOD)

          AMATRX_h(i,j) = - (RHS_1(i,1) - RHS_2(i,1))/(2*delta_h)
        end do
      end do


      write(fileUnit,'(A)') 'RHS: '
      do i = 1,NDOFEL
        write(fileUnit,*) RHS(i,1)
      enddo

      write(fileUnit,'(A)') 'AMATRX: '
      do i = 1, NDOFEL
        write(fileUnit, '(*(F24.14))') AMATRX(i,:)
      end do

      write(fileUnit,'(A)') 'AMATRX_h: '
      do i = 1, NDOFEL
        write(fileUnit, '(*(F24.14))') AMATRX_h(i,:)
      end do

      ! close the debugging file
      close(unit=fileUnit)

      END PROGRAM NLMECHANICAL

! **********************************************************************

      ! theses subroutines emulate the utility subroutines from ABAQUS
      SUBROUTINE XIT
        stop
      END SUBROUTINE XIT

      SUBROUTINE GETJOBNAME(jobName,lenJobName)

        integer, intent(inout)            :: lenJobName
        character(len=256), intent(inout) :: jobName

        RETURN
      END SUBROUTINE GETJOBNAME

      SUBROUTINE GETOUTDIR(outDir,lenOutDir)

        integer, intent(inout)            :: lenOutDir
        character(len=256), intent(inout) :: outDir

        RETURN
      END SUBROUTINE GETOUTDIR

! **********************************************************************