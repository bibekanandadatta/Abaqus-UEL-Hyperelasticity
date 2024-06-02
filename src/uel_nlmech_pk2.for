! **********************************************************************
! ********* ABAQUS/ STANDARD USER ELEMENT SUBROUTINE (UEL) *************
! **********************************************************************
* large strain displacement element + Neo-Hookean & Arruda-Boyce model *
! **********************************************************************
!                     BIBEKANANDA DATTA (C) MAY 2024
!                 JOHNS HOPKINS UNIVERSITY, BALTIMORE, MD
! **********************************************************************
! **********************************************************************
!                       JTYPE DEFINITION
!
!     U1                THREE-DIMENSIONAL TET4 ELEMENT
!     U2                THREE-DIMENSIONAL TET10 ELEMENT
!     U3                THREE-DIMENSIONAL HEX8 ELEMENT
!     U4                THREE-DIMENSIONAL HEX20 ELEMENT
!
!     U5                PLANE STRAIN TRI3 ELEMENT
!     U6                PLANE STRAIN TRI6 ELEMENT
!     U7                PLANE STRAIN QUAD4 ELEMENT
!     U8                PLANE STRAIN QUAD8 ELEMENT
!
! **********************************************************************
!          VOIGT NOTATION CONVENTION FOR STRESS/ STRAIN TENSORS
!
!       In this subroutine we adopted the following convention for
!       symmetric stress and strain tensor following Voigt notation
!       This is different than what is followed by Abaqus/ Standard
!
!         sigma11, sigma22, sigma33, sigma23, sigma13, sigma12
!       strain11, strain22, strain33, strain23, strain13, strain12
!
! **********************************************************************
!                       LIST OF MATERIAL PROPERTIES
!
!     G           = props(1)        Shear modulus
!     kappa       = props(2)        Bulk modulus
!     lam_L       = props(3)        Locking stretch for AB model (only)
!
! **********************************************************************
!                        LIST OF ELEMENT PROPERTIES
!
!     jprops(1)   = nInt            no of integration points in element
!     jprops(2)   = matID           constitutive relation for material
!     jprops(3)   = nPostVars       no of local (int pt) post-processing variables
!
! **********************************************************************
!                        POST-PROCESSED VARIABLES
!                     (follows the convention above)
!
!     uvar(1:nStress)                 Cauchy stress tensor components
!     uvar(nStress+1:2*nStress)         Euler strain tensor components
!
! **********************************************************************
!               VARIABLES TO BE UPDATED WITHIN THE SUBROUTNE
!
!     RHS(i)                        Right hand side vector.
!     AMATRX(i,j)                   Stiffness matrix
!     SVARS(1:NSVARS)               Element state variables.  Must be updated in this routine
!     ENERGY(1:8)                   Energy(1) Kinetic Energy
!                                   Energy(2) Elastic Strain Energy
!                                   Energy(3) Creep Dissipation
!                                   Energy(4) Plastic Dissipation
!                                   Energy(5) Viscous Dissipation
!                                   Energy(6) Artificial strain energy
!                                   Energy(7) Electrostatic energy
!                                   Energy(8) Incremental work done by loads applied to the element
!     PNEWDT                        Allows user to control ABAQUS time increments.
!                                   If PNEWDT<1 then time step is abandoned and computation is restarted with
!                                   a time increment equal to PNEWDT*DTIME
!                                   If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
!
!                       VARIABLES PROVIDED FOR INFORMATION
!
!     NDOFEL                        Total # DOF for the element
!     NRHS                          Dimension variable
!     NSVARS                        Total # element state variables
!     PROPS(1:NPROPS)               User-specified properties of the element
!     NPROPS                        No. properties
!     JPROPS(1:NJPROPS)             Integer valued user specified properties for the element
!     NJPROPS                       No. integer valued properties
!     COORDS(i,N)                   ith coordinate of Nth node on element
!     MCRD                          Maximum of (# coords,minimum of (3,#DOF)) on any node
!     Uall                          Vector of DOF at the end of the increment
!     DUall                         Vector of DOF increments
!     Vel                           Vector of velocities (defined only for implicit dynamics)
!     Accn                          Vector of accelerations (defined only for implicit dynamics)
!     JTYPE                         Integer identifying element type (the number n in the Un specification in the input file)
!     TIME(1:2)                     TIME(1)   Current value of step time
!                                   TIME(2)   Total time
!     DTIME                         Time increment
!     KSTEP                         Current step number
!     KINC                          Increment number
!     JELEM                         User assigned element number in ABAQUS
!     PARAMS(1:3)                   Time increment parameters alpha, beta, gamma for implicit dynamics
!     NDLOAD                        Number of user-defined distributed loads defined for this element
!     JDLTYP(1:NDLOAD)              Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
!     ADLMAG(1:NDLOAD)              Distributed load magnitudes
!     DDLMAG(1:NDLOAD)              Increment in distributed load magnitudes
!     PREDEF(1:2,1:NPREDF,1:NNODE)  Predefined fields.
!     PREDEF(1,...)                 Value of predefined field
!     PREDEF(2,...)                 Increment in predefined field
!     PREDEF(1:2,1,k)               Value of temperature/temperature increment at kth node
!     PREDEF(1:2,2:NPREDF,k)        Value of user defined field/field increment at kth node
!     NPREDF                        Number of predefined fields
!     LFLAGS                        Load type control variable
!     MLVARX                        Dimension variable
!     PERIOD                        Time period of the current step
! **********************************************************************

      ! make sure to have the correct directory
      include '../module/global_parameters.for'     ! parameter module
      include '../module/error_logging.for'         ! error/ debugging module
      include '../module/linear_algebra.for'        ! linear algebra module
      include '../module/nonlinear_solver.for'      ! Newton-Raphson solver module
      include '../module/lagrange_element.for'      ! module for Lagrange elements
      include '../module/gauss_quadrature.for'      ! Guassian quadrature module
      include '../module/solid_mechanics.for'       ! solid mechanics module
      include '../module/post_processing.for'       ! post-processing module

! **********************************************************************
! ****************** ABAQUS USER ELEMENT SUBROUTINE ********************
! **********************************************************************

      subroutine UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)

      use global_parameters
      use error_logging
      use post_processing

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     & DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      ! type specification of UEL arguments
      integer             :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD
      integer             :: NNODE, JTYPE, KSTEP, KINC, JELEM
      integer             :: NDLOAD, JDLTYP, NPREDF, LFLAGS
      integer             :: MLVARX, MDLOAD, JPROPS, NJPROPS

      real(kind=wp)       :: RHS, AMATRX, SVARS, ENERGY, PNEWDT
      real(kind=wp)       :: PROPS, COORDS, DUall, Uall, Vel, Accn
      real(kind=wp)       :: TIME, DTIME, PARAMS, ADLMAG, PREDEF
      real(kind=wp)       :: DDLMAG, PERIOD

      integer             :: nInt, matID, nPostVars
      integer             :: nDim, nStress, uDOF, uDOFEL

      integer             :: lenJobName,lenOutDir
      character(len=256)  :: outDir
      character(len=256)  :: jobName
      character(len=512)  :: errFile, dbgFile
      character(len=8)    :: analysis, abqProcedure
      logical             :: nlgeom
      type(logger)        :: msg


      !! open a debug file for the current job
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      errFile = trim(outDir)//'\aaERR_'//trim(jobName)//'.dat'
      dbgFile = trim(outDir)//'\aaDBG_'//trim(jobName)//'.dat'
      call msg%fopen( errfile=errFile, dbgfile=dbgFile )


      !! change the LFLAGS criteria as needed (check abaqus UEL manual)
      if((lflags(1).eq.1).or.(lflags(1).eq.2)) then
        abqProcedure = 'STATIC'
      else
       call msg%ferror(flag=error, src='UEL',
     &                msg='Incorrect Abaqus procedure.', ia=lflags(1))
        call xit
      end if

      !! check if the procedure is linear or nonlinear
      if (lflags(2) .eq. 0) then
        nlgeom = .false.
      else if (lflags(2) .eq. 1) then
        nlgeom = .true.
      end if

      ! check to see if it's a general step or a linear purturbation step
      if(lflags(4).eq.1) then
         call msg%ferror(flag=error, src='UEL',
     &      msg='The step should be a GENERAL step.', ia=lflags(4))
        call xit
      end if


      !! assign parameter specific to analysis and element types
      if ((JTYPE.ge.1).and.(JTYPE.le.4)) then
        nDim      = 3
        analysis  = '3D'            ! three-dimensional analysis
        nStress   = 6
        uDOF      = nDim            ! displacement degrees of freedom of a node
        uDOFEL    = nNode*uDOF      ! total displacement degrees of freedom in element
      else if ((JTYPE.ge.5).and.(JTYPE.le.8)) then
        nDim      = 2
        analysis  = 'PE'            ! plane strain analysis
        nStress   = 3
        uDOF      = nDim            ! displacement degrees of freedom of a node
        uDOFEL    = nNode*uDOF      ! total displacement degrees of freedom in element
      else
        call msg%ferror( flag=error, src='UEL',
     &            msg='Element is unavailable.', ia=jtype )
        call xit
      end if


      nInt      = jprops(1)
      matID     = jprops(2)
      nPostVars = jprops(3)


      ! array containing variables for post-processing
      if (.not. allocated(globalPostVars)) then
        allocate(globalPostVars(numElem,nInt,nPostVars))

        ! print job-related information the first time
        call msg%finfo('---------------------------------------')
        call msg%finfo('------- ABAQUS FINITE STRAIN UEL ------')
        call msg%finfo('---------------------------------------')
        call msg%finfo('--- Abaqus Job: ', ch=trim(jobName))
        call msg%finfo('---------------------------------------')
        call msg%finfo('------- PROCEDURE       = ', ch=abqProcedure)
        call msg%finfo('------- ANALYSIS TYPE   = ', ch=analysis)
        call msg%finfo('---------- NLGEOM       = ', la=nlgeom)
        call msg%finfo('------- MODEL DIMENSION = ', ia=nDim)
        call msg%finfo('------- ELEMENT NODES   = ', ia=nNode)
        call msg%finfo('---------------------------------------')
        call msg%finfo('-------- INTEGRATION SCHEME -----------')
        call msg%finfo('----------- NINT   = ', ia=nInt)
        call msg%finfo('---------------------------------------')
        call msg%finfo('---------- POST-PROCESSING ------------')
        call msg%finfo('--- NO OF ELEMENTS            = ', ia=numElem)
        call msg%finfo('--- OVERLAY ELEMENT OFFSET    = ',ia=elemOffset)
        call msg%finfo('--- NO OF VARIABLES AT INT PT = ', ia=nPostVars)
        call msg%finfo('---------------------------------------')
      end if

       ! call your UEL subroutine
       call uelNLMECH(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL)


      end subroutine UEL

! **********************************************************************
! **********************************************************************

      subroutine uelNLMECH(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL)

      use global_parameters
      use error_logging
      use linear_algebra
      use lagrange_element
      use gauss_quadrature
      use solid_mechanics
      use post_processing

      implicit none

     !!!!!!!!!!!!! VARIABLE DECLARATION AND INITIALIZATION !!!!!!!!!!!!!

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     & DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      !! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      integer           :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD
      integer           :: NNODE, JTYPE, KSTEP, KINC, JELEM
      integer           :: NDLOAD, JDLTYP, NPREDF, LFLAGS
      integer           :: MLVARX, MDLOAD, JPROPS, NJPROPS

      real(kind=wp)     :: RHS, AMATRX, SVARS, ENERGY, PNEWDT
      real(kind=wp)     :: PROPS, COORDS, DUall, Uall, Vel, Accn, TIME
      real(kind=wp)     :: DTIME, PARAMS, ADLMAG, PREDEF, DDLMAG, PERIOD

      character(len=8)  :: analysis
      integer           :: nDim, nStress, nInt
      integer           :: uDOF, uDOFEL, nlocalSdv
      logical           :: nlgeom

      real(kind=wp)     :: uNode(nDim,nNode), duNode(nDim,nNode), F(3,3)
      real(kind=wp)     :: ID(nDim,nDim), w(nInt), xi(nInt,nDim)
      real(kind=wp)     :: Nxi(nNode), dNdxi(nNode,nDim)
      real(kind=wp)     :: dxdxi(nDim,nDim), dxidx(nDim,nDim)
      real(kind=wp)     :: dNdx(nNode,nDim), detJ
      real(kind=wp)     :: Na(nDim,nDim), Nmat(nDim,uDOFEl)
      real(kind=wp)     :: Ba(nStress,nDim), Bmat(nStress,uDOFEl)
      real(kind=wp)     :: Ga(nDim*nDim,nDim), Gmat(nDim*nDim,uDOFEl)
      real(kind=wp)     :: fieldNode(npredf,nNode)
      real(kind=wp)     :: dfieldNode(npredf,nNode), fieldVar(npredf)
      real(kind=wp)     :: dfieldVar(npredf)
      real(kind=wp)     :: stressTensorPK2(nStress,nStress)
      real(kind=wp)     :: SIGMA_S(nDim**2,nDim**2)
      real(kind=wp)     :: SIGMA_F(nDim*nNode,nDim*nNode)
      real(kind=wp)     :: strainLagrange(nStress,1)
      real(kind=wp)     :: strainEuler(nStress,1)
      real(kind=wp)     :: stressCauchy(nStress,1)
      real(kind=wp)     :: stressPK1(nDim*nDim,1)
      real(kind=wp)     :: stressPK2(nStress,1), Dmat(nStress,nStress)
      real(kind=wp)     :: Kuu(uDOFEl,uDOFEl), Ru(uDOFEl,1)

      integer           :: i, j, intPt, matID
      type(element)     :: solidFiniteStrain
      type(logger)      :: msg

      !! set the element parameters
      solidFiniteStrain = element(nDim=nDim,analysis=analysis,
     &                            nNode=nNode,nInt=nInt)

      !! initialize the matrices and vectors
      F       = zero
      Na      = zero
      Ba      = zero
      Ga      = zero
      Nmat    = zero
      Bmat    = zero
      Gmat    = zero
      Dmat    = zero
      SIGMA_F = zero
      SIGMA_S = zero
      Kuu     = zero
      Ru      = zero

      ENERGY                    = zero
      AMATRX(1:NDOFEL,1:NDOFEL) = zero
      RHS(1:MLVARX,1)           = zero

      matID = jprops(2)

     !!!!!!!!!!! END VARIABLE DECLARATION AND INITIALIZATION !!!!!!!!!!!

      ! reshape the displacement vectors into matrix forms
      uNode  = reshape(UAll,[nDim,nNode])
      duNode = reshape(DUAll(:,1),[nDim,nNode])

      ! if applicable gather the prescribed field variables in a matrix
      ! such as temperature/ something (as shown below - not yet tested)
      ! do k = 1 , npredf
      !   fieldNode(k,1:nNode) = predef(1,k,1:nNode)
      !   dfieldNode(k,1:nNode) = predef(2,k,1:nNode)
      ! end do

     !!!!!!!!!!!!!!!!! ELEMENT RELATED OPERATIONS !!!!!!!!!!!!!!!!!!!!!!

      ! obtain gauss quadrature points and weights
      call eyeMat(ID)
      call getGaussQuadrtr(solidFiniteStrain,w,xi)

      ! loop through all the integration points (main/ external loop)
      do intPt = 1, nInt

        call evalInterpFunc(solidFiniteStrain,xi(intPt,:),Nxi,dNdxi)

        ! calculate element jacobian and global shape func gradient
        dXdxi = matmul(coords,dNdxi)      ! calculate dXdxi
        call detMat(dXdxi,detJ)
        call inverse(dXdxi,dxidX)         ! calculate inverse
        dNdX  = matmul(dNdxi,dxidX)       ! calculate dNdX

        if (detJ .le. 0) then
          call msg%ferror( flag=warn, src='uelNLMech',
     &     msg='Negative element jacobian.', ivec=[jelem, intpt])
        end if

        ! loop over all the nodes (internal loop)
        do i=1,nNode

          ! form the nodal-level matrices: [Na], [Ga], [Ba]
          do j = 1, nDim
            Na(j,j) = Nxi(i)
            Ga(nDim*(j-1)+1:nDim*j,1:nDim) = dNdX(i,j)*ID
          end do

          ! form [Ba] matrix: plane stress/ plane strain case
          if (analysis.eq.'PE') then
            Ba(1,1)       = dNdx(i,1)
            Ba(2,2)       = dNdx(i,2)
            Ba(3,1:nDim)  = [dNdx(i,2), dNdx(i,1)]

          ! form [Ba] matrix: 3D case
          else if (analysis.eq.'3D') then
            Ba(1,1)       = dNdx(i,1)
            Ba(2,2)       = dNdx(i,2)
            Ba(3,3)       = dNdx(i,3)
            Ba(4,1:nDim)  = [  zero,      dNdx(i,3),  dNdx(i,2)]
            Ba(5,1:nDim)  = [dNdx(i,3),     zero,     dNdx(i,1)]
            Ba(6,1:nDim)  = [dNdx(i,2),   dNdx(i,1),    zero   ]

          else
            call msg%ferror( flag=error, src='uelMech',
     &                  msg='Wrong analysis.', ch=analysis )
            call xit
          end if

          ! form the [N], [B], and [G] matrix
          Nmat(1:nDim,nDim*(i-1)+1:nDim*i)    = Na(1:nDim,1:nDim)
          Bmat(1:nStress,nDim*(i-1)+1:nDim*i) = Ba(1:nStress,1:nDim)
          Gmat(1:nDim**2,nDim*(i-1)+1:nDim*i) = Ga(1:nDim**2,1:nDim)
        end do                             ! end of nodal point loop

      !!!!!!!!!!!!!! COMPLETE ELEMENT RELATED OPERATIONS !!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!! CONSTITUTIVE MODEL !!!!!!!!!!!!!!!!!!!!!!!

        ! calculate deformation gradient and deformation tensors
        F(1:nDim,1:nDim) = ID + matmul(uNode,dNdX)

        if (analysis .eq. 'PE') then
          F(3,3) = one
        end if

        !! interpolate the field variables at the integration point
        !! (this not yet tested or used)
    !     do k = 1, npredf
    !       fieldVar(k)   = dot_product( Nxi, 
    !  &                    reshape( fieldNode(k,1:nNode), [nNode] ) )
    !       dfieldVar(k)  = dot_product( Nxi, 
    !  &                    reshape( dfieldNode(k,1:nNode), [nNode] ) )
    !     end do


        ! call material point subroutine (UMAT) for specific material
        if (matID .eq. 1) then
          call umatNeoHookean(stressCauchy,stressPK1,stressPK2,
     &          Dmat,F,svars,nsvars,strainLagrange,strainEuler,time,
     &          dtime,fieldVar,npredf,nDim,analysis,nStress,jelem,intpt,
     &          coords,nnode,kstep,kinc,props,nprops,jprops,njprops,
     &          nlocalSdv)

        else if (matID .eq. 2) then
          call umatArrudaBoyce(stressCauchy,stressPK1,stressPK2,
     &          Dmat,F,svars,nsvars,strainLagrange,strainEuler,time,
     &          dtime,fieldVar,npredf,nDim,analysis,nStress,jelem,intpt,
     &          coords,nnode,kstep,kinc,props,nprops,jprops,njprops,
     &          nlocalSdv)
        end if
        ! can add more constitutive models using else if construct here

      !!!!!!!!!!!!!!!!!!!! END CONSTITUTIVE MODEL !!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!! TANGENT MATRIX AND RESIDUAL VECTOR !!!!!!!!!!!!!!!!!
        call vector2symtensor(stressPK2, stressTensorPK2)

        ! form the [SIGMA_S] matrix for geometric stiffness
        do i = 1, nDim
          do j = 1, nDim
              SIGMA_S(nDim*(i-1)+1:nDim*i,nDim*(j-1)+1:nDim*j) =
     &                    stressTensorPK2(i,j)*ID
          end do
        end do

        ! form [SIGMA_F] matrix for material stiffness
        do i = 1, nNode
          do j = 1, nNode
              if (i.eq.j) then                    ! banded diagonal
                  SIGMA_F(nDim*(i-1)+1:nDim*i,nDim*(j-1)+1:nDim*j)
     &                              = F(1:nDim,1:nDim)
              end if
          end do
        end do

        ! form the stiffness matrix and residual vector
        Kuu = Kuu + w(intpt)*detJ*
     &      ( matmul( transpose(matmul(Bmat,SIGMA_F)),
     &      matmul (Dmat,matmul(Bmat,SIGMA_F)) ) +
     &      matmul( transpose(Gmat), matmul(SIGMA_S,Gmat)) )

        Ru  = Ru - w(intpt)*detJ*
     &      matmul( transpose(matmul(Bmat,SIGMA_F)), stressPK2 )

      !!!!!!!!!!!!!! TANGENT MATRIX AND RESIDUAL VECTOR !!!!!!!!!!!!!!!!

      end do                           ! end of integration point loop

      ! body force and surface load can be added using dummy elements

      ! assign the element stiffness matrix to abaqus-defined variable
      AMATRX(1:NDOFEL,1:NDOFEL) = Kuu(1:uDOFEl,1:uDOFEl)
      RHS(1:MLVARX,1)           = Ru(1:uDOFEl,1)

      end subroutine uelNLMECH

! **********************************************************************
! **********************************************************************

      subroutine umatNeoHookean(stressCauchy,stressPK1,stressPK2,
     &      Dmat,F,svars,nsvars,strainLagrange,strainEuler,time,
     &      dtime,fieldVar,npredf,nDim,analysis,nStress,jelem,intpt,
     &      coords,nnode,kstep,kinc,props,nprops,jprops,njprops,
     &      nlocalSdv)

      use global_parameters
      use error_logging
      use linear_algebra
      use lagrange_element
      use solid_mechanics
      use post_processing

      implicit none

      integer             :: nsvars, npredf, nDim, nStress
      integer             :: jelem, intpt, nNode, kstep, kinc
      integer             :: nprops, njprops

      real(kind=wp)       :: F(3,3), stressCauchy(nStress,1)
      real(kind=wp)       :: stressPK1(nDim*nDim,1)
      real(kind=wp)       :: stressPK2(nStress,1)
      real(kind=wp)       :: Dmat(nStress,nStress)
      real(kind=wp)       :: strainLagrange(nStress,1)
      real(kind=wp)       :: strainEuler(nStress,1)
      real(kind=wp)       :: props(1:nprops), svars(1:nsvars), time(2)
      real(kind=wp)       :: dtime, coords(nDim,nNode), fieldVar(npredf)

      integer             :: jprops(1:njprops)
      character(len=8)    :: analysis

      real(kind=wp)       :: detF, C(3,3), Cinv(3,3), detC, B(3,3)
      real(kind=wp)       :: Binv(3,3), detB, straintensorEuler(3,3)
      real(kind=wp)       :: stressTensorPK2(3,3)
      real(kind=wp)       :: stressTensorCauchy(3,3), Cmat(3,3,3,3)
      real(kind=wp)       :: VoigtMat(nSymm,nSymm)
      real(kind=wp)       :: strainVoigtEuler(nSymm,1)
      real(kind=wp)       :: stressVoigtPK1(nUnsymmm,1)
      real(kind=wp)       :: stressVoigtPK2(nSymm,1)
      real(kind=wp)       :: stressVoigtCauchy(nSymm,1)
      real(kind=wp)       :: Gshear, kappa, lam_L

      integer             :: nInt, nlocalSdv
      integer             :: i, j, k, l
      type(logger)        :: msg

      ! initialize matrial stiffness tensors
      Cmat      = zero
      Dmat      = zero

      ! assign material properties to variables
      Gshear    = props(1)        ! Shear modulus
      kappa     = props(2)        ! Bulk modulus
      lam_L     = props(3)        ! locking stretch

      nInt      = jprops(1)
      nlocalSdv = NSVARS/nInt

      ! locking stretch should be infinity (0 as input) for NH model
      if (lam_L .ne. zero) then
        call msg%ferror(flag=error, src='umatNeohookean',
     &       msg='Incorrect material parameter (lam_L).', ra=lam_l)
        call xit
      end if


      ! perform all the constitutitve relations in 3D
      call detMat(F,detF)

      if (detF .le. zero) then
        call msg%ferror(flag=error, src='umatNeoHookean',
     &          msg='Issue with volume change (detF)', 
     &          ivec=[jelem, intpt], ra= detF)
        call xit
      end if

      B = matmul(F,transpose(F))
      C = matmul(transpose(F),F)

      call inverse(B,Binv)
      call inverse(C,Cinv)

      ! calculate Euler-Almansi strain tensor
      straintensorEuler = half*(ID3-Binv)

      ! calculate material tangent, C_ijkl
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              Cmat(i,j,k,l) = Cmat(i,j,k,l)
     &            + kappa*Cinv(i,j)*Cinv(k,l)
     &            + (Gshear-kappa*log(detF)) *( Cinv(i,k)*Cinv(j,l)
     &            + Cinv(i,l)*Cinv(j,k) )
            end do
          end do
        end do
      end do

      ! calculate stress tensors
      stressTensorCauchy  = (one/detF)*(Gshear*(B-ID3) +
     & 											kappa*log(detF)*ID3)
      stressTensorPK2     = Gshear*(ID3-Cinv) + kappa*log(detF)*Cinv


      ! transforms the stiffness tensor 3x3x3x3 to a 6x6 matrix
      call symtangent2matrix(Cmat,VoigtMat)

      ! transform the stress tensor (3x3) to Voigt vector form (6x1)
      call symtensor2vector(straintensorEuler, strainVoigtEuler)
      call symtensor2vector(stressTensorCauchy, stressVoigtCauchy)
      call symtensor2vector(stressTensorCauchy, stressVoigtPK2)

     !!!!!!!!!!!!!!! END OF CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!!


      ! reshape the Voigt matrix and tensor based on analysis
      if (analysis .eq. 'PE')  then
        Dmat(1:ndim,1:ndim)     = VoigtMat(1:ndim,1:ndim)
        Dmat(1:ndim,nStress)    = VoigtMat(1:ndim,nSymm)
        Dmat(nStress,1:ndim)    = VoigtMat(nSymm,1:ndim)
        Dmat(nStress,nStress)   = VoigtMat(nSymm,nSymm)

        strainEuler(1:ndim,1)   = strainVoigtEuler(1:ndim,1)
        strainEuler(nStress,1)  = strainVoigtEuler(nSymm,1)

        stressPK2(1:nDim,1)     = stressVoigtPK2(1:nDim,1)
        stressPK2(nStress,1)    = stressVoigtPK2(nSymm,1)

        stressCauchy(1:nDim,1)  = stressVoigtCauchy(1:nDim,1)
        stressCauchy(nStress,1) = stressVoigtCauchy(nSymm,1)

      else if (analysis .eq. '3D') then
        Dmat = VoigtMat
        strainEuler   = strainVoigtEuler
        stressPK2     = stressVoigtPK2
        stressCauchy  = stressVoigtCauchy

      else
        call msg%ferror(flag=error, src='umatNeohookean',
     &            msg='Wrong analysis.', ch=analysis)
        call xit
      end if

      ! save the variables to be post-processed in globalPostVars
      globalPostVars(jelem,intpt,1:nStress) = stressCauchy(1:nStress,1)
      globalPostVars(jelem,intpt,nStress+1:2*nStress) =
     &                     strainEuler(1:nStress,1)

      end subroutine umatNeoHookean

! **********************************************************************

      subroutine umatArrudaBoyce(stressCauchy,stressPK1,stressPK2,
     &      Dmat,F,svars,nsvars,strainLagrange,strainEuler,time,
     &      dtime,fieldVar,npredf,nDim,analysis,nStress,jelem,intpt,
     &      coords,nnode,kstep,kinc,props,nprops,jprops,njprops,
     &      nlocalSdv)

      use global_parameters
      use linear_algebra
      use lagrange_element
      use solid_mechanics
      use post_processing
      use error_logging

      implicit none

      integer             :: nsvars, npredf, nDim, nStress
      integer             :: jelem, intpt, nNode, kstep, kinc
      integer             :: nprops, njprops

      real(kind=wp)       :: F(3,3), stressCauchy(nStress,1)
      real(kind=wp)       :: stressPK1(nDim*nDim,1)
      real(kind=wp)       :: stressPK2(nStress,1), Dmat(nStress,nStress)
      real(kind=wp)       :: strainLagrange(nStress,1)
      real(kind=wp)       :: strainEuler(nStress,1)
      real(kind=wp)       :: props(1:nprops), svars(1:nsvars), time(2)
      real(kind=wp)       :: dtime, coords(nDim,nNode), fieldVar(npredf)

      integer             :: jprops(1:njprops)
      character(len=8)    :: analysis

      real(kind=wp)       :: detF, C(3,3), Cinv(3,3), detC, B(3,3)
      real(kind=wp)       :: Binv(3,3), detB, straintensorEuler(3,3)
      real(kind=wp)       :: trC, lam_c, lam_r, beta_c, dBeta_c
      real(kind=wp)       :: stressTensorPK2(3,3)
      real(kind=wp)       :: stressTensorCauchy(3,3), Cmat(3,3,3,3)
      real(kind=wp)       :: VoigtMat(nSymm,nSymm)
      real(kind=wp)       :: strainVoigtEuler(nSymm,1)
      real(kind=wp)       :: stressVoigtPK1(nUnsymmm,1)
      real(kind=wp)       :: stressVoigtPK2(nSymm,1)
      real(kind=wp)       :: stressVoigtCauchy(nSymm,1)
      real(kind=wp)       :: Gshear, kappa, lam_L

      integer             :: nInt, nlocalSdv
      integer             :: i, j, k, l
      type(logger)        :: msg

      ! initialize matrial stiffness tensors
      Cmat   = zero
      Dmat   = zero

      ! assign material properties to variables
      Gshear= props(1)        ! Shear modulus
      kappa = props(2)        ! Bulk modulus
      lam_L = props(3) 				! Locking stretch

      nInt   = jprops(1)
      nlocalSdv = NSVARS/nInt

      if (lam_L .le. zero) then
        call msg%ferror(flag=error, src='umatArrudaBoyce',
     &       msg='Incorrect material parameter (lam_L).', ra=lam_l)
        call xit
      end if

      ! perform all the constitutitve relations in 3D
      call detMat(F,detF)

      if (detF .le. zero) then
        call msg%ferror(flag=error, src='umatNeoHookean',
     &        msg='Issue with volume change (detF)', 
     &        ivec=[jelem, intpt], ra= detF)
        call xit
      end if

      B = matmul(F,transpose(F))
      C = matmul(transpose(F),F)

      call inverse(B,Binv)
      call inverse(C,Cinv)

      call traceMat(C,trC)

      ! calculate Euler-Almansi strain tensor
      straintensorEuler = half*(ID3-Binv)

      lam_c = sqrt(trC/three)
      lam_r = lam_c/lam_L
      beta_c = InvLangevin(lam_r)
      dBeta_c = DInvLangevin(lam_r)

      ! form material tangent, C_ijkl
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              Cmat(i,j,k,l) = Cmat(i,j,k,l) + Gshear/(nine*lam_c**two)
     &            * (dBeta_c- lam_r*beta_c) * ID3(i,j)*ID3(k,l)
     &            + kappa*Cinv(i,j)*Cinv(k,l)
     &            + ( Gshear/three*lam_r*beta_c-kappa*log(detF) )
     &            * ( Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k) )
            end do
          end do
        end do
      end do

      stressTensorCauchy = (1/detF)*( (Gshear/three)*lam_r*beta_c*B -
     &      (Gshear*lam_L/three - kappa*log(detF))*ID3 )
      stressTensorPK2 = Gshear/three*lam_r*beta_c*ID3 -
     &      (Gshear*lam_L/three - kappa*log(detF))*Cinv


      ! transforms the stiffness tensor 3x3x3x3 to a 6x6 matrix
      call symtangent2matrix(Cmat,VoigtMat)

      ! transform the stress tensor (3x3) to Voigt vector form (6x1)
      call symtensor2vector(straintensorEuler, strainVoigtEuler)
      call symtensor2vector(stressTensorCauchy, stressVoigtCauchy)
      call symtensor2vector(stressTensorCauchy, stressVoigtPK2)

   !!!!!!!!!!!!!!!!! END OF CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!!


      ! reshape the Voigt matrix and tensor based on analysis
      if (analysis .eq. 'PE')  then
        Dmat(1:ndim,1:ndim)     = VoigtMat(1:ndim,1:ndim)
        Dmat(1:ndim,nStress)    = VoigtMat(1:ndim,nSymm)
        Dmat(nStress,1:ndim)    = VoigtMat(nSymm,1:ndim)
        Dmat(nStress,nStress)   = VoigtMat(nSymm,nSymm)

        strainEuler(1:ndim,1)   = strainVoigtEuler(1:ndim,1)
        strainEuler(nStress,1)  = strainVoigtEuler(nSymm,1)

        stressPK2(1:ndim,1)     = stressVoigtPK2(1:ndim,1)
        stressPK2(nStress,1)    = stressVoigtPK2(nSymm,1)

        stressCauchy(1:ndim,1)  = stressVoigtCauchy(1:ndim,1)
        stressCauchy(nStress,1) = stressVoigtCauchy(nSymm,1)

      else if (analysis .eq. '3D') then
        Dmat = VoigtMat
        strainEuler = strainVoigtEuler
        stressPK2 = stressVoigtPK2
        stressCauchy = stressVoigtCauchy

      else
        call msg%ferror(flag=error, src='umatArrudaBoyce',
     &            msg='Wrong analysis.', ch=analysis)
        call xit
      end if

        ! save the variables to be post-processed in globalPostVars
        globalPostVars(jelem,intpt,1:nStress) = stressCauchy(1:nStress,1)
        globalPostVars(jelem,intpt,nStress+1:2*nStress) =
     & 								  strainEuler(1:nStress,1)


      contains

      ! approximation of inverse Langevin function
      ! reference: Bergstorm (PhD thesis, MIT, 1999)
      function InvLangevin(x)

      implicit none

      real(kind=wp),intent(in)  :: x
      real(kind=wp)             :: InvLangevin

      if (abs(x) .lt. 0.84136_wp) then
        InvLangevin = 1.31446_wp*tan(1.58986_wp*x) + 0.91209_wp*x
      else if ((abs(x) .ge. 0.84136_wp) .and. (abs(x) .lt. one)) then
        InvLangevin = one/(sign(one,x)-x)
      else
        call msg%ferror(flag=error, src='umatArrudaBoyce:InvLangevin',
     &                  msg='Unbound argument.', ra = x)
        call xit
      end if

      return
      end function InvLangevin


      ! derivative of inverse Langevin function
      function DInvLangevin(x)

      implicit none

      real(kind=wp), intent(in)   :: x
      real(kind=wp)               :: DInvLangevin, sec

      if (abs(x) .lt. 0.84136_wp) then
        DInvLangevin = 2.0898073756_wp*(tan(1.58986_wp*x))**two
     &                + 3.0018973756_wp
      else if ((abs(x) .ge. 0.84136_wp) .and. (abs(x) .lt. one)) then
        DInvLangevin = one/( (sign(one,x)-x)**two )
      else
       call msg%ferror(flag=error, src='umatArrudaBoyce:InvLangevin',
     &                  msg='Unbound argument.', ra = x)
        call xit
      end if

      return

      end function DInvLangevin

      end subroutine umatArrudaBoyce

! **********************************************************************
! **********************************************************************

       subroutine UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     & NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     & JMAC,JMATYP,MATLAYO,LACCFLA)
       ! this subroutine is used to transfer postVars from the UEL
       ! onto the dummy mesh for viewing. Note that an offset of
       ! elemOffset is used between the real mesh and the dummy mesh.

      use global_parameters
      use post_processing

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      ! the dimensions of the variables FLGRAY, ARRAY and JARRAY
      ! must be set equal to or greater than 15.
      ! explicityly define the type for uvar to avoid issues

      real(kind=wp)   :: uvar

      uvar(1:nuvarm)  = globalPostVars(noel-elemOffset,npt,1:nuvarm)

      end subroutine UVARM

! **********************************************************************