! **********************************************************************
! ********** ABAQUS/STANDARD USER ELEMENT (UEL) SUBROUTINE *************
! **********************************************************************
!                   BIBEKANANDA DATTA (C) FEBRUARY 2024
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
!     uvar(1:ntens)                 Cauchy stress tensor components 
!     uvar(ntens+1:2*ntens)         Euler strain tensor components
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
!     NSEND MODULEVARS                        Total # element state variables
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


! **********************************************************************
!                        GLOBAL PARAMETERS MODULE
!
!     defines real(kind=wp) real global parameters to be used in 
!     various functions and subroutines within and outside of UEL
! **********************************************************************

      module global_parameters

      integer, parameter :: wp = kind(0.0d0)

      real(kind=wp) :: zero, one, two, three, four, five,
     &        six, eight, nine, half, third, fourth, 
     &        fifth, sixth, eighth, eps, pi

      parameter( zero = 0.0_wp, one =  1.0_wp, two = 2.0_wp,
     &        three = 3.0_wp, four = 4.0_wp, five = 5.0_wp,
     &        six = 6.0_wp, eight = 8.0_wp, nine = 9.0_wp,
     &        half= 0.5_wp, third = 1.0_wp/3.0_wp, 
     &        fourth = 0.25_wp, fifth = 1.0_wp/5.0_wp, 
     &        sixth = 1.0_wp/6.0_wp, eighth = 0.125_wp,
     &        eps = 1.d-09,
     &        pi = 3.14159265358979323846264338327950_wp)


      ! no of symmetric and unSymmmetric tensor components
      integer, parameter  :: nSymm = 6, nUnsymmm = 9

      ! identity matrices in 2- and 3-dimensions
      real(kind=wp) :: ID2(2,2), ID3(3,3)
      parameter(  ID2 = reshape([ one, zero,
     &                            zero, one ], shape(ID2) ),
     &            ID3 = reshape([ one, zero, zero,
     &                            zero, one, zero,
     &                            zero, zero, one ], shape(ID3) ) )

      !! debug file unit
      integer, parameter :: dFile = 15


      !! post processing variables
      ! no of total user elements and offset for overlaying elements
      integer, parameter:: numElem = 50000, elemOffset = 100000
      ! numElem should be same (recommended) or higher than total
      ! number of elements in the overlaying element set
      ! elemOffset should match with the number used in input file

      real(kind=wp), allocatable:: globalPostVars(:,:,:)

      end module global_parameters

! **********************************************************************
! ****************** ABAQUS USER ELEMENT SUBROUTINE ********************
! **********************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)

      use global_parameters

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     & DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      real(kind=wp):: RHS, AMATRX, SVARS, ENERGY, PNEWDT,
     &      PROPS, COORDS, DUall, Uall, Vel, Accn, TIME, DTIME,
     &      PARAMS, ADLMAG, PREDEF, DDLMAG, PERIOD

      integer:: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE,
     &      KSTEP, KINC, JELEM, NDLOAD, JDLTYP, NPREDF, LFLAGS,
     &      MLVARX, MDLOAD, JPROPS, NJPROPS

      integer:: nInt, matID, nPostVars
      integer:: nDim, ndi, nshr, ntens, uDOF, uDOFEL

      integer:: lenJobName,lenOutDir
      character(len=256)  :: outDir, fileName
      character(len=64)   :: jobName
      character(len=8)    :: analysis, abq_procedure
      logical :: nlgeom

      ! open a debug file for the current job
      call getJobName(jobName,lenJobName)
      jobName = trim(jobName)
      call getOutDir(outDir,lenOutDir)
      fileName = outDir(1:lenOutDir)//'\aaMSG_'//
     &     jobName(1:lenJobName)//'.dat'
      open(unit=dFile,file=fileName,status='unknown')

      ! assign parameter specific to analysis and element types
      if ((JTYPE.ge.1).and.(JTYPE.le.4)) then
        analysis = '3D'         ! three-dimensional analysis
        nDim = 3
        ndi = 3
        nshr = 3
        ntens = 6
      elseif ((JTYPE.ge.5).and.(JTYPE.le.8)) then
        analysis = 'PE'         ! plane strain analysis
        nDim = 2
        ndi = 2
        nshr = 1
        ntens = 3
      else
        write(dFile,*) 'ERROR > UEL: Element type is available: ', JTYPE
        write(*,*) 'ERROR > UEL: Element type is available: ', JTYPE
        call xit
      endif

      ! change the LFLAGS criteria as needed (check abaqus UEL manual)
      if((lflags(1).eq.1).or.(lflags(1).eq.2)) then
        ABQ_PROCEDURE = 'STATIC'
      else
        write(dFile,*) 'ERROR > UEL: Incorrect procedure: ', lflags(1)
        write(*,*) 'ERROR > UEL: Incorrect procedure: ', lflags(1)
        call xit
      endif

      ! check if the procedure is linear or nonlinear
      if (lflags(2).eq.0) then
        nlgeom = .false.
      elseif (lflags(2).eq.1) then
        nlgeom = .true.
      endif

      ! check to see if it's a general step or a linear purturbation step
      if(lflags(4).eq.1) then
        write(dFile,*) 'ERROR > UEL: ',
     &     'The step should be a GENERAL step: ', lflags(4)
        write(*,*) 'ERROR > UEL: ',
     &      'The step should be a GENERAL step: ', lflags(4)
        call xit
      endif

      ! for mixed or coupled problem, add other DOF counts as needed
      uDOF = nDim             ! displacement degrees of freedom of a node
      uDOFEL = nNode*uDOF     ! total displacement degrees of freedom in element

      nInt = jprops(1)
      matID = jprops(2)
      nPostVars = jprops(3)

      ! array containing variables for post-processing
      if (.not. allocated(globalPostVars)) then
        allocate(globalPostVars(numElem,nInt,nPostVars))

        ! print necessary information to the debug file (one time)
        write(dFile,*) '---------------------------------------'
        write(dFile,*) '------- ABAQUS FINITE STRAIN UEL ------'
        write(dFile,*) '---------------------------------------'
        write(dFile,*) 'Abaqus Job: ', jobName
        write(dFile,*) '---------------------------------------'
        write(dFile,*) '------- PROCEDURE = ', ABQ_PROCEDURE
        write(dFile,*) '------- ANALYSIS TYPE   = ', analysis
        write(dFile,*) '---------- NLGEOM = ', NLGEOM
        write(dFile,*) '------- MODEL DIMENSION = ', nDim
        write(dFile,*) '------- ELEMENT NODES   = ', NNODE
        write(dFile,*) '---------------------------------------'
        write(dFile,*) '-------- INTEGRATION SCHEME -----------'
        write(dFile,*) '----------- NINT  = ', nInt
        write(dFile,*) '---------------------------------------'
        write(dFile,*) '---------- POST-PROCESSING ------------'
        write(dFile,*) '--- NO OF ELEMENTS            = ', numElem
        write(dFile,*) '--- OVERLAY ELEMENT OFFSET    = ', ElemOffset
        write(dFile,*) '--- NO OF VARIABLES AT INT PT = ', nPostVars
        write(dFile,*) '---------------------------------------'

        ! print necessary information to the screen now (one time)
        write(*,*) '---------------------------------------'
        write(*,*) '------- ABAQUS FINITE STRAIN UEL ------'
        write(*,*) '---------------------------------------'
        write(*,*) 'Abaqus Job: ', jobName
        write(*,*) '---------------------------------------'
        write(*,*) '------- PROCEDURE = ', ABQ_PROCEDURE
        write(*,*) '------- ANALYSIS TYPE   = ', analysis
        write(*,*) '---------- NLGEOM = ', NLGEOM
        write(*,*) '------- MODEL DIMENSION = ', nDim
        write(*,*) '------- ELEMENT NODES   = ', NNODE
        write(*,*) '---------------------------------------'
        write(*,*) '-------- INTEGRATION SCHEME -----------'
        write(*,*) '----------- NINT  = ', nInt
        write(*,*) '---------------------------------------'
        write(*,*) '---------- POST-PROCESSING ------------'
        write(*,*) '--- NO OF ELEMENTS            = ', numElem
        write(*,*) '--- OVERLAY ELEMENT OFFSET    = ', ElemOffset
        write(*,*) '--- NO OF VARIABLES AT INT PT = ', nPostVars
        write(*,*) '---------------------------------------'

      endif

       ! call your UEL subroutine
       call uelNLMECH(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & ANALYSIS,NDIM,NDI,NSHR,NTENS,NINT,UDOF,UDOFEL)


      RETURN 

      END SUBROUTINE UEL

! **********************************************************************

      subroutine uelNLMECH(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & ANALYSIS,NDIM,NDI,NSHR,NTENS,NINT,UDOF,UDOFEL)

      use global_parameters

      implicit none

     !!!!!!!!!!!!! VARIABLE DECLARATION AND INITIALIZATION !!!!!!!!!!!!!
     
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     & DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      real(kind=wp):: RHS, AMATRX, SVARS, ENERGY, PNEWDT,
     &      PROPS, COORDS, DUall, Uall, Vel, Accn, TIME, DTIME,
     &      PARAMS, ADLMAG, PREDEF, DDLMAG, PERIOD

      integer:: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE,
     &      KSTEP, KINC, JELEM, NDLOAD, JDLTYP, NPREDF, LFLAGS,
     &      MLVARX, MDLOAD, JPROPS, NJPROPS

      character(len=8):: analysis
      integer:: nDim, ndi, nshr, ntens, nInt, uDOF, uDOFEL, nlocalSdv
      logical:: nlgeom

      real(kind=wp):: uNode(nDim,nNode), duNode(nDim,nNode), F(3,3),
     &      ID(nDim,nDim), w(nInt), xi(nInt,nDim),
     &      Nxi(nNode,1), dNdxi(nNode,nDim), dxdxi(nDim,nDim),
     &      dxidx(nDim,nDim), dNdx(nNode,nDim), detJ,
     &      Na(nDim,nDim), Nmat(nDim,uDOFEl),
     &      Ba(ntens,nDim), Bmat(ntens,uDOFEl),
     &      Ga(nDim*nDim,nDim), Gmat(nDim*nDim,uDOFEl),
     &      fieldNode(npredf,nNode), dfieldNode(npredf,nNode),
     &      fieldVar(npredf), dfieldVar(npredf),
     &      stressTensorPK2(ntens,ntens), SIGMA_S(nDim**2,nDim**2),
     &      SIGMA_F(nDim*nNode,nDim*nNode),
     &      kuu(uDOFEl,uDOFEl), Ru(uDOFEl,1)

      real(kind=wp):: stranLagrange(ntens,1), stranEuler(ntens,1),
     &      stressCauchy(ntens,1), stressPK1(nDim*nDim,1),
     &      stressPK2(ntens,1), Dmat(ntens,ntens)

      integer:: i, j, intPt, matID, info

      ! initialize the matrices and vectors
      F  = zero
      Na = zero
      Ba = zero
      Ga = zero
      Nmat = zero
      Bmat  = zero
      Gmat = zero
      Dmat = zero
      SIGMA_F = zero
      SIGMA_S = zero
      kuu   = zero
      Ru    = zero

      ENERGY = zero
      AMATRX(1:NDOFEL,1:NDOFEL) = zero
      RHS(1:MLVARX,1) = zero

      matID = jprops(2)

     !!!!!!!!!!! END VARIABLE DECLARATION AND INITIALIZATION !!!!!!!!!!!

      ! reshape the displacement vectors into matrix forms
      uNode  = reshape(UAll,(/nDim,nNode/))
      duNode = reshape(DUAll(:,1),(/nDim,nNode/))

      ! if applicable gather the prescribed field variables in a vector
      ! such as temperature (as shown below in commented line - not tested)
      ! fieldNode(1,1:nNode) = predef(1,1,1:nNode)
      ! dfieldNode(1,1:nNode) = predef(2,1,1:nNode)


     !!!!!!!!!!!!!!!!! ELEMENT RELATED OPERATIONS !!!!!!!!!!!!!!!!!!!!!!

      ! obtain gauss quadrature points and weights
      if (nDim.eq.2) then
        ID = ID2
        call gaussQuadrtr2(nNode,nInt,w,xi)
      elseif (nDim.eq.3) then
        ID = ID3
        call gaussQuadrtr3(nNode,nInt,w,xi)
      else
        write(dFile,*) 'ERROR > uelNLMech: Incorrect dimension: ', nDim
        write(*,*) 'ERROR > uelNLMech: Incorrect dimension: ', nDim
        call xit
      endif


      ! loop through all the integration points (main/ external loop)
      do intPt = 1, nInt

        ! obtain the shape functions and their gradient for the element
        if (nDim.eq.2) then
          call interpFunc2(nNode,nInt,intPt,xi,Nxi,dNdxi)
        elseif (nDim.eq.3) then
          call interpFunc3(nNode,nInt,intPt,xi,Nxi,dNdxi)
        else
          write(dFile,*)'ERROR > uelNLMech: Incorrect dimension: ',nDim
          write(*,*) 'ERROR > uelNLMech: Incorrect dimension: ',nDim
          call xit
        endif

        ! calculate element jacobian and global gradient
        dXdxi   = matmul(coords,dNdxi)      ! calculate dXdxi

        if (nDim.eq.2) then
          call inverseMat2(dXdxi,dxidX,detJ,info)
        elseif(nDim.eq.3) then
          call inverseMat3(dXdxi,dxidX,detJ,info)
        endif

        if (info .eq. 0) then
          write(dFile,*) 'WARNING > uelNLMech: Element jacobian: ', 
     &                    jElem, intPt
          write(*,*) 'WARNING > uelNLMech: Element jacobian: ', 
     &                    jElem, intPt
        endif

        dNdX    = matmul(dNdxi,dxidX)       ! calculate dNdX

        ! loop over all the nodes (internal loop)
        do i=1,nNode

          ! form the nodal-level matrices: [Na], [Ga], [Ba]
          do j = 1, nDim
            Na(j,j) = Nxi(i,1)
            Ga(nDim*(j-1)+1:nDim*j,1:nDim) = dNdX(i,j)*ID
          enddo

          ! form [Ba] matrix: plane stress/ plane strain case
          if (analysis.eq.'PE') then
            Ba(1,1) = dNdx(i,1)
            Ba(2,2) = dNdx(i,2)
            Ba(3,1:nDim) = [dNdx(i,2), dNdx(i,1)]

          ! form [Ba] matrix: 3D case
          elseif (analysis.eq.'3D') then
            Ba(1,1) = dNdx(i,1)
            Ba(2,2) = dNdx(i,2)
            Ba(3,3) = dNdx(i,3)
            Ba(4,1:nDim) = [zero, dNdx(i,3), dNdx(i,2)]
            Ba(5,1:nDim) = [dNdx(i,3), zero,  dNdx(i,1)]
            Ba(6,1:nDim) = [dNdx(i,2), dNdx(i,1), zero]

          else
            write(dFile,*)'ERROR > uelNLMech: Wrong analysis: ',analysis
            write(*,*) 'ERROR > uelNLMech: Wrong analysis: ', analysis
            call xit
          endif

          ! form the [N], [B], and [G] matrix
          Nmat(1:nDim,nDim*(i-1)+1:nDim*i) = Na(1:nDim,1:nDim)
          Bmat(1:nTens,nDim*(i-1)+1:nDim*i) = Ba(1:nTens,1:nDim)
          Gmat(1:nDim**2,nDim*(i-1)+1:nDim*i) = Ga(1:nDim**2,1:nDim)
        enddo                             ! end of nodal point loop

      !!!!!!!!!!!!!! COMPLETE ELEMENT RELATED OPERATIONS !!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!! CONSTITUTIVE MODEL !!!!!!!!!!!!!!!!!!!!!!!

        ! calculate deformation gradient and deformation tensors
        F(1:nDim,1:nDim) = ID + matmul(uNode,dNdX)
        
        if (analysis .eq. 'PE') then
          F(3,3) = one
        endif

        ! interpolate the field variable at the integration point
        ! (as shown below - not tested)
    !     fieldVar = dot_product(reshape(Nxi,(/nNode/)),
    !  &                        reshape(fieldNode,(/nNode/)))
    !     dfieldVar = dot_product(reshape(Nxi,(/nNode/)),
    !  &                        reshape(dfieldNode,(/nNode/)))


        ! call material point subroutine (UMAT) for specific material
        if (matID .eq. 1) then
          call umatNeoHookean(stressCauchy,stressPK1,stressPK2,
     &          Dmat,F,svars,nsvars,stranLagrange,stranEuler,time,
     &          dtime,fieldVar,npredf,nDim,ndi,nshr,ntens,jelem,intpt,
     &          coords,nnode,kstep,kinc,props,nprops,jprops,njprops,
     &          analysis,nlocalSdv)

        elseif (matID .eq. 2) then
          call umatArrudaBoyce(stressCauchy,stressPK1,stressPK2,
     &          Dmat,F,svars,nsvars,stranLagrange,stranEuler,time,
     &          dtime,fieldVar,npredf,nDim,ndi,nshr,ntens,jelem,intPt,
     &          coords,nnode,kstep,kinc,props,nprops,jprops,njprops,
     &          analysis,nlocalSdv)
        endif
        ! can add more constitutive models using elseif construct here

      !!!!!!!!!!!!!!!!!!!! END CONSTITUTIVE MODEL !!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!! TANGENT MATRIX AND RESIDUAL VECTOR !!!!!!!!!!!!!!!!!
        if (nDim .eq. 2) then
          call vector2symtensor2(stressPK2,stressTensorPK2)
        elseif (nDim .eq. 3) then
          call vector2symtensor3(stressPK2,stressTensorPK2)
        endif

        ! form the [SIGMA_S] matrix for geometric stiffness
        do i = 1, nDim
          do j = 1, nDim
              SIGMA_S(nDim*(i-1)+1:nDim*i,nDim*(j-1)+1:nDim*j) =
     &                    stressTensorPK2(i,j)*ID
          enddo
        enddo

        ! form [SIGMA_F] matrix for material stiffness
        do i = 1, nNode
          do j = 1, nNode
              if (i.eq.j) then                    ! banded diagonal
                  SIGMA_F(nDim*(i-1)+1:nDim*i,nDim*(j-1)+1:nDim*j)
     &                              = F(1:nDim,1:nDim)
              endif
          enddo
        enddo

        ! form the stiffness matrix and residual vector
        kuu = kuu + w(intpt)*detJ*
     &      ( matmul( transpose(matmul(Bmat,SIGMA_F)),
     &      matmul (Dmat,matmul(Bmat,SIGMA_F)) ) +
     &      matmul( transpose(Gmat), matmul(SIGMA_S,Gmat)) )

        Ru  = Ru - w(intpt)*detJ*
     &            matmul( transpose(matmul(Bmat,SIGMA_F)), stressPK2 )

      !!!!!!!!!!!!!! TANGENT MATRIX AND RESIDUAL VECTOR !!!!!!!!!!!!!!!!

      enddo                         ! end of integration point loop

      ! body force and surface load can be added using dummy elements

      ! assign the element stiffness matrix to abaqus-defined variable
      AMATRX(1:NDOFEL,1:NDOFEL) = kuu(1:uDOFEl,1:uDOFEl)
      RHS(1:MLVARX,1) = Ru(1:uDOFEl,1)

      return

      end subroutine uelNLMECH

! **********************************************************************

      subroutine umatNeoHookean(stressCauchy,stressPK1,stressPK2,
     &      Dmat,F,svars,nsvars,stranLagrange,stranEuler,time,
     &      dtime,fieldVar,npredf,nDim,ndi,nshr,ntens,jelem,intpt,
     &      coords,nnode,kstep,kinc,props,nprops,jprops,njprops,
     &      analysis,nlocalSdv)

      use global_parameters

      implicit none

      integer:: nsvars, npredf, nDim, ndi, nshr, ntens,
     &    jelem, intpt, nNode, kstep, kinc, nprops, njprops

      real(kind=wp):: stressCauchy(ntens,1), stressPK1(nDim*nDim,1),
     &    stressPK2(ntens,1), Dmat(ntens,ntens), F(3,3),
     &    stranLagrange(ntens,1), stranEuler(ntens,1), props(1:nprops),
     &    svars(1:nsvars), coords(nDim,nNode), time(2), dtime,
     &    fieldVar(npredf)

      integer :: jprops(1:njprops)
      character(len=8)  :: analysis

      real(kind=wp) :: detF, C(3,3), Cinv(3,3), detC, B(3,3),
     &    Binv(3,3), detB, strantensorEuler(3,3), stressTensorPK2(3,3),
     &    stressTensorCauchy(3,3), Cmat(3,3,3,3), VoigtMat(nSymm,nSymm),
     &    stranVoigtEuler(nSymm,1), stressVoigtPK1(nUnsymmm,1),
     &    stressVoigtPK2(nSymm,1), stressVoigtCauchy(nSymm,1),
     &    Gshear, kappa, lam_L

      integer:: nInt, nlocalSdv

      ! loop counters
      integer:: i, j, k, l, info

      ! initialize matrial stiffness tensors
      Cmat   = zero
      Dmat   = zero

      ! assign material properties to variables
      Gshear= props(1)        ! Shear modulus
      kappa = props(2)        ! Bulk modulus
      lam_L = props(3)        ! locking stretch

      nInt   = jprops(1)
      nlocalSdv = NSVARS/nInt

      ! locking stretch should be infinity (0 as input) for NH model
      if (lam_L .ne. zero) then
        write(dFile,*)  'ERROR > umatNeoHookean: ', 
     &                  'Incorrect material parameter (lam_L)', lam_L
        write(*,*)  'ERROR > umatNeoHookean: ', 
     &              'Incorrect material parameter (lam_L)', lam_L
        call xit
      endif

      ! perform all the constitutitve relations in 3D
      call detMat3(F,detF)

      if (detF .le. zero) then
        write(dFile,*) 'WARRNING > umatNeoHookean: ', 
     &    'Issue with volume change(detF): ', jelem, intPt, detF

        write(*,*) 'WARRNING > umatNeoHookean: ', 
     &    'Issue with volume change (detF): ', jelem, intPt, detF
      endif


      B = matmul(F,transpose(F))
      C = matmul(transpose(F),F)

      call inverseMat3(B,Binv,detB,info)
      call inverseMat3(C,Cinv,detC,info)

      ! calculate Euler-Almansi strain tensor
      strantensorEuler = half*(ID3-Binv)

      ! calculate material tangent, C_ijkl
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              Cmat(i,j,k,l) = Cmat(i,j,k,l)
     &            + kappa*Cinv(i,j)*Cinv(k,l)
     &            + (Gshear-kappa*log(detF)) * ( Cinv(i,k)*Cinv(j,l)
     &            + Cinv(i,l)*Cinv(j,k) )
            enddo
          enddo
        enddo
      enddo

      ! calculate stress tensors
      stressTensorCauchy = (one/detF)*(Gshear*(B-ID3) +
     & 											kappa*log(detF)*ID3)
      stressTensorPK2 = Gshear*(ID3-Cinv) + kappa*log(detF)*Cinv


      ! transforms the stiffness tensor 3x3x3x3 to a 6x6 matrix
      call symtangent2matrix(Cmat,VoigtMat)

      ! transform the stress tensor (3x3) to Voigt vector form (6x1)
      call symtensor2vector3(strantensorEuler,stranVoigtEuler)
      call symtensor2vector3(stressTensorCauchy,stressVoigtCauchy)
      call symtensor2vector3(stressTensorCauchy,stressVoigtPK2)

     !!!!!!!!!!!!!!! END OF CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!!


      ! reshape the Voigt matrix and tensor based on analysis
      if (analysis .eq. 'PE')  then
        Dmat(1:ndi,1:ndi) = VoigtMat(1:ndi,1:ndi)
        Dmat(1:ndi,ntens) = VoigtMat(1:ndi,nSymm)
        Dmat(ntens,1:ndi) = VoigtMat(nSymm,1:ndi)
        Dmat(ntens,ntens) = VoigtMat(nSymm,nSymm)

        stranEuler(1:ndi,1) = stranVoigtEuler(1:ndi,1)
        stranEuler(ntens,1) = stranVoigtEuler(nSymm,1)

        stressPK2(1:ndi,1) = stressVoigtPK2(1:ndi,1)
        stressPK2(ntens,1) = stressVoigtPK2(nSymm,1)

        stressCauchy(1:ndi,1) = stressVoigtCauchy(1:ndi,1)
        stressCauchy(ntens,1) = stressVoigtCauchy(nSymm,1)

      elseif (analysis .eq. '3D') then
        Dmat = VoigtMat
        stranEuler = stranVoigtEuler
        stressPK2 = stressVoigtPK2
        stressCauchy = stressVoigtCauchy
      endif
    

      ! save the variables to be post-processed in globalPostVars
      globalPostVars(jelem,intpt,1:ntens) = stressCauchy(1:ntens,1)
      globalPostVars(jelem,intpt,ntens+1:2*ntens) =
     &                     stranEuler(1:ntens,1)

      return

      end subroutine umatNeoHookean

! **********************************************************************

      subroutine umatArrudaBoyce(stressCauchy,stressPK1,stressPK2,
     &      Dmat,F,svars,nsvars,stranLagrange,stranEuler,time,
     &      dtime,fieldVar,npredf,nDim,ndi,nshr,ntens,jelem,intpt,
     &      coords,nnode,kstep,kinc,props,nprops,jprops,njprops,
     &      analysis,nlocalSdv)

      use global_parameters

      implicit none

      integer:: nsvars, npredf, nDim, ndi, nshr, ntens,
     &    jelem, intpt, nNode, kstep, kinc, nprops, njprops

      real(kind=wp):: stressCauchy(ntens,1), stressPK1(nDim*nDim,1),
     &    stressPK2(ntens,1), Dmat(ntens,ntens), F(3,3),
     &    stranLagrange(ntens,1), stranEuler(ntens,1), props(1:nprops),
     &    svars(1:nsvars), coords(nDim,nNode), time(2), dtime,
     &    fieldVar(npredf)

      integer:: jprops(1:njprops)
      character(len=8):: analysis

      real(kind=wp):: detF, C(3,3), Cinv(3,3), detC, B(3,3), 
     & 		Binv(3,3), detB, trC, lam_c, lam_r, beta_c, dBeta_c,
     &    strantensorEuler(3,3), stressTensorPK2(3,3),
     &    stressTensorCauchy(3,3), Cmat(3,3,3,3), VoigtMat(nSymm,nSymm),
     &    stranVoigtEuler(nSymm,1), stressVoigtPK1(nUnsymmm,1),
     &    stressVoigtPK2(nSymm,1), stressVoigtCauchy(nSymm,1),
     &    Gshear, kappa, lam_L

      integer:: nInt, nlocalSdv
      ! loop counters
      integer:: i, j, k, l, info

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
         write(dFile,*)  'ERROR > umatArrudaBoyce: ', 
     &                   'Incorrect material parameter (lam_L)', lam_L
        write(*,*)  'ERROR > umatArrudaBoyce ', 
     &              'Incorrect material parameter (lam_L)', lam_L
        call xit
      endif

      ! perform all the constitutitve relations in 3D
      call detMat3(F,detF)
      if (detF .le. zero) then
        write(dFile,*) 'WARRNING > umatArrudaBoyce: ', 
     &    'Issue with volume change(detF): ', jelem, intPt, detF

        write(*,*) 'WARRNING > umatArrudaBoyce: ', 
     &    'Issue with volume change (detF): ', jelem, intPt, detF
      endif

      B = matmul(F,transpose(F))
      C = matmul(transpose(F),F)

      call inverseMat3(B,Binv,detB,info)
      call inverseMat3(C,Cinv,detC,info)

      call traceMat(C,trC,size(C,1))

      ! calculate Euler-Almansi strain tensor
      strantensorEuler = half*(ID3-Binv)

      lam_c = sqrt(trC/3.)
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
     &            + ( Gshear/three*lam_r*beta_c-kappa*log(detF) )*
     &              ( Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k) )
            enddo
          enddo
        enddo
      enddo

      stressTensorCauchy = (one/detF)*( (Gshear/three)*lam_r*beta_c*B -
     &      (Gshear*lam_L/three - kappa*log(detF))*ID3 )
      stressTensorPK2 = Gshear/three*lam_r*beta_c*ID3 -
     &      (Gshear*lam_L/three - kappa*log(detF))*Cinv


      ! transforms the stiffness tensor 3x3x3x3 to a 6x6 matrix
      call symtangent2matrix(Cmat,VoigtMat)

      ! transform the stress tensor (3x3) to Voigt form (6x1)
      call symtensor2vector3(strantensorEuler,stranVoigtEuler)
      call symtensor2vector3(stressTensorCauchy,stressVoigtCauchy)
      call symtensor2vector3(stressTensorCauchy,stressVoigtPK2)

   !!!!!!!!!!!!!!!!! END OF CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!!


      ! reshape the Voigt matrix and tensor based on analysis
      if (analysis .eq. 'PE')  then
        Dmat(1:ndi,1:ndi) = VoigtMat(1:ndi,1:ndi)
        Dmat(1:ndi,ntens) = VoigtMat(1:ndi,nSymm)
        Dmat(ntens,1:ndi) = VoigtMat(nSymm,1:ndi)
        Dmat(ntens,ntens) = VoigtMat(nSymm,nSymm)

        stranEuler(1:ndi,1) = stranVoigtEuler(1:ndi,1)
        stranEuler(ntens,1) = stranVoigtEuler(nSymm,1)

        stressPK2(1:ndi,1) = stressVoigtPK2(1:ndi,1)
        stressPK2(ntens,1) = stressVoigtPK2(nSymm,1)

        stressCauchy(1:ndi,1) = stressVoigtCauchy(1:ndi,1)
        stressCauchy(ntens,1) = stressVoigtCauchy(nSymm,1)

      elseif (analysis .eq. '3D') then
        Dmat = VoigtMat
        stranEuler = stranVoigtEuler
        stressPK2 = stressVoigtPK2
        stressCauchy = stressVoigtCauchy
      endif

        ! save the variables to be post-processed in globalPostVars
        globalPostVars(jelem,intpt,1:ntens) = stressCauchy(1:ntens,1)
        globalPostVars(jelem,intpt,ntens+1:2*ntens) =
     & 				  stranEuler(1:ntens,1)

      return

      contains


      ! approximation of inverse Langevin function
      ! reference: Bergstorm (PhD thesis, MIT, 1999)
      function InvLangevin(x)

      implicit none
      
      real(kind=wp), intent(in) :: x
      real(kind=wp):: InvLangevin


      if (abs(x) .lt. 0.84136_wp) then
        InvLangevin = 1.31446_wp*tan(1.58986_wp*x) + 0.91209_wp*x
      elseif ((abs(x) .ge. 0.84136_wp) .and. (abs(x) .lt. one)) then
        InvLangevin = one/(sign(one,x)-x)
      else
        write(dFile,*) 'ERROR > InvLangevin: Unbound argument: ', x
        write(*,*) 'ERROR > InvLangevin: Unbound argument: ', x
        call xit
      endif

      return

      end function InvLangevin


      ! derivative of inverse Langevin function
      function DInvLangevin(x)

      implicit none

      real(kind=wp), intent(in) :: x
      real(kind=wp) :: DInvLangevin, sec

      if (abs(x) .lt. 0.84136_wp) then
        DInvLangevin =  2.0898073756_wp*(tan(1.58986_wp*x))**two 
     &                  + 3.0018973756_wp
      elseif ((abs(x) .ge. 0.84136_wp) .and. (abs(x) .lt. one)) then
        DInvLangevin = one/((sign(one,x)-x)**two)
      else
        write(dFile,*) 'ERROR > InvLangevin: Unbound argument: ', x
        write(*,*) 'ERROR > InvLangevin: Unbound argument: ', x
        call xit
      endif

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

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      ! the dimensions of the variables FLGRAY, ARRAY and JARRAY
      ! must be set equal to or greater than 15.  

      ! explicityly define the type for uvar to avoid issues
      real(kind=wp):: uvar

      uvar(1:nuvarm) = globalPostVars(noel-elemOffset,npt,1:nuvarm)

      return

      end subroutine UVARM

! **********************************************************************
! **********************************************************************




! **********************************************************************
! **************** INTERPOLATION FUNCTION SUBROUTINES ******************
! **********************************************************************
!  available elements:    (a) 1D bar/truss element (2 and 3 nodes)
!                         (b) 2D tri elements (3 and 6 nodes) 
!                         (c) 2D quad elements (4 and 8 nodes)
!                         (d) 3D tet elements (4 and 10 nodes)
!                         (e) 3D hex elements (8 and 20 nodes)
! **********************************************************************

      subroutine interpFunc1(nNode,nInt,intPt,xi_int,Nxi,dNdxi)
      ! this subroutine evaluates 1D shape functions
      ! available 1D element is: 2 and 3 noded bar

      ! Nxi(i)          = shape function of node i at the intpt
      ! dNdxi(i,j)      = derivative wrt j direction of shape fn of node i

      use global_parameters

      implicit none
      
      integer:: nNode, nInt, intpt
      
      real(kind=wp) :: xi, xi_int(nInt,1), Nxi(nNode), dNdxi(nNode,1)

      xi    = xi_int(intpt,1)

      if (nNode .eq. 2) then      ! 2 node linear bar element
        ! shape functions
        Nxi(1) = half*(one - xi)
        Nxi(2) = half*(one + xi)

        ! the first derivatives of the shape functions dN/dxi (2x1)
        dNdxi(1,1)  = -half
        dNdxi(2,1)  = half

      elseif (nNode .eq. 3) then  ! 3 node quadratic bar element
        ! shape functions
        Nxi(1)  = -half*xi*(one - xi)
        Nxi(2)  = one-xi**two
        Nxi(3)  = half*xi*(one + xi)

        ! the first derivatives of the shape functions dN/dxi (3x1)
        dNdxi(1,1)  = -half+xi
        dNdxi(2,1)  = -two*xi
        dNdxi(3,1)  = half+xi

      else
        write(dFile,*) 'ERROR > interpFunc1: Element unavailable.'
        write(*,*) 'ERROR > interpFunc1: Element unavailable.'
        return
      endif

      return

      end subroutine interpFunc1

! **********************************************************************

      subroutine interpFunc2(nNode,nInt,intPt,xi_int,Nxi,dNdxi)
      ! this subroutine calculates shape function of 2D elements
      ! available 2D elements are: 3 and 6 node tri, 4 and 8 node quad
        
      ! Nxi(i)          = shape function of node i at the intpt
      ! dNdxi(i,j)      = derivative wrt j direction of shape fn of node i

      use global_parameters

      implicit none

      integer:: nNode, nInt, intpt

      real(kind=wp):: xi_int(nInt,2), Nxi(nNode), dNdxi(nNode,2)
      real(kind=wp):: xi, eta, lam

      ! location in the master element
      xi    = xi_int(intpt,1)
      eta   = xi_int(intpt,2)

      !              A eta (=xi_2)
      !              |
      !              3
      !              |\
      !              | \
      !              |  \
      !              |   \
      !              |    \
      !              |     \
      !              |      \
      !              |       \
      !              1--------2--> xi (=xi_1)

      if (nNode.eq.3) then        ! 3-noded tri3 linear element
        ! shape functions
        Nxi(1) = xi
        Nxi(2) = eta
        Nxi(3) = one - xi - eta

        ! the first derivatives of the shape functions dN/dxi (3x2)
        dNdxi(1, 1) = one
        dNdxi(1, 2) = zero
        dNdxi(2, 1) = zero
        dNdxi(2, 2) = one
        dNdxi(3, 1) = -one
        dNdxi(3, 2) = -one

      !              A eta (=xi_2)
      !              |
      !              3
      !              |\
      !              | \
      !              |  \
      !              |   5
      !              6    \
      !              |     \
      !              |      \
      !              |       \
      !              1---4----2--> xi (=xi_1)

      elseif (nNode.eq.6) then    ! 6-noded quadratic tri6 element
        ! shape functions
        lam = one - xi - eta
        Nxi(1) = lam*(two*lam - one)
        Nxi(2) = xi*(two*xi - one)
        Nxi(3) = eta*(two*eta - one)
        Nxi(4) = four*xi*lam
        Nxi(5) = four*xi*eta
        Nxi(6) = four*eta*lam

        ! the first derivatives of the shape functions dN/dxi (6x2)
        dNdxi(1,1) = one - four*lam
        dNdxi(1,2) = one - four*lam
        dNdxi(2,1) = four*xi - one
        dNdxi(2,2) = zero
        dNdxi(3,1) = zero
        dNdxi(3,2) = four*eta - one
        dNdxi(4,1) = four*(lam - xi)
        dNdxi(4,2) = -four*xi
        dNdxi(5,1) = four*eta
        dNdxi(5,2) = four*xi
        dNdxi(6,1) = -four*eta
        dNdxi(6,2) = four*(lam - eta)

      !                          eta
      !   4-----------3          |
      !   |           |          |
      !   |           |          |
      !   |           |          O--------- xi
      !   |           |
      !   |           |        origin at center
      !   1-----------2

      elseif (nNode.eq.4) then    ! 4-noded bilinear quad4 element
        ! shape functions
        Nxi(1) = fourth*(one - xi)*(one - eta)
        Nxi(2) = fourth*(one + xi)*(one - eta)
        Nxi(3) = fourth*(one + xi)*(one + eta)
        Nxi(4) = fourth*(one - xi)*(one + eta)

        ! the first derivatives of the shape functions dN/dxi (4x2)
        dNdxi(1,1) = -fourth*(one - eta)
        dNdxi(1,2) = -fourth*(one - xi)
        dNdxi(2,1) = fourth*(one - eta)
        dNdxi(2,2) = -fourth*(one + xi)
        dNdxi(3,1) = fourth*(one + eta)
        dNdxi(3,2) = fourth*(one + xi)
        dNdxi(4,1) = -fourth*(one + eta)
        dNdxi(4,2) = fourth*(one - xi)

      !         A eta
      !         |
      !         |
      !   4-----7-----3
      !   |     |     |
      !   |     |     |
      !   8     ------6---> xi
      !   |           |
      !   |           |
      !   1-----5-----2

      elseif (nNode.eq.8) then    ! 8-noded serendipity quad8 element
        ! shape functions
        Nxi(1) = -fourth*(one - xi)*(one - eta)*(one + xi + eta)
        Nxi(2) = -fourth*(one + xi)*(one - eta)*(one - xi + eta)
        Nxi(3) = -fourth*(one + xi)*(one + eta)*(one - xi - eta)
        Nxi(4) = -fourth*(one - xi)*(one + eta)*(one + xi - eta)
        Nxi(5) = half*(one - xi)*(one + xi)*(one - eta)
        Nxi(6) = half*(one - eta)*(one + eta)*(one + xi)
        Nxi(7) = half*(one - xi)*(one + xi)*(one + eta)
        Nxi(8) = half*(one - eta)*(one + eta)*(one - xi)

        ! the first derivatives of the shape functions dN/dxi (8x2)
        dNdxi(1,1) = fourth*(one - eta)*(two*xi + eta)
        dNdxi(1,2) = fourth*(one - xi)*(two*eta + xi)
        dNdxi(2,1) = fourth*(one - eta)*(two*xi - eta)
        dNdxi(2,2) = fourth*(one + xi)*(two*eta - xi)
        dNdxi(3,1) = fourth*(one + eta)*(two*xi + eta)
        dNdxi(3,2) = fourth*(one + xi)*(two*eta + xi)
        dNdxi(4,1) = fourth*(one + eta)*(two*xi - eta)
        dNdxi(4,2) = fourth*(one - xi)*(two*eta - xi)
        dNdxi(5,1) = -xi*(one - eta)
        dNdxi(5,2) = -half*(one - xi)*(one + xi)
        dNdxi(6,1) = half*(one - eta)*(one + eta)
        dNdxi(6,2) = -eta*(one + xi)
        dNdxi(7,1) = -xi*(one + eta)
        dNdxi(7,2) = half*(one - xi)*(one + xi)
        dNdxi(8,1) = -half*(one - eta)*(one + eta)
        dNdxi(8,2) = -eta*(one - xi)

      else
        write(dFile,*) 'ERROR > interpFunc2: Element unavailable.'
        write(*,*) 'ERROR > interpFunc2: Element unavailable.'
        return
      endif

      return

      end subroutine interpFunc2

! **********************************************************************

      subroutine interpFunc3(nNode,nInt,intPt,xi_int,Nxi,dNdxi)
      ! this subroutine calculates shape function of 3D elements
      ! available 3D elements are: 4 and 10 node tet, 8 and 20 node hex.

      ! Nxi(i)          = shape function of node i at the intpt
      ! dNdxi(i,j)      = derivative wrt j direction of shape fn of node i

      use global_parameters

      implicit none

      integer:: nNode, nInt, intpt

      real(kind=wp):: xi_int(nInt,3), Nxi(nNode), dNdxi(nNode,3)
      real(kind=wp):: xi, eta, zeta, lam

      ! Nxi(i)          = shape function of node i at the intpt.
      ! dNdxi(i,j)      = derivative wrt j direction of shape fn of node i

      ! location in the master element
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)

      Nxi = zero
      dNdxi = zero

      if (nNode.eq.4) then      ! 4-noded linear tet4 element
        ! shape functions
        Nxi(1) = one-xi-eta-zeta
        Nxi(2) = xi
        Nxi(3) = eta
        Nxi(4) = zeta
        ! the first derivatives of the shape functions dN/dxi (4x3)
        dNdxi(1,1) = -one
        dNdxi(1,2) = -one
        dNdxi(1,3) = -one
        dNdxi(2,1) = one
        dNdxi(3,2) = one
        dNdxi(4,3) = one

      elseif (nNode.eq.10) then  ! 10-noded quadratic tet10 element
    
        ! shape functions
        lam = one-xi-eta-zeta

        Nxi(1) = (two*lam-one)*lam
        Nxi(2) = (two*xi-one)*xi
        Nxi(3) = (two*eta-one)*eta
        Nxi(4) = (two*zeta-one)*zeta
        Nxi(5) = four*xi*lam
        Nxi(6) = four*xi*eta
        Nxi(7) = four*eta*lam
        Nxi(8) = four*zeta*lam
        Nxi(9) = four*zeta*xi
        Nxi(10) = four*eta*zeta
        
        ! first derivative of shape functions dN/dxi (10x3)
        dNdxi(1,1) = -(four*lam-one)
        dNdxi(1,2) = -(four*lam-one)
        dNdxi(1,3) = -(four*lam-one)
        dNdxi(2,1) = (four*xi-one)
        dNdxi(3,2) = (four*eta-one)
        dNdxi(4,3) = (four*zeta-one)
        dNdxi(5,1) = four*(lam-xi)
        dNdxi(5,2) = -four*xi
        dNdxi(5,3) = -four*xi
        dNdxi(6,1) = four*eta
        dNdxi(6,2) = four*xi
        dNdxi(7,1) = -four*eta
        dNdxi(7,2) = four*(lam-eta)
        dNdxi(7,3) = -four*eta
        dNdxi(8,1) = -four*zeta
        dNdxi(8,2) = -four*zeta
        dNdxi(8,3) = four*(lam-zeta)
        dNdxi(9,1) = four*zeta
        dNdxi(9,3) = four*xi
        dNdxi(10,2) = four*zeta
        dNdxi(10,3) = four*eta

      !      8-----------7
      !     /|          /|
      !    / |         / |      zeta
      !   5-----------6  |       |   eta
      !   |  |        |  |       |   /
      !   |  |        |  |       |  /
      !   |  4--------|--3       | /
      !   | /         | /        |/
      !   |/          |/         O--------- xi
      !   1-----------2        origin at cube center

      elseif(nNode.eq.8) then   ! 8-noded trilinear hex8 element

        ! shape functions
        Nxi(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
        Nxi(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
        Nxi(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
        Nxi(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
        Nxi(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
        Nxi(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
        Nxi(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
        Nxi(8) = eighth*(one - xi)*(one + eta)*(one + zeta)

        ! the first derivatives: dN/dxi (8x3)
        dNdxi(1,1) = -eighth*(one - eta)*(one - zeta)
        dNdxi(1,2) = -eighth*(one - xi)*(one - zeta)
        dNdxi(1,3) = -eighth*(one - xi)*(one - eta)
        dNdxi(2,1) = eighth*(one - eta)*(one - zeta)
        dNdxi(2,2) = -eighth*(one + xi)*(one - zeta)
        dNdxi(2,3) = -eighth*(one + xi)*(one - eta)
        dNdxi(3,1) = eighth*(one + eta)*(one - zeta)
        dNdxi(3,2) = eighth*(one + xi)*(one - zeta)
        dNdxi(3,3) = -eighth*(one + xi)*(one + eta)
        dNdxi(4,1) = -eighth*(one + eta)*(one - zeta)
        dNdxi(4,2) = eighth*(one - xi)*(one - zeta)
        dNdxi(4,3) = -eighth*(one - xi)*(one + eta)
        dNdxi(5,1) = -eighth*(one - eta)*(one + zeta)
        dNdxi(5,2) = -eighth*(one - xi)*(one + zeta)
        dNdxi(5,3) = eighth*(one - xi)*(one - eta)
        dNdxi(6,1) = eighth*(one - eta)*(one + zeta)
        dNdxi(6,2) = -eighth*(one + xi)*(one + zeta)
        dNdxi(6,3) = eighth*(one + xi)*(one - eta)
        dNdxi(7,1) = eighth*(one + eta)*(one + zeta)
        dNdxi(7,2) = eighth*(one + xi)*(one + zeta)
        dNdxi(7,3) = eighth*(one + xi)*(one + eta)
        dNdxi(8,1) = -eighth*(one + eta)*(one + zeta)
        dNdxi(8,2) = eighth*(one - xi)*(one + zeta)
        dNdxi(8,3) = eighth*(one - xi)*(one + eta)

      !
      !       8-----15----- 7
      !      /|            /|
      !    16 |          14 |
      !    /  20        /   |     zeta
      !   5-----13-----6   19      |     eta
      !   |   |        |    |      |    /
      !   |   |        |    |      |   /
      !   17  4-----11-|----3      |  /
      !   |  /         18  /       | /
      !   | 12         | 10        |/
      !   |/           |/          O--------- xi
      !   1-----9------2        origin at cube center
      !
      ! mid-side nodes are not properly illustrated

      elseif (nNode.eq.20) then   ! 20-noded serendipity hex20 element
        ! shape functions
        Nxi(1) = (one-xi)*(one-eta)*(one-zeta)*(-xi-eta-zeta-two)/eight
        Nxi(2) = (one+xi)*(one-eta)*(one-zeta)*(xi-eta-zeta-two)/eight
        Nxi(3) = (one+xi)*(one+eta)*(one-zeta)*(xi+eta-zeta-two)/eight
        Nxi(4) = (one-xi)*(one+eta)*(one-zeta)*(-xi+eta-zeta-two)/eight
        Nxi(5) = (one-xi)*(one-eta)*(one+zeta)*(-xi-eta+zeta-two)/eight
        Nxi(6) = (one+xi)*(one-eta)*(one+zeta)*(xi-eta+zeta-two)/eight
        Nxi(7) = (one+xi)*(one+eta)*(one+zeta)*(xi+eta+zeta-two)/eight
        Nxi(8) = (one-xi)*(one+eta)*(one+zeta)*(-xi+eta+zeta-two)/eight
        Nxi(9)  = (one-xi**two)*(one-eta)*(one-zeta)/four
        Nxi(10) = (one+xi)*(one-eta**two)*(one-zeta)/four
        Nxi(11) = (one-xi**two)*(one+eta)*(one-zeta)/four
        Nxi(12) = (one-xi)*(one-eta**two)*(one-zeta)/four
        Nxi(13) = (one-xi**two)*(one-eta)*(one+zeta)/four
        Nxi(14) = (one+xi)*(one-eta**two)*(one+zeta)/four
        Nxi(15) = (one-xi**two)*(one+eta)*(one+zeta)/four
        Nxi(16) = (one-xi)*(one-eta**two)*(one+zeta)/four
        Nxi(17) = (one-xi)*(one-eta)*(one-zeta**two)/four
        Nxi(18) = (one+xi)*(one-eta)*(one-zeta**two)/four
        Nxi(19) = (one+xi)*(one+eta)*(one-zeta**two)/four
        Nxi(20) = (one-xi)*(one+eta)*(one-zeta**two)/four

        ! the first derivatives: dN/dxi (8x3)
        dNdxi(1,1) = (-(one-eta)*(one-zeta)*(-xi-eta-zeta-two)-
     &                (one-xi)*(one-eta)*(one-zeta))/eight
        dNdxi(1,2) = (-(one-xi)*(one-zeta)*(-xi-eta-zeta-two)-
     &                (one-xi)*(one-eta)*(one-zeta))/eight
        dNdxi(1,3) = (-(one-xi)*(one-eta)*(-xi-eta-zeta-two)-
     &                (one-xi)*(one-eta)*(one-zeta))/eight
        dNdxi(2,1) = ((one-eta)*(one-zeta)*(xi-eta-zeta-two)+
     &                (one+xi)*(one-eta)*(one-zeta))/eight
        dNdxi(2,2) = (-(one+xi)*(one-zeta)*(xi-eta-zeta-two)-
     &                (one+xi)*(one-eta)*(one-zeta))/eight
        dNdxi(2,3) = (-(one+xi)*(one-eta)*(xi-eta-zeta-two)-
     &                (one+xi)*(one-eta)*(one-zeta))/eight
        dNdxi(3,1) = ((one+eta)*(one-zeta)*(xi+eta-zeta-two)+
     &                (one+xi)*(one+eta)*(one-zeta))/eight
        dNdxi(3,2) = ((one+xi)*(one-zeta)*(xi+eta-zeta-two)+
     &                (one+xi)*(one+eta)*(one-zeta))/eight
        dNdxi(3,3) = (-(one+xi)*(one+eta)*(xi+eta-zeta-two)-
     &                (one+xi)*(one+eta)*(one-zeta))/eight
        dNdxi(4,1) = (-(one+eta)*(one-zeta)*(-xi+eta-zeta-two)-
     &                (one-xi)*(one+eta)*(one-zeta))/eight
        dNdxi(4,2) = ((one-xi)*(one-zeta)*(-xi+eta-zeta-two)+
     &                (one-xi)*(one+eta)*(one-zeta))/eight
        dNdxi(4,3) = (-(one-xi)*(one+eta)*(-xi+eta-zeta-two)-
     &              (one-xi)*(one+eta)*(one-zeta))/eight
        dNdxi(5,1) = (-(one-eta)*(one+zeta)*(-xi-eta+zeta-two)-
     &                (one-xi)*(one-eta)*(one+zeta))/eight
        dNdxi(5,2) = (-(one-xi)*(one+zeta)*(-xi-eta+zeta-two)-
     &                (one-xi)*(one-eta)*(one+zeta))/eight
        dNdxi(5,3) = ((one-xi)*(one-eta)*(-xi-eta+zeta-two)+
     &                (one-xi)*(one-eta)*(one+zeta))/eight
        dNdxi(6,1) = ((one-eta)*(one+zeta)*(xi-eta+zeta-two)+
     &                (one+xi)*(one-eta)*(one+zeta))/eight
        dNdxi(6,2) = (-(one+xi)*(one+zeta)*(xi-eta+zeta-two)-
     &                (one+xi)*(one-eta)*(one+zeta))/eight
        dNdxi(6,3) = ((one+xi)*(one-eta)*(xi-eta+zeta-two)+
     &                (one+xi)*(one-eta)*(one+zeta))/eight
        dNdxi(7,1) = ((one+eta)*(one+zeta)*(xi+eta+zeta-two)+
     &                (one+xi)*(one+eta)*(one+zeta))/eight
        dNdxi(7,2) = ((one+xi)*(one+zeta)*(xi+eta+zeta-two)+
     &                (one+xi)*(one+eta)*(one+zeta))/eight
        dNdxi(7,3) = ((one+xi)*(one+eta)*(xi+eta+zeta-two)+
     &                (one+xi)*(one+eta)*(one+zeta))/eight
        dNdxi(8,1) = (-(one+eta)*(one+zeta)*(-xi+eta+zeta-two)-
     &                (one-xi)*(one+eta)*(one+zeta))/eight
        dNdxi(8,2) = ((one-xi)*(one+zeta)*(-xi+eta+zeta-two)+
     &                (one-xi)*(one+eta)*(one+zeta))/eight
        dNdxi(8,3) = ((one-xi)*(one+eta)*(-xi+eta+zeta-two)+
     &                (one-xi)*(one+eta)*(one+zeta))/eight
        dNdxi(9,1)  = -two*xi*(one-eta)*(one-zeta)/four
        dNdxi(9,2)  = -(one-xi**two)*(one-zeta)/four
        dNdxi(9,3)  = -(one-xi**two)*(one-eta)/four
        dNdxi(10,1)  = (one-eta**two)*(one-zeta)/four
        dNdxi(10,2)  = -two*eta*(one+xi)*(one-zeta)/four
        dNdxi(10,3)  = -(one-eta**two)*(one+xi)/four
        dNdxi(11,1)  = -two*xi*(one+eta)*(one-zeta)/four
        dNdxi(11,2)  = (one-xi**two)*(one-zeta)/four
        dNdxi(11,3)  = -(one-xi**two)*(one+eta)/four
        dNdxi(12,1)  = -(one-eta**two)*(one-zeta)/four
        dNdxi(12,2)  = -two*eta*(one-xi)*(one-zeta)/four
        dNdxi(12,3)  = -(one-eta**two)*(one-xi)/four
        dNdxi(13,1)  = -two*xi*(one-eta)*(one+zeta)/four
        dNdxi(13,2)  = -(one-xi**two)*(one+zeta)/four
        dNdxi(13,3)  = (one-xi**two)*(one-eta)/four
        dNdxi(14,1)  = (one-eta**two)*(one+zeta)/four
        dNdxi(14,2)  = -two*eta*(one+xi)*(one+zeta)/four
        dNdxi(14,3)  = (one-eta**two)*(one+xi)/four
        dNdxi(15,1)  = -two*xi*(one+eta)*(one+zeta)/four
        dNdxi(15,2)  =  (one-xi**two)*(one+zeta)/four
        dNdxi(15,3)  = (one-xi**two)*(one+eta)/four
        dNdxi(16,1)  = -(one-eta**two)*(one+zeta)/four
        dNdxi(16,2)  = -two*eta*(one-xi)*(one+zeta)/four
        dNdxi(16,3)  = (one-eta**two)*(one-xi)/four
        dNdxi(17,1) = -(one-eta)*(one-zeta**two)/four
        dNdxi(17,2) = -(one-xi)*(one-zeta**two)/four
        dNdxi(17,3) = -zeta*(one-xi)*(one-eta)/two
        dNdxi(18,1) = (one-eta)*(one-zeta**two)/four
        dNdxi(18,2) = -(one+xi)*(one-zeta**two)/four
        dNdxi(18,3) = -zeta*(one+xi)*(one-eta)/two
        dNdxi(19,1) = (one+eta)*(one-zeta**two)/four
        dNdxi(19,2) = (one+xi)*(one-zeta**two)/four
        dNdxi(19,3) = -zeta*(one+xi)*(one+eta)/two
        dNdxi(20,1) = -(one+eta)*(one-zeta**two)/four
        dNdxi(20,2) = (one-xi)*(one-zeta**two)/four
        dNdxi(20,3) = -zeta*(one-xi)*(one+eta)/two
      else
        write(dFile,*) 'ERROR > interpFunc3: Element unavailable.'
        write(*,*) 'ERROR > interpFunc3: Element unavailable.'
        return
      endif

      return
      end subroutine interpFunc3

! **********************************************************************

      subroutine faceNodes(nDim,nNode,face,list,nFaceNodes)
      ! this subroutine returns the list of nodes on an
      ! element face for standard 2D and 3D Lagrangian elements
      ! this subroutine is useful for applying traction-type BC

      implicit none

      integer, intent (in) :: nDim, nNode, face
      integer, intent (out):: list(*)
      integer, intent (out):: nFaceNodes

      integer :: list3(3), list4(4)

      if (nDim.eq.2) then
        list3(1:3) = [2,3,1]
        list4(1:4) = [2,3,4,1]

        if (nNode.eq.3) then
          nFaceNodes = 2
          list(1) = face
          list(2) = list3(face)
        else if (nNode.eq.4) then
          nFaceNodes = 2
          list(1) = face
          list(2) = list4(face)
        else if (nNode.eq.6) then
          nFaceNodes = 3
          list(1) = face
          list(2) = list3(face)
          list(3) = face+3
        else if (nNode.eq.8) then
          nFaceNodes = 3
          list(1) = face
          list(2) = list4(face)
          list(3) = face+4
        endif

      else if (nDim.eq.3) then

        if (nNode.eq.4) then
          nFaceNodes = 3
          if (face.eq.1) list(1:3) = [1,2,3]
          if (face.eq.2) list(1:3) = [1,4,2]
          if (face.eq.3) list(1:3) = [2,4,3]
          if (face.eq.4) list(1:3) = [3,4,1]
        else if (nNode .eq.6) then
          nFaceNodes = 3
          if (face.eq.1) list(1:3) = [1,2,3]
          if (face.eq.2) list(1:3) = [6,5,4]
          if (face.eq.3) list(1:4) = [1,2,5,4]
          if (face.eq.4) list(1:4) = [2,3,6,5]
          if (face.eq.5) list(1:4) = [4,6,3,1]
          if (face>2) nFaceNodes = 4
        else if (nNode.eq.10) then
          nFaceNodes = 6
          if (face.eq.1) list(1:6) = [1,2,3,5,6,7]
          if (face.eq.2) list(1:6) = [1,4,2,8,9,5]
          if (face.eq.3) list(1:6) = [2,4,3,9,10,6]
          if (face.eq.4) list(1:6) = [3,4,1,10,8,7]
        else if (nNode.eq.8) then
          nFaceNodes = 4
          if (face.eq.1) list(1:4) = [1,2,3,4]
          if (face.eq.2) list(1:4) = [5,8,7,6]
          if (face.eq.3) list(1:4) = [1,5,6,2]
          if (face.eq.4) list(1:4) = [2,6,7,3]
          if (face.eq.5) list(1:4) = [3,7,8,4]
          if (face.eq.6) list(1:4) = [4,8,5,1]
        else  if (nNode.eq.20) then
          nFaceNodes = 8
          if (face.eq.1) list(1:8) = [1,2,3,4,9,10,11,12]
          if (face.eq.2) list(1:8) = [5,8,7,6,16,15,14,13]
          if (face.eq.3) list(1:8) = [1,5,6,2,17,13,18,9]
          if (face.eq.4) list(1:8) = [2,6,7,3,18,14,19,10]
          if (face.eq.5) list(1:8) = [3,7,8,4,19,15,6,11]
          if (face.eq.6) list(1:8) = [4,8,5,1,20,16,17,12]
        endif
      endif

      return

      end subroutine faceNodes

! **********************************************************************

! **********************************************************************
! ***************** GAUSSIAN QUADRATURE SUBROUTINES ********************
! **********************************************************************
! integration schemes:  (a) reduced and full integration: bar elements
!                       (b) full integration: tri and tet elements
!                       (c) reduced and full integration: quad and hex elements
! **********************************************************************
      subroutine gaussQuadrtr1(nNode,nInt,w,xi)

      use global_parameters

      implicit none

      integer:: nNode, nInt
      real(kind=wp):: w(nInt), xi(nInt,1)

      w = zero
      xi = zero

      if (nNode.eq.2) then
        if (nInt.eq.1) then       ! full integration for bar2
          w = zero
          xi = two
        else
          write(dFile,*) 'ERROR > gaussQuadrtr1: wrong Gauss points.'
          write(*,*) 'ERROR > gaussQuadrtr1: wrong Gauss points.'
          return
        endif           ! end int for 2-noded bar

      elseif (nNode.eq.3) then
        if (nInt .eq. 2) then     ! reduced integration for bar3
          w(1:2)  = one
          xi(1,1) = -sqrt(third)
          xi(2,1) = sqrt(third)
        elseif (nInt.eq.3) then   ! full integration for bar3
          w(1) = five/nine
          w(2) = eight/nine
          w(3) = five/nine
          xi(1,1) = -sqrt(three/five)
          xi(2,1) = zero
          xi(3,1) = sqrt(three/five)
        else
          write(dFile,*) 'ERROR > gaussQuadrtr1: wrong Gauss points.'
          write(*,*) 'ERROR > gaussQuadrtr1: wrong Gauss points.'
          return
        endif           ! end int for 3-noded bar

      else
        write(dFile,*) 'ERROR > gaussQuadrtr1: Element unavailable.'
        write(*,*) 'ERROR > gaussQuadrtr1: Element unavailable.'
        return
      endif           

      return

      end subroutine gaussQuadrtr1

! **********************************************************************

      subroutine gaussQuadrtr2(nNode,nInt,w,xi)
      ! this subroutines returns the weights and
      ! gauss point coordinate of 2D Lagrangian elements
      ! currently supports: nInt = 1 (tri3) and 3 (tri6)
      !                     nInt = 1, 4 (quad4) and 4, 9 (quad8)

      use global_parameters
      implicit none

      integer:: nNode, nInt
      real(kind=wp):: x1D(4), w1D(4)
      real(kind=wp):: w(nInt), xi(nInt,2)

      w  = zero
      xi = zero

      if (nNode.eq.3) then      ! plane tri3 elements (full integration)
        if (nInt.eq.1) then
          w(1) = half
          xi(1,1) = third
          xi(2,1) = third
        else
          write(dFile,*) 'ERROR > gaussQuadrtr2: Wrong Gauss points.'
          write(*,*) 'ERROR > gaussQuadrtr2: Wrong Gauss points.'
          return
        endif

      elseif (nNode.eq.6) then  ! plane tri6 elements (full integration)
        if (nInt.eq.3) then
          w(1:3) = sixth

          xi(1,1) = half
          xi(1,2) = half
          xi(2,1) = zero
          xi(2,2) = half
          xi(3,1) = half
          xi(3,2) = zero
        else
          write(dFile,*) 'ERROR > gaussQuadrtr2: Wrong Gauss points.'
          write(*,*) 'ERROR > gaussQuadrtr2: Wrong Gauss points.'
          return
        endif
        
      elseif((nNode.eq.4)) then ! plane quad4 element

        if (nInt.eq.1) then     ! reduced integration for quad4
          w = four

          xi(1,1) = zero
          xi(1,2) = zero

        elseif (nInt.eq.4) then ! full integration for quad4

          w(1:4) = one

          x1D(1) = sqrt(third)
          xi(1,1) = -x1D(1)
          xi(1,2) = -x1D(1)
          xi(2,1) = x1D(1)
          xi(2,2) = -x1D(1)
          xi(3,1) = -x1D(1)
          xi(3,2) = x1D(1)
          xi(4,1) = x1D(1)
          xi(4,2) = x1D(1)
        else
          write(dFile,*) 'ERROR > gaussQuadrtr2: Wrong Gauss points.'
          write(*,*) 'ERROR > gaussQuadrtr2: Wrong Gauss points.'
          return
        endif

      elseif (nNode.eq.8) then  ! plane quad8 element
        
        if (nInt.eq.4) then     ! reduced integration for quad8

          w(1:4) = one

          x1D(1) = sqrt(third)
          xi(1,1) = -x1D(1)
          xi(1,2) = -x1D(1)
          xi(2,1) = x1D(1)
          xi(2,2) = -x1D(1)
          xi(3,1) = -x1D(1)
          xi(3,2) = x1D(1)
          xi(4,1) = x1D(1)
          xi(4,2) = x1D(1)

        elseif(nInt.eq.9) then  ! full integration for quad8
          w1D(1) = (five/nine)*(five/nine)
          w1D(2) = (five/nine)*(eight/nine)
          w1D(3) = (eight/nine)*(eight/nine)

          w(1) = w1D(1)
          w(2) = w1D(2)
          w(3) = w1D(1)
          w(4) = w1D(2)
          w(5) = w1D(3)
          w(6) = w1D(2)
          w(7) = w1D(1)
          w(8) = w1D(2)
          w(9) = w1D(1)

          x1D(1) = sqrt(three/five)
          xi(1,1) = -x1D(1)
          xi(1,2) = -x1D(1)
          xi(2,1) = zero
          xi(2,2) = -x1D(1)
          xi(3,1) = x1D(1)
          xi(3,2) = -x1D(1)
          xi(4,1) = -x1D(1)
          xi(4,2) = zero
          xi(5,1) = zero
          xi(5,2) = zero
          xi(6,1) = x1D(1)
          xi(6,2) = zero
          xi(7,1) = -x1D(1)
          xi(7,2) = x1D(1)
          xi(8,1) = zero
          xi(8,2) = x1D(1)
          xi(9,1) = x1D(1)
          xi(9,2) = x1D(1)

        else
          write(dFile,*) 'ERROR > gaussQuadrtr2: Wrong Gauss points.'
          write(*,*) 'ERROR > gaussQuadrtr2: Wrong Gauss points.'
          return
        endif

      else
        write(dFile,*) 'ERROR > gaussQuadrtr2: Element unavailable.'
        write(*,*) 'ERROR > gaussQuadrtr2: Element unavailable.'
        return
      endif


      return
      end subroutine gaussQuadrtr2

! **********************************************************************

      subroutine gaussQuadrtr3(nNode,nInt,w,xi)
      ! this subroutines returns the weights and
      ! gauss point coordinate of 3D Lagrangian elements
      ! currently supports: nInt = 1 (tet4) and 4 (tet10)
      !                     nInt = 1, 8 (hex8) and 8, 27 (hex20)

      use global_parameters

      implicit none

      integer:: nNode, nInt
      integer:: i, j, k, n
      real(kind=wp):: x1D(4), w1D(4)
      real(kind=wp):: w(nInt), xi(nInt,3)

      w  = zero
      xi = zero
      
      if(nNode.eq.4) then       ! 3D tet4 element (full integration)
        if(nInt.eq.1) then
          w(1) = sixth
          xi(1,1:3) = fourth

        else
          write(dFile,*) 'ERROR > gaussQuadrtr3: Wrong Gauss points.'
          write(*,*) 'ERROR > gaussQuadrtr3: Wrong Gauss points.'
          return
        endif

      elseif(nNode.eq.10) then  ! 3D tet10 element (full integration)

        if (nInt.eq.4) then
          w(1:4) = one/24.0_wp

          x1D(1) = 0.58541020_wp
          x1D(2) = 0.13819660_wp

          xi(1,1) = x1D(1)
          xi(2,1) = x1D(2)
          xi(3,1) = x1D(2)
          xi(1,2) = x1D(2)
          xi(2,2) = x1D(1)
          xi(3,2) = x1D(2)
          xi(1,3) = x1D(2)
          xi(2,3) = x1D(2)
          xi(3,3) = x1D(1)
          xi(4,1) = x1D(2)
          xi(4,2) = x1D(2)
          xi(4,3) = x1D(2)

        else
          write(dFile,*) 'ERROR > gaussQuadrtr3: Wrong Gauss points.'
          write(*,*) 'ERROR > gaussQuadrtr3: Wrong Gauss points.'
          return
        endif

      elseif(nNode.eq.8) then   ! 3D hex8 element

        if(nInt.eq.1) then      ! reduced integration for hex8
          w(1) = eight
          xi(1,1:3) = zero

        elseif(nInt.eq.8) then  ! full-integration for hex8
          w(1:8) = one

          x1D(1) = -sqrt(third)
          x1D(2) = sqrt(third)
         
          do k = 1,2
            do j = 1,2
              do i = 1,2
                n = 4*(k-1) + 2*(j-1) + i
                xi(n,1) = x1D(i)
                xi(n,2) = x1D(j)
                xi(n,3) = x1D(k)
              end do
            end do
          end do

        else
          write(dFile,*) 'ERROR > gaussQuadrtr3: Wrong Gauss points.'
          write(*,*) 'ERROR > gaussQuadrtr3: Wrong Gauss points.'
          return
        endif

      elseif(nNode.eq.20) then  ! 3D hex20 element
        
        if (nInt.eq.8) then     ! reduced integration for hex20
          w(1:8) = one

          x1D(1) = -sqrt(third)
          x1D(2) = sqrt(third)
         
          do k = 1,2
            do j = 1,2
              do i = 1,2
                n = 4*(k-1) + 2*(j-1) + i
                xi(n,1) = x1D(i)
                xi(n,2) = x1D(j)
                xi(n,3) = x1D(k)
              end do
            end do
          end do
        
        elseif(nInt.eq.27) then ! full integration for hex20
          w1D(1) = five/nine
          w1D(2) = eight/nine
          w1D(3) = w1D(1)

          x1D(1) = -sqrt(0.6_wp)
          x1D(2) = zero
          x1D(3) = sqrt(0.6_wp)
          do k = 1,3
            do j = 1,3
              do i = 1,3
                n = 9*(k-1) + 3*(j-1) + i
                xi(n,1) = x1D(i)
                xi(n,2) = x1D(j)
                xi(n,3) = x1D(k)
                w(n) = w1D(i)*w1D(j)*w1D(k)
              end do
            end do
          end do

        else
          write(dFile,*) 'ERROR > gaussQuadrtr3: Wrong Gauss points.'
          write(*,*) 'ERROR > gaussQuadrtr3: Wrong Gauss points.'
          return
        endif

      else
        write(dFile,*) 'ERROR > gaussQuadrtr3: Element unavailable.'
        write(*,*) 'ERROR > gaussQuadrtr3: Element unavailable.'
        return
      endif

      return
      end subroutine gaussQuadrtr3

! **********************************************************************


! **********************************************************************
*********************** MATRIX ALGEBRA SECTION *************************
! **********************************************************************
!
!     SUBROUTINE to create identity matrix of any dimention
!     SUBROUTINE to calculate determinant of matrix
!     SUBROUTINE to calculate direct inverse of 2x2 and 3x3 matrix
!     SUBROUTINE to map symmetric tensor to a vector
!     SUBROUTINE to map 4th order tensor to 2D Voigt matrix
! **********************************************************************

      subroutine crossProduct(a,b,c)
      ! this subroutine computes the cross product of two 3 dimensional vectors

      use global_parameters, only: wp
      implicit none

      real(kind=wp):: a(3), b(3), c(3)

      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = b(1)*a(3)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)

      end subroutine crossProduct

! **********************************************************************

      subroutine traceMat(A,trA,nDim)
      ! this subroutine calculates the trace of a square matrix

      use global_parameters

      implicit none

      integer:: nDim, i
      real(kind=wp):: A(nDim,nDim), trA

      trA = zero

      do i = 1, nDim
        trA = trA + A(i,i)
      enddo

      return
      end subroutine traceMat

! **********************************************************************

      subroutine detMat2(A,detA)
      ! this subroutine calculates the determinant of a 2x2 or 3x3 matrix [A]

      use global_parameters, only: wp

      implicit none

      real(kind=wp):: A(2,2), detA

      detA = A(1,1)*A(2,2) - A(1,2)*A(2,1)

      return
      end subroutine detMat2

! **********************************************************************

      subroutine detMat3(A,detA)
      ! this subroutine calculates the determinant of a 2x2 or 3x3 matrix [A]

      use global_parameters, only: wp

      implicit none

      real(kind=wp):: A(3,3), detA

      detA = A(1,1)*A(2,2)*A(3,3)
     &     + A(1,2)*A(2,3)*A(3,1)
     &     + A(1,3)*A(2,1)*A(3,2)
     &     - A(3,1)*A(2,2)*A(1,3)
     &     - A(3,2)*A(2,3)*A(1,1)
     &     - A(3,3)*A(2,1)*A(1,2)

      return
      end subroutine detMat3

! **********************************************************************

      subroutine inverseMat2(A,Ainv,detA,info)
      ! this subroutine returns Ainv and detA for a 2D matrix A

      use global_parameters

      implicit none

      integer:: info
      real(kind=wp):: A(2,2),Ainv(2,2), detA, detAinv

      info = 1

      call detMat2(A,detA)

      if (detA .le. zero) then
          write(dFile,*) 'WARNING > inverseMAt2: det of mat= ', detA
          write(*,*) 'WARNING > inverseMAt2: det of mat= ', detA
          info = 0
      end if

      detAinv = one/detA

      Ainv(1,1) =  detAinv*A(2,2)
      Ainv(1,2) = -detAinv*A(1,2)
      Ainv(2,1) = -detAinv*A(2,1)
      Ainv(2,2) =  detAinv*A(1,1)

      return

      end subroutine inverseMat2

! **********************************************************************

      subroutine inverseMat3(A,Ainv,detA,info)
      ! this subroutine returns Ainv and detA for a 3D matrix A

      use global_parameters

      implicit none

      integer:: info
      real(kind=wp):: A(3,3),Ainv(3,3), detA, detAinv

      info = 1

      call detMat3(A,detA)

      if (detA .le. zero) then
        write(dFile,*) 'WARNING > inverseMAt3: det of mat= ', detA
        write(*,*) 'WARNING > inverseMAt3: det of mat= ', detA
        info = 0
      end if

      detAinv = one/detA

      Ainv(1,1) = detAinv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      Ainv(1,2) = detAinv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      Ainv(1,3) = detAinv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      Ainv(2,1) = detAinv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      Ainv(2,2) = detAinv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      Ainv(2,3) = detAinv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      Ainv(3,1) = detAinv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      Ainv(3,2) = detAinv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      Ainv(3,3) = detAinv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))

      return
      end subroutine inverseMat3

! **********************************************************************

      subroutine inverseMat(A,Ainv,n)
      ! this subroutine computes the inverse of an arbitrary
      ! square matrix of size nxn using LU decomposition

      use global_parameters

      implicit none

      integer,intent(in)   :: n

      real(kind=wp),intent(inout) :: A(n,n)
      real(kind=wp),intent(out)   :: Ainv(n,n)

      real(kind=wp):: L(n,n), U(n,n), b(n), d(n), x(n)
      real(kind=wp):: coeff
      integer :: i, j, k

      L = zero
      U = zero
      b = zero

      do k=1, n-1
        do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          A(i,k+1:n) = A(i,k+1:n)-coeff*A(k,k+1:n)
        end do
      end do

      forall (i=1:n)  L(i,i) = one
      forall (j=1:n) U(1:j,j) = A(1:j,j)

      do k=1,n
        b(k)=one
        d(1) = b(1)
        do i=2,n
          d(i)=b(i)
          d(i) = d(i) - dot_product(L(i,1:i-1),d(1:i-1))
        end do

        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
          x(i) = d(i)
          x(i)=x(i)-dot_product(U(i,i+1:n),x(i+1:n))
          x(i) = x(i)/U(i,i)
        end do
        Ainv(1:n,k) = x(1:n)
        b(k)=zero
      end do

      end subroutine inverseMat

! **********************************************************************

      subroutine voigtAugment(vect2D, vect3D, ntens)
      ! this subroutine augemnts a 3x1 (plane) or 4x1 (axisymmetric)
      ! array to a 6x1 Voigt array of 3D dimensional case

      use global_parameters

      implicit none

      integer:: ntens
      real(kind=wp):: vect2D(ntens,1), vect3D(nSymm,1)

      vect3D = zero

      vect3D(1,1) = vect2D(1,1)
      vect3D(2,1) = vect2D(2,1)

      if (ntens.eq.3) then      ! plane strain/stress
        vect3D(6,1) = vect2D(3,1)
      elseif (ntens.eq.4) then  ! axisymmetry
        vect3D(3,1) = vect2D(3,1)
        vect3D(6,1) = vect2D(4,1)
      elseif(ntens.eq.nSymm) then
        vect3D = vect2D
      endif

      return
      end subroutine voigtAugment

! **********************************************************************

      subroutine voigtTruncate(vect3D, vect2D, ntens)
      ! this subroutine truncates a 6x1 Voigt array
      ! to a 3x1 (plane) or 4x1 (axisymmetry) Voigt array

      use global_parameters

      implicit none

      integer:: ntens
      real(kind=wp):: vect3D(nSymm,1), vect2D(ntens,1)

      vect2D = zero

      vect2D(1,1) = vect3D(1,1)
      vect2D(2,1) = vect3D(2,1)

      if (ntens.eq. 3) then
        vect2D(3,1) = vect3D(6,1)
      elseif (ntens.eq.4) then
        vect2D(3,1) = vect3D(3,1)
        vect2D(4,1) = vect3D(6,1)
      endif

      return

      end subroutine voigtTruncate

! **********************************************************************

      subroutine vector2symtensor2(Avect,Atens)
      ! this subroutine transforms a 4x1 Voigt vector to 2x2 symmetric tensor

      use global_parameters

      implicit none

      integer:: i
      real(kind=wp):: ATens(2,2), AVect(3,1)

      do i = 1, 2
        ATens(i,i) = AVect(i,1)
      enddo

      ATens(1,2) = AVect(3,1)
      ATens(2,1) = ATens(1,2)


      return
      end subroutine vector2symtensor2

! **********************************************************************
      subroutine vector2symtensor3(Avect,Atens)
      ! this subroutine transforms a 6x1 Voigt vector to 3x3 symmetric tensor

      use global_parameters

      implicit none

      integer:: i
      real(kind=wp):: AVect(6,1), ATens(3,3)

      do i = 1, 3
        ATens(i,i) = AVect(i,1)
      enddo

      ATens(2,3) = AVect(4,1)
      ATens(1,3) = AVect(5,1)
      ATens(1,2) = AVect(6,1)
      ATens(2,1) = ATens(1,2)
      ATens(3,1) = ATens(1,3)
      ATens(3,2) = ATens(2,3)

      return
      end subroutine vector2symtensor3

! **********************************************************************

      subroutine symtangent2matrix(C,Dmat)

      ! this subroutine maps the fourth order material/spatial tangent
      ! tensor (3x3x3x3) to a 2nd order stiffness tensor (6x6) using
      ! voigt notation: 11> 1, 22> 2, 33> 3, 23/32> 4, 13/31> 5, 12/21> 6

      use global_parameters

      implicit none

      integer:: i, j, k, l, rw, cl
      integer:: Voigt(nSymm,2)
      real(kind=wp):: C(3,3,3,3), Dmat(nSymm,nSymm)

      ! Voigt convetion: (1,1) (2,2) (3,3) (2,3) (1,3) (1,2)
      Voigt = reshape( [  1, 2, 3, 2, 1, 1,  
     &                    1, 2, 3, 3, 3, 2 ], shape(Voigt))
  
        do rw = 1, nSymm
          do cl = 1, nSymm
            i = Voigt(rw,1)
            j = Voigt(rw,2)
            k = Voigt(cl,1)
            l = Voigt(cl,2)
  
            Dmat(rw,cl) = C(i,j,k,l)
          enddo
        enddo
  
        return

        end subroutine symtangent2matrix
  
! **********************************************************************

      subroutine symtensor2vector2(ATens,AVect)
      ! this subroutine maps a symmetric tensor to a vector
      ! for unSymmmetric tensor you can use "reshape" function
  
      use global_parameters

      implicit none

      integer:: i
      real(kind=wp):: ATens(2,2), AVect(3,1)

      do i = 1, 2
        AVect(i,1) = ATens(i,i)
      enddo
      AVect(3,1) = ATens(1,2)

      return
      end subroutine symtensor2vector2
  
! **********************************************************************
  
      subroutine symtensor2vector3(ATens,AVect)
      ! this subroutine maps a symmetric tensor to a vector
      ! for unSymmmetric tensor you can use "reshape" function

      use global_parameters

      implicit none

      integer:: i
      real(kind=wp):: ATens(3,3), AVect(nSymm,1)

      do i = 1, 3
        AVect(i,1) = ATens(i,i)
      enddo

      AVect(4,1) = ATens(2,3)
      AVect(5,1) = ATens(1,3)
      AVect(6,1) = ATens(1,2)

      return
      end subroutine symtensor2vector3

! **********************************************************************