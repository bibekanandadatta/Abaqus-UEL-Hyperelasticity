************************************************************************
*********** ABAQUS/ STANDARD USER ELEMENT SUBROUTINE (UEL) *************
************************************************************************
!                   BIBEKANANDA DATTA (C) FEBRUARY 2024
!                 JOHNS HOPKINS UNIVERSITY, BALTIMORE, MD
************************************************************************
************************************************************************
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
************************************************************************
!                       LIST OF MATERIAL PROPERTIES
!
!     G           = props(1)        Shear modulus
!     kappa       = props(2)        Bulk modulus
!     lam_L       = props(3)        Locking stretch for AB model (only)
************************************************************************
!                        LIST OF ELEMENT PROPERTIES
!
!     jprops(1)   = nInt            no of integration points in element
!     jprops(2)   = matID           constitutive relation for material
!     jprops(3)   = nPostVars       no of local (int pt) post-processing variables
************************************************************************
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
************************************************************************
************************************************************************
!                          PARAMETERS MODULE
!
!     defines double precision real parameters to be used in various
!     functions and subroutines within and outside of UEL
************************************************************************


      MODULE PARAMETERS

      real*8 :: zero, one, two, three, four, five, six, eight, nine,
     &      half, third, fourth, fifth, sixth, eighth, eps, pi

      parameter( zero = 0.d0, one =  1.d0, two = 2.d0, three = 3.d0,
     &      four = 4.d0, five = 5.d0, six = 6.d0, eight = 8.d0,
     &      nine = 9.d0, half= 0.5d0, third = 1.d0/3.d0,
     &      fourth = 0.25d0, fifth = 1.d0/5.d0, sixth = 1.d0/6.d0,
     &      eighth = 0.125d0, eps = 1.d-08,
     &      pi = 3.14159265358979323846264338327950d0)

      ! no of symmetric and unSymmmetric tensor components
      integer, parameter:: nSymm = 6, nUnsymmm = 9

      ! no of total user elements and offset for dummy elements
      integer, parameter:: numElem = 50000, elemOffset = 100000
      ! numElem should be same (recommended) or higher than total
      ! number of elements in the DUMMY element set
      ! elemOffset should match with the number used in input file

      real*8, allocatable:: globalPostVars(:,:,:)

      real*8 :: ID2(2,2), ID3(3,3)
      parameter(ID2 = reshape((/1.d0,0.d0,0.d0,1.d0/),(/2,2/)),
     &          ID3 = reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,
     &                        0.d0,0.d0,1.d0/),(/3,3/)))

      END MODULE

************************************************************************
******************** ABAQUS USER ELEMENT SUBROUTINE ********************
************************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)

      USE PARAMETERS

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     & DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      real*8 :: RHS, AMATRX, SVARS, ENERGY, PNEWDT,
     &      PROPS, COORDS, DUall, Uall, Vel, Accn, TIME, DTIME,
     &      PARAMS, ADLMAG, PREDEF, DDLMAG, PERIOD

      integer:: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE,
     &      KSTEP, KINC, JELEM, NDLOAD, JDLTYP, NPREDF, LFLAGS,
     &      MLVARX, MDLOAD, JPROPS, NJPROPS

      integer:: nInt, matID, nPostVars
      integer:: nDim, ndi, nshr, ntens, uDOF, uDOFEL

      integer:: lenJobName,lenOutDir
      character*256:: outDir, fileName
      character*64 :: jobName
      character*8 :: analysis, abq_procedure
      logical:: nlgeom


      ! open a debug file for the current job
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      fileName = outDir(1:lenOutDir)//'\aaMSG_'//
     &     jobName(1:lenJobName)//'.dat'
      open(unit=80,file=fileName,status='unknown')


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
        write(80,*) 'element type is not supported', JTYPE
        write(*,*) 'element type is not supported', JTYPE
        call xit
      endif


      ! change the LFLAGS criteria as needed (check abaqus UEL manual)
      if((lflags(1).eq.1).or.(lflags(1).eq.2)) then
        ABQ_PROCEDURE = 'STATIC'
      else
        write(80,*) 'Incorrect STEP procedure in Abaqus', lflags(1)
        write(*,*) 'Incorrect STEP procedure in Abaqus', lflags(1)
        call xit
      endif

      ! check if the procedure is linear or nonlinear
      if (lflags(2).eq.0) then
        NLGEOM = .false.
      elseif (lflags(2).eq.1) then
        NLGEOM = .true.
      endif

      ! check to see if it's a general step or a linear purturbation step
      if(lflags(4).eq.1) then
        write(80,*) 'The load step should be a GENERAL step'
        write(*,*) 'The load step should be a GENERAL step'
        call xit
      endif


      ! for mixed or coupled problem, add other DOF counts as needed
      uDOF = nDim             ! displacement degrees of freedom of a node
      uDOFEL = nDim*nNode     ! total displacement degrees of freedom in element

      nInt = jprops(1)
      matID = jprops(2)
      nPostVars = jprops(3)


      ! array containing variables for post-processing
      if (.not. allocated(globalPostVars)) then
        allocate(globalPostVars(numElem,nInt,nPostVars))

        ! print necessary information to the screen
        write(80,*) '---------------------------------------'
        write(80,*) '------- ABAQUS FINITE STRAIN UEL ------'
        write(80,*) '---------------------------------------'
        write(80,*) jobName
        write(80,*) '---------------------------------------'
        write(80,*) '------- PROCEDURE = ', ABQ_PROCEDURE
        write(80,*) '------- ANALYSIS TYPE   = ', analysis
        write(80,*) '---------- NLGEOM = ', NLGEOM
        write(80,*) '------- MODEL DIMENSION = ', nDim
        write(80,*) '------- ELEMENT NODES   = ', NNODE
        write(80,*) '------- MATERIAL ID     = ', matID
        write(80,*) '---------------------------------------'
        write(80,*) '-------- INTEGRATION SCHEME -----------'
        write(80,*) '----------- NINT  = ', nInt
        write(80,*) '---------------------------------------'
        write(80,*) '---------- POST-PROCESSING ------------'
        write(80,*) '--- NO OF ELEMENTS            = ', numElem
        write(80,*) '--- DUMMY ELEMENT OFFSET      = ', ElemOffset
        write(80,*) '--- NO OF VARIABLES AT INT PT = ', nPostVars
        write(80,*) '---------------------------------------'

        ! print necessary information to the screen
        write(*,*) '---------------------------------------'
        write(*,*) '------- ABAQUS FINITE STRAIN UEL ------'
        write(*,*) '---------------------------------------'
        write(*,*) jobName
        write(*,*) '---------------------------------------'
        write(*,*) '------- PROCEDURE = ', ABQ_PROCEDURE
        write(*,*) '------- ANALYSIS TYPE   = ', analysis
        write(*,*) '---------- NLGEOM = ', NLGEOM
        write(*,*) '------- MODEL DIMENSION = ', nDim
        write(*,*) '------- ELEMENT NODES   = ', NNODE
        write(*,*) '------- MATERIAL ID     = ', matID
        write(*,*) '---------------------------------------'
        write(*,*) '-------- INTEGRATION SCHEME -----------'
        write(*,*) '----------- NINT  = ', nInt
        write(*,*) '---------------------------------------'
        write(*,*) '---------- POST-PROCESSING ------------'
        write(*,*) '--- NO OF ELEMENTS            = ', numElem
        write(*,*) '--- DUMMY ELEMENT OFFSET      = ', ElemOffset
        write(*,*) '--- NO OF VARIABLES AT INT PT = ', nPostVars
        write(*,*) '---------------------------------------'

      endif

       ! call your UEL subroutine
       call uel_NLMECH(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & ANALYSIS,NDIM,NDI,NSHR,NTENS,NINT,UDOF,UDOFEL)


      RETURN
      END SUBROUTINE UEL

************************************************************************
************************************************************************

      SUBROUTINE uel_NLMECH(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & ANALYSIS,NDIM,NDI,NSHR,NTENS,NINT,UDOF,UDOFEL)

      USE PARAMETERS

     !!!!!!!!!!!!! VARIABLE DECLARATION AND INITIALIZATION !!!!!!!!!!!!!
     
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     & DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      real*8 :: RHS, AMATRX, SVARS, ENERGY, PNEWDT,
     &      PROPS, COORDS, DUall, Uall, Vel, Accn, TIME, DTIME,
     &      PARAMS, ADLMAG, PREDEF, DDLMAG, PERIOD

      integer:: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE,
     &      KSTEP, KINC, JELEM, NDLOAD, JDLTYP, NPREDF, LFLAGS,
     &      MLVARX, MDLOAD, JPROPS, NJPROP

      character*8:: analysis
      integer:: nDim, ndi, nshr, ntens, nInt, uDOF, uDOFEL
      logical:: nlgeom

      real*8 :: uNode(nDim,nNode), duNode(nDim,nNode), F(3,3),
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

      real*8 :: stranLagrange(ntens,1), stranEuler(ntens,1),
     &      stressCauchy(ntens,1), stressPK1(nDim*nDim,1),
     &      stressPK2(ntens,1), Dmat(ntens,ntens)

      integer:: i, j, intPt, matID, istat

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
      RHS(1:MLVARX,1) = zeros

      matID = jprops(2)
      nlocalSdv = NSVARS/NINT

     !!!!!!!!!!! END VARIABLE DECLARATION AND INITIALIZATION !!!!!!!!!!!

      ! reshape the displacement vectors into matrix forms
      uNode  = reshape(UAll,(/nDim,nNode/))
      duNode = reshape(DUAll(:,1),(/nDim,nNode/))

      ! if applicable gather the prescribed field variables in a vector
      ! such as temperature (as shown below in commented line - not tested)
      ! fieldNode(1,1:nNode) = predef(1,1,1:nNode)
      ! dfieldNode(1,1:nNode) = predef(2,1,1:nNode)


     !!!!!!!!!!!!!!!!! ELEMENT RELATED OPERATIONS !!!!!!!!!!!!!!!!!!!!!!

      ! obtain gauss quadrature points and weight
      if (nDim.eq.2) then
        ID = ID2
        call gaussQuadrtr2(nNode,nInt,w,xi)
      elseif (nDim.eq.3) then
        ID = ID3
        call gaussQuadrtr3(nNode,nInt,w,xi)
      else
        write(80,*) 'incorrect model dimension, nDim =', nDim
        write(*,*) 'incorrect model dimension, nDim =', nDim
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
          write(80,*) 'incorrect model dimension, nDim =', nDim
          write(*,*) 'incorrect model dimension, nDim =', nDim
          call xit
        endif

        ! calculate element jacobian and global gradient
        dXdxi   = matmul(coords,dNdxi)      ! calculate dXdxi

        if (nDim.eq.2) then
          call inverseMat2(dXdxi,dxidX,detJ,istat)
        elseif(nDim.eq.3) then
          call inverseMat3(dXdxi,dxidX,detJ,istat)
        endif

        if (istat .eq. 0) then
          write(80,*) 'ill-condiitoned element jacobian', jElem, intPt
          write(*,*) 'ill-condiitoned element jacobian', jElem, intPt
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
            Ba(3,1) = dNdx(i,2)
            Ba(3,2) = dNdx(i,1)

          ! form [Ba] matrix: 3D case
          elseif (analysis.eq.'3D') then
            Ba(1,1) = dNdx(i,1)
            Ba(2,2) = dNdx(i,2)
            Ba(3,3) = dNdx(i,3)
            Ba(4,2) = dNdx(i,3)
            Ba(4,3) = dNdx(i,2)
            Ba(5,1) = dNdx(i,3)
            Ba(5,3) = dNdx(i,1)
            Ba(6,1) = dNdx(i,2)
            Ba(6,2) = dNdx(i,1)

          else
            write(80,*) 'wrong analysis type: ', analysis
            write(*,*) 'wrong analysis type: ', analysis
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
     &          coords,nnode,kstep,kinc,props,nprops,nlocalSdv,
     &          analysis)

        elseif (matID .eq. 2) then
          call umatArrudaBoyce(stressCauchy,stressPK1,stressPK2,
     &          Dmat,F,svars,nsvars,stranLagrange,stranEuler,time,
     &          dtime,fieldVar,npredf,nDim,ndi,nshr,ntens,jelem,intPt,
     &          coords,nnode,kstep,kinc,props,nprops,nlocalSdv,
     &          analysis)
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


      RETURN
      END SUBROUTINE uel_NLMECH

************************************************************************
************************************************************************

      SUBROUTINE umatNeoHookean(stressCauchy,stressPK1,stressPK2,
     &          Dmat,F,svars,nsvars,stranLagrange,stranEuler,time,
     &          dtime,fieldVar,npredf,nDim,ndi,nshr,ntens,jelem,intpt,
     &          coords,nnode,kstep,kinc,props,nprops,nlocalSdv,
     &          analysis)

      USE PARAMETERS

      IMPLICIT NONE

      integer:: nsvars, npredf, nDim, ndi, nshr, ntens,
     &    jelem, intpt, nNode, kstep, kinc, nprops

      real*8 :: stressCauchy(ntens,1), stressPK1(nDim*nDim,1),
     &    stressPK2(ntens,1), Dmat(ntens,ntens), F(3,3),
     &    stranLagrange(ntens,1), stranEuler(ntens,1), props(1:nprops),
     &    svars(1:nsvars), coords(nDim,nNode), time(2), dtime,
     &    fieldVar(npredf)

      character*8:: analysis

      real*8 :: detF, C(3,3), Cinv(3,3), detC, B(3,3), Binv(3,3), detB,
     &    strantensorEuler(3,3), stressTensorPK2(3,3),
     &    stressTensorCauchy(3,3), Cmat(3,3,3,3), VoigtMat(nSymm,nSymm),
     &    stranVoigtEuler(nSymm,1), stressVoigtPK1(nUnsymmm,1),
     &    stressVoigtPK2(nSymm,1), stressVoigtCauchy(nSymm,1),
     &    Gshear, kappa, lam_L

      integer:: nInt, nlocalSdv

      ! loop counters
      integer:: i, j, k, l, istat

      ! initialize matrial stiffness tensors
      Cmat   = zero
      Dmat   = zero

      ! assign material properties to variables
      Gshear= props(1)        ! Shear modulus
      kappa = props(2)        ! Bulk modulus
      lam_L = props(3)        ! locking stretch

      ! locking stretch should be infinity (0 as input) for NH model
      if (lam_L .ne. zero) then
        write(80,*) 'Incorrect material parameter (lam_L)', lam_L
        write(*,*) 'Incorrect material parameter (lam_L)', lam_L
        call xit
      endif

      ! perform all the constitutitve relations in 3D
      call detMat3(F,detF)
      if (detF .lt. zero) then
        write(80,*) 'check result: detF.lt.zero', jelem, intpt, detF
        write(*,*) 'check result: detF.lt.zero', jelem, intpt, detF
      endif

      B = matmul(F,transpose(F))
      C = matmul(transpose(F),F)

      call inverseMat3(B,Binv,detB,istat)
      call inverseMat3(C,Cinv,detC,istat)

      ! calculate Euler-Almansi strain tensor
      strantensorEuler = half*(ID3-Binv)

      ! calculate material tangent, C_ijkl
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              Cmat(i,j,k,l) = Cmat(i,j,k,l)
     &            + kappa*Cinv(i,j)*Cinv(k,l)
     &            + (Gshear-kappa*dlog(detF)) *( Cinv(i,k)*Cinv(j,l)
     &            + Cinv(i,l)*Cinv(j,k) )
            enddo
          enddo
        enddo
      enddo

      ! calculate stress tensors
      stressTensorCauchy = (1/detF)*(Gshear*(B-ID3) +
     & 											kappa*dlog(detF)*ID3)
      stressTensorPK2 = Gshear*(ID3-Cinv) + kappa*dlog(detF)*Cinv


      ! transforms the stiffness tensor 3x3x3x3 to a 6x6 matrix
      call tangent2matrix(Cmat,VoigtMat)

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


      RETURN

      END SUBROUTINE umatNeoHookean

************************************************************************
************************************************************************
      SUBROUTINE umatArrudaBoyce(stressCauchy,stressPK1,stressPK2,
     &          Dmat,F,svars,nsvars,stranLagrange,stranEuler,time,
     &          dtime,fieldVar,npredf,nDim,ndi,nshr,ntens,jelem,intpt,
     &          coords,nnode,kstep,kinc,props,nprops,nlocalSdv,
     &          analysis)

      USE PARAMETERS

      IMPLICIT NONE

      integer:: nsvars, npredf, nDim, ndi, nshr, ntens,
     &    jelem, intpt, nNode, kstep, kinc, nprops

      real*8 :: stressCauchy(ntens,1), stressPK1(nDim*nDim,1),
     &    stressPK2(ntens,1), Dmat(ntens,ntens), F(3,3),
     &    stranLagrange(ntens,1), stranEuler(ntens,1), props(1:nprops),
     &    svars(1:nsvars), coords(nDim,nNode), time(2), dtime,
     &    fieldVar(npredf)

      character*8:: analysis

      real*8 :: detF, C(3,3), Cinv(3,3), detC, B(3,3), Binv(3,3), detB,
     & 		trC, lam_c, lam_r, beta_c, dBeta_c,
     &    strantensorEuler(3,3), stressTensorPK2(3,3),
     &    stressTensorCauchy(3,3), Cmat(3,3,3,3), VoigtMat(nSymm,nSymm),
     &    stranVoigtEuler(nSymm,1), stressVoigtPK1(nUnsymmm,1),
     &    stressVoigtPK2(nSymm,1), stressVoigtCauchy(nSymm,1),
     &    Gshear, kappa, lam_L

      integer:: nInt, nlocalSdv
      ! loop counters
      integer:: i, j, k, l, istat

      ! initialize matrial stiffness tensors
      Cmat   = zero
      Dmat   = zero

      ! assign material properties to variables
      Gshear= props(1)        ! Shear modulus
      kappa = props(2)        ! Bulk modulus
      lam_L = props(3) 				! Locking stretch

      if (lam_L .eq. zero) then
        write(80,*) 'Incorrect material parameter for AB model', lam_L
        write(*,*) 'Incorrect material parameter for AB model', lam_L
        call xit
      endif

      ! perform all the constitutitve relations in 3D
      call detMat3(F,detF)
      if (detF .lt. zero) then
        write(*,*) 'Check result: detF.lt.zero', jelem, intPt, detF
        write(*,*) 'Check result: detF.lt.zero', jelem, intPt, detF
      endif

      B = matmul(F,transpose(F))
      C = matmul(transpose(F),F)

      call inverseMat3(B,Binv,detB,istat)
      call inverseMat3(C,Cinv,detC,istat)

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
              Cmat(i,j,k,l) = Cmat(i,j,k,l) + Gshear/(nine*lam_c**2)*
     &              (dBeta_c- lam_r*beta_c)*ID3(i,j)*ID3(k,l)
     &            + kappa*Cinv(i,j)*Cinv(k,l)
     &            + ( Gshear/3.0*lam_r*beta_c-kappa*dlog(detF) )*
     &              ( Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k) )
            enddo
          enddo
        enddo
      enddo

      stressTensorCauchy = (1/detF)*( (Gshear/three)*lam_r*beta_c*B -
     &      (Gshear*lam_L/three - kappa*dlog(detF))*ID3 )
      stressTensorPK2 = Gshear/three*lam_r*beta_c*ID3 -
     &      (Gshear*lam_L/three - kappa*dlog(detF))*Cinv


      ! transforms the stiffness tensor 3x3x3x3 to a 6x6 matrix
      call tangent2matrix(Cmat,VoigtMat)

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
     & 																		stranEuler(1:ntens,1)

      RETURN

      contains

      ! approximation of inverse Langevin function
      ! reference: Bergstorm (PhD thesis, MIT, 1999)
      FUNCTION InvLangevin(x)

      IMPLICIT NONE
      real*8,intent(in) :: x
      real*8 :: InvLangevin


      if (dabs(x) .lt. 0.84136d0) then
        InvLangevin = 1.31446*dtan(1.58986d0*x) + 0.91209d0*x
      elseif ((dabs(x) .ge. 0.84136d0) .and. (dabs(x) .lt. one)) then
        InvLangevin = one/(dsign(one,x)-x)
      else
        write(80,*) 'unbound argument for inverse Langevin function', x
        write(*,*) 'unbound argument for inverse Langevin function', x
        call xit
      endif

      return
      END FUNCTION InvLangevin


      ! derivative of inverse Langevin function
      FUNCTION DInvLangevin(x)

      IMPLICIT NONE

      real*8,intent(in) :: x
      real*8 :: DInvLangevin, sec

      if (dabs(x) .lt. 0.84136d0) then
        DInvLangevin = 2.0898073756d0*(dtan(1.58986d0*x))**two 
     &                + 3.0018973756d0
      elseif ((dabs(x) .ge. 0.84136) .and. (dabs(x) .lt. one)) then
        DInvLangevin = one/((dsign(one,x)-x)**two)
      else
        write(80,*) 'unbound argument for inverse Langevin function', x
        write(*,*) 'unbound argument for inverse Langevin function', x
        call xit
      endif

      return
      END FUNCTION DInvLangevin

      END SUBROUTINE umatArrudaBoyce

************************************************************************
************************************************************************

       SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     & NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     & JMAC,JMATYP,MATLAYO,LACCFLA)
       ! this subroutine is used to transfer postVars from the UEL
       ! onto the dummy mesh for viewing. Note that an offset of
       ! elemOffset is used between the real mesh and the dummy mesh.

      USE PARAMETERS

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      ! the dimensions of the variables FLGRAY, ARRAY and JARRAY
      ! must be set equal to or greater than 15.  

      ! explicityly define the type for uvar to avoid issues
      real*8 :: uvar

      uvar(1:nuvarm) = globalPostVars(noel-elemOffset,npt,1:nuvarm)

      RETURN

      END SUBROUTINE UVARM

************************************************************************
************************************************************************




************************************************************************
!                      ELEMENT UTILITY SECTION
!
!     SUBROUTINE to evaluate shape functions at the gauss pts
!     SUBROUTINE to calculate Gauss Pts in 2D and 3D elements
!     ----- currently supports:    (a) tri elements   (3 and 6 nodes) : full integration
!                                  (b) quad elements  (4 and 8 nodes) : reduced and full integration
!                                  (c) tet elements   (4 and 10 nodes): full integration
!                                  (d) hex elements   (8 and 20 nodes): reduced and full integration
************************************************************************

      SUBROUTINE interpFunc2(nNode,nInt,intPt,xi_int,Nxi,dNdxi)

      ! this subroutine calculates shape function of 2D elements at gauss pts
      ! available 2D elements are: 3 node tri, 6 node tri, 4 node quad, 8 node quad
        
      ! Nxi(i)          = shape function of node i at the intpt.
      ! dNdxi(i,j)      = derivative wrt j direction of shape fn of node i

      USE PARAMETERS

      integer:: nNode, nInt, intpt

      real*8 :: xi_int(nInt,2), Nxi(nNode), dNdxi(nNode,2)
      real*8 :: xi, eta, lam

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
        write(80,*) 'element is not supported for 2D analysis', nNode
        write(*,*) 'element is not supported for 2D analysis', nNode
        call xit
      endif

      RETURN
      END SUBROUTINE interpFunc2

************************************************************************

      subroutine interpFunc3(nNode,nInt,intPt,xi_int,Nxi,dNdxi)

      ! this subroutine calculates shape function of 3D elements at gauss pts
      ! available 3D elements are: 4 node tet, 10 node ter, 8 node hex.

      ! Nxi(i)          = shape function of node i at the intpt.
      ! dNdxi(i,j)      = derivative wrt j direction of shape fn of node i

      USE PARAMETERS

      IMPLICIT NONE

      integer:: nNode, nInt, intpt

      real*8 :: xi_int(nInt,3), Nxi(nNode), dNdxi(nNode,3)
      real*8 :: xi, eta, zeta, lam

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
        Nxi(1) = xi
        Nxi(2) = eta
        Nxi(3) = zeta
        Nxi(4) = one-xi-eta-zeta
        ! the first derivatives of the shape functions dN/dxi (4x3)
        dNdxi(1,1) = one
        dNdxi(2,2) = one
        dNdxi(3,3) = one
        dNdxi(4,1) = -one
        dNdxi(4,2) = -one
        dNdxi(4,3) = -one

      elseif (nNode.eq.10) then  ! 10-noded quadratic tet10 element
    
        ! shape functions
        lam = one-xi-eta-zeta
        Nxi(1) = (two*xi-one)*xi
        Nxi(2) = (two*eta-one)*eta
        Nxi(3) = (two*zeta-one)*zeta
        Nxi(4) = (two*lam-one)*lam
        Nxi(5) = four*xi*eta
        Nxi(6) = four*eta*zeta
        Nxi(7) = four*zeta*xi
        Nxi(8) = four*xi*lam
        Nxi(9) = four*eta*lam
        Nxi(10) = four*zeta*lam

        dNdxi(1,1) = (four*xi-one)
        dNdxi(2,2) = (four*eta-one)
        dNdxi(3,3) = (four*zeta-one)
        dNdxi(4,1) = -(four*lam-one)
        dNdxi(4,2) = -(four*lam-one)
        dNdxi(4,3) = -(four*lam-one)
        dNdxi(5,1) = four*eta
        dNdxi(5,2) = four*xi
        dNdxi(6,2) = four*zeta
        dNdxi(6,3) = four*eta
        dNdxi(7,1) = four*zeta
        dNdxi(7,3) = four*xi
        dNdxi(8,1) = four*(lam-xi)
        dNdxi(8,2) = -four*xi
        dNdxi(8,3) = -four*xi
        dNdxi(9,1) = -four*eta
        dNdxi(9,2) = four*(lam-eta)
        dNdxi(9,3) = -four*eta
        dNdxi(10,1) = -four*zeta*lam
        dNdxi(10,2) = -four*zeta
        dNdxi(10,3) = four*(lam-zeta)

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
        write(80,*) 'element is not supported for 3D analysis', nNode
        write(*,*) 'element is not supported for 3D analysis', nNode
        call xit
      endif

      RETURN
      END SUBROUTINE interpFunc3

************************************************************************

      SUBROUTINE gaussQuadrtr2(nNode,nInt,w,xi)

      ! this subroutines returns the weights and
      ! gauss point coordinate of 2D Lagrangian elements
      ! currently supports: nInt = 1 (tri3) and 3 (tri6)
      !                     nInt = 1, 4 (quad4) and 4, 9 (quad8)

      USE PARAMETERS

      IMPLICIT NONE
      integer:: nNode, nInt
      real*8 :: x1D(4), w1D(4)
      real*8 :: w(nInt), xi(nInt,2)

      w  = zero
      xi = zero

      if (nNode.eq.3) then      ! plane tri3 elements (full integration)
        if (nInt.eq.1) then
          w(1) = half
          xi(1,1) = third
          xi(2,1) = third
        else
          write(80,*) 'wrong gauss points for tri3 element', nInt
          write(*,*) 'wrong gauss points for tri3 element', nInt
          call xit
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
          write(80,*) 'wrong gauss points for tri6 element', nInt
          write(*,*) 'wrong gauss points for tri6 element', nInt
          call xit
        endif
        
      elseif((nNode.eq.4)) then ! plane quad4 element

        if (nInt.eq.1) then     ! reduced integration for quad4
          w = four

          xi(1,1) = zero
          xi(1,2) = zero

        elseif (nInt.eq.4) then ! full integration for quad4

          w(1:4) = one

          x1D(1) = dsqrt(third)
          xi(1,1) = -x1D(1)
          xi(1,2) = -x1D(1)
          xi(2,1) = x1D(1)
          xi(2,2) = -x1D(1)
          xi(3,1) = -x1D(1)
          xi(3,2) = x1D(1)
          xi(4,1) = x1D(1)
          xi(4,2) = x1D(1)
        else
          write(80,*) 'wrong gauss points for quad4 element', nInt
          write(*,*) 'wrong gauss points for quad4 element', nInt
          call xit
        endif

      elseif (nNode.eq.8) then  ! plane quad8 element
        
        if (nInt.eq.4) then     ! reduced integration for quad8

          w(1:4) = one

          x1D(1) = dsqrt(third)
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

          x1D(1) = dsqrt(three/five)
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
          write(80,*) 'wrong gauss points for quad8 element', nInt
          write(*,*) 'wrong gauss points for quad8 element', nInt
          call xit
        endif

      else
        write(80,*) 'elements not supported for 2D analysis', nNode
        write(*,*) 'elements not supported for 2D analysis', nNode
        call xit
      endif


      RETURN
      END SUBROUTINE gaussQuadrtr2

************************************************************************

      SUBROUTINE gaussQuadrtr3(nNode,nInt,w,xi)

      ! this subroutines returns the weights and
      ! gauss point coordinate of 3D Lagrangian elements
      ! currently supports: nInt = 1 (tet4) and 4 (tet10)
      !                     nInt = 1, 8 (hex8) and 8, 27 (hex20)

      USE PARAMETERS

      IMPLICIT NONE

      integer:: nNode, nInt
      integer:: i, j, k, n
      real*8 :: x1D(4), w1D(4)
      real*8 :: w(nInt), xi(nInt,3)

      w  = zero
      xi = zero
      
      if(nNode.eq.4) then       ! 3D tet4 element (full integration)
        if(nInt.eq.1) then
          w(1) = sixth
          xi(1:3,1) = fourth

        else
          write(80,*) 'wrong gauss points for tet4 element', nInt
          write(*,*) 'wrong gauss points for tet4 element', nInt
          call xit
        endif

      elseif(nNode.eq.10) then  ! 3D tet10 element (full integration)

        if (nInt.eq.4) then
          w(1:4) = one/24.d0

          xi(1,1) = 0.58541020d0
          xi(2,1) = 0.13819660d0
          xi(3,1) = xi(2,1)
          xi(1,2) = xi(2,1)
          xi(2,2) = xi(1,1)
          xi(3,2) = xi(2,1)
          xi(1,3) = xi(2,1)
          xi(2,3) = xi(2,1)
          xi(3,3) = xi(1,1)
          xi(4,1) = xi(2,1)
          xi(4,2) = xi(2,1)
          xi(4,3) = xi(2,1)

        else
          write(80,*) 'wrong gauss points for tet10 element', nInt
          write(*,*) 'wrong gauss points for tet10 element', nInt
          call xit
        endif

      elseif(nNode.eq.8) then   ! 3D hex8 element

        if(nInt.eq.1) then      ! reduced integration for hex8
          w(1) = eight
          xi(1,1:3) = zero

        elseif(nInt.eq.8) then  ! full-integration for hex8
          w(1:8) = one

          x1D(1) = -dsqrt(third)
          x1D(2) = dsqrt(third)
         
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
          write(80,*) 'wrong gauss points for hex8 element', nInt
          write(*,*) 'wrong gauss points for hex8 element', nInt
          call xit
        endif

      elseif(nNode.eq.20) then  ! 3D hex20 element
        
        if (nInt.eq.8) then     ! reduced integration for hex20
          w(1:8) = one

          x1D(1) = -dsqrt(third)
          x1D(2) = dsqrt(third)
         
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

          x1D(1) = -dsqrt(0.6d0)
          x1D(2) = zero
          x1D(3) = dsqrt(0.6d0)
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
          write(80,*) 'wrong gauss points for hex20 element', nInt
          write(*,*) 'wrong gauss points for hex20 element', nInt
          call xit
        endif

      else
        write(80,*) 'element is not supported for 3D analysis', nNode
        write(*,*) 'element is not supported for 3D analysis', nNode
        call xit
      endif

      RETURN
      END SUBROUTINE gaussQuadrtr3


************************************************************************

      SUBROUTINE faceNodes(nDim,nNode,face,list,nFaceNodes)
      ! this subroutine RETURNs the list of nodes on an
      ! element face for standard 2D and 3D Lagrangian elements

      IMPLICIT NONE

      integer, intent (in) :: nDim, nNode, face
      integer, intent (out):: list(*)
      integer, intent (out):: nFaceNodes

      integer :: list2(3), list3(4)

      if (nDim.eq.2) then
        list2(1:3) = [2,3,1]
        list3(1:4) = [2,3,4,1]

        if (nNode.eq.3) then
          nFaceNodes = 2
          list(1) = face
          list(2) = list2(face)
        else if (nNode.eq.4) then
          nFaceNodes = 2
          list(1) = face
          list(2) = list3(face)
        else if (nNode.eq.6) then
          nFaceNodes = 3
          list(1) = face
          list(2) = list2(face)
          list(3) = face+3
        else if (nNode.eq.8) then
          nFaceNodes = 3
          list(1) = face
          list(2) = list3(face)
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

      RETURN

      END SUBROUTINE faceNodes

************************************************************************
*********************** MATRIX ALGEBRA SECTION *************************
************************************************************************
!
!     SUBROUTINE to create identity matrix of any dimention
!     SUBROUTINE to calculate determinant of matrix
!     SUBROUTINE to calculate direct inverse of 2x2 and 3x3 matrix
!     SUBROUTINE to map symmetric tensor to a vector
!     SUBROUTINE to map 4th order tensor to 2D Voigt matrix
************************************************************************

      SUBROUTINE crossProduct(a,b,c)
      ! this subroutine computes the cross product of two 3 dimensional vectors

      IMPLICIT NONE

      real*8 :: a(3), b(3), c(3)

      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = b(1)*a(3)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)

      END SUBROUTINE crossProduct

************************************************************************

      SUBROUTINE traceMat(A,trA,nDim)
      ! this subroutine calculates the trace of a square matrix

      USE PARAMETERS
      IMPLICIT NONE

      integer:: nDim, i
      real*8 :: A(nDim,nDim), trA

      trA = zero

      do i = 1, nDim
        trA = trA + A(i,i)
      enddo

      RETURN
      END SUBROUTINE traceMat

************************************************************************

      SUBROUTINE detMat2(A,detA)
      ! this subroutine calculates the determinant of a 2x2 or 3x3 matrix [A]

      IMPLICIT NONE

      real*8 :: A(2,2), detA

      detA = A(1,1)*A(2,2) - A(1,2)*A(2,1)

      RETURN
      END SUBROUTINE detMat2
************************************************************************

      SUBROUTINE detMat3(A,detA)
      ! this subroutine calculates the determinant of a 2x2 or 3x3 matrix [A]

      IMPLICIT NONE

      real*8 :: A(3,3), detA

      detA = A(1,1)*A(2,2)*A(3,3)
     &     + A(1,2)*A(2,3)*A(3,1)
     &     + A(1,3)*A(2,1)*A(3,2)
     &     - A(3,1)*A(2,2)*A(1,3)
     &     - A(3,2)*A(2,3)*A(1,1)
     &     - A(3,3)*A(2,1)*A(1,2)

      RETURN
      END SUBROUTINE detMat3

************************************************************************

      SUBROUTINE inverseMat2(A,Ainv,detA,istat)
      ! this subroutine returns Ainv and detA for a 2D matrix A

      USE PARAMETERS

      IMPLICIT NONE

      integer:: istat
      real*8 :: A(2,2),Ainv(2,2), detA, detAinv

      istat = 1

      call detMat2(A,detA)

      if (detA .le. zero) then
        write(80,*) 'WARNING: subroutine inverseMat2:'
        write(80,*) 'WARNING: det of mat= ', detA
        write(*,*) 'WARNING: subroutine inverseMat2:'
        write(*,*) 'WARNING: det of mat= ', detA
        istat = 0
        RETURN
      end if

      detAinv = one/detA

      Ainv(1,1) =  detAinv*A(2,2)
      Ainv(1,2) = -detAinv*A(1,2)
      Ainv(2,1) = -detAinv*A(2,1)
      Ainv(2,2) =  detAinv*A(1,1)

      RETURN

      END SUBROUTINE inverseMat2
************************************************************************

      SUBROUTINE inverseMat3(A,Ainv,detA,istat)
      ! this subroutine returns Ainv and detA for a 3D matrix A

      USE PARAMETERS

      IMPLICIT NONE

      integer:: istat
      real*8 :: A(3,3),Ainv(3,3), detA, detAinv

      istat = 1

      call detMat3(A,detA)

      if (detA .le. zero) then
        write(80,*) 'WARNING: subroutine inverseMat3:'
        write(80,*) 'WARNING: det of mat= ', detA
        write(*,*) 'WARNING: subroutine inverseMat3:'
        write(*,*) 'WARNING: det of mat= ', detA
        istat = 0
        RETURN
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

      RETURN
      END SUBROUTINE inverseMat3
************************************************************************

      subroutine inverseMat(A,Ainv,n)
      ! this subroutine computes the inverse of an arbitrary
      ! square matrix (nxn) by LU decomposition approach

      USE PARAMETERS

      IMPLICIT NONE

      integer,intent(in)   :: n

      real*8,intent(inout) :: A(n,n)
      real*8,intent(out)   :: Ainv(n,n)

      real*8 :: L(n,n), U(n,n), b(n), d(n), x(n)
      real*8 :: coeff
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

      END SUBROUTINE inverseMat

************************************************************************

      subroutine eigenSym3(A,eigenvalues,eigenvectors)
      ! this subroutine computes eigenvals and eigenvectors of symmetric 3x3 matrix

      USE PARAMETERS

      IMPLICIT NONE

      real*8, intent(in)  :: A(3,3)                   ! input matrix
      real*8, intent(out) :: eigenvalues(3)           ! eigenvalues
      real*8, intent(out) :: eigenvectors(3,3)        ! ith eigenvector is eigenvectors(1:3,i)

      real*8 :: B(3,3), C(3,3), D(3,3)
      real*8 :: p1, p2, p, q, r, phi, evnorm, tol
      integer:: i

      eigenvectors = zero

      p1 = A(1,2)*A(1,2) + A(1,3)*A(1,3) + A(2,3)*A(2,3)
      tol = eps*(A(1,1)*A(1,1) + A(2,2)*A(2,2) + A(3,3)*A(3,3) + p1)

      if (p1 == zero) then
        eigenvalues(1) = A(1,1)
        eigenvalues(2) = A(2,2)
        eigenvalues(3) = A(3,3)

        eigenvectors(1,1) = one
        eigenvectors(2,2) = one
        eigenvectors(3,3) = one

      else
          q = (A(1,1)+A(2,2)+A(3,3))/three
          p2 = (A(1,1) - q)*(A(1,1)-q) + (A(2,2) - q)*(A(2,2) - q)
     &         + (A(3,3) - q)*(A(3,3) - q) + two * p1
          p = dsqrt(p2 / six)
          B = (one / p) * (A - q * ID3)
          r =   half*( B(1,1)*B(2,2)*B(3,3)
     &                - B(1,1)*B(2,3)*B(3,2)
     &                - B(1,2)*B(2,1)*B(3,3)
     &                + B(1,2)*B(2,3)*B(3,1)
     &                + B(1,3)*B(2,1)*B(3,2)
     &                - B(1,3)*B(2,2)*B(3,1) )

            if (r < -one) then
                phi = pi / three
            else if (r > one) then
                phi = zero
            else
                phi = dacos(r) / three
            end if

            ! the eigenvalues satisfy eig3 <= eig2 <= eig1
            eigenvalues(1) = q + two*p*dcos(phi)
            eigenvalues(3) = q + two*p*dcos(phi + (two*pi/three))
            eigenvalues(2) = three*q - eigenvalues(1) - eigenvalues(3)

            do i = 1,3
              B = A - eigenvalues(i)*ID3
              C = A - eigenvalues(mod(i,3)+1)*ID3
              D = matmul(B,C)
              eigenvectors(1:3,mod(i+1,3)+1) = matmul(D,(/one,one,one/))
            end do

            do i = 1,3
              evnorm = dsqrt(dot_product(eigenvectors(1:3,i),
     &                 eigenvectors(1:3,i)))
              if (evnorm>tol) then
                  eigenvectors(1:3,i) = eigenvectors(1:3,i)/evnorm
              else
                  call crossProduct(
     &            eigenvectors(1:3,mod(i,3)+1),
     &            eigenvectors(1:3,mod(i+1,3)+1),
     &            eigenvectors(1:3,i))

                  evnorm = dsqrt(dot_product(eigenvectors(1:3,i),
     &                        eigenvectors(1:3,i)))

                  eigenvectors(1:3,i) = eigenvectors(1:3,i)/evnorm
              endif
            end do
      end if

      end subroutine eigenSym3

************************************************************************

      SUBROUTINE sqrtMat3(A, B)
      ! this subroutines computes square root of a symmetric 3x3 matrix

      USE PARAMETERS

      IMPLICIT NONE

      real*8, intent(in) :: A(3,3)

      real*8 :: D(3,3), V(3,3), eig(3)
      real*8 :: B(3,3)

      call eigenSym3(A,eig,V)

      D = zero

      D(1,1) = dsqrt(eig(1))
      D(2,2) = dsqrt(eig(2))
      D(3,3) = dsqrt(eig(3))

      B = matmul(V,matmul(D,transpose(V)))

      END SUBROUTINE sqrtMat3

************************************************************************

      SUBROUTINE polarDecomp3(A,V,U,R)
      ! this subroutine compute left and right polar decompositions
      ! of a 3x3 matrix A (used for deformation gradient)

      USE PARAMETERS

      IMPLICIT NONE

      real*8, intent (in)   :: A(3,3)
      real*8, intent (out)  :: R(3,3), U(3,3), V(3,3)
      real*8 :: Vinv(3,3), detV
      integer:: nDim, istat

      !  Decompose A into A=RU=VR  where U,V are symmetric and R is orthogonal

      R = matmul(A,transpose(A))                      ! R is just temporary variable here
      call sqrtMat3(R,V)                              ! V= sqrt(A*A^T)
      call inverseMat3(V,Vinv,detV,istat)
      R = matmul(Vinv,A)                              ! R = V^-1*A
      U = matmul(transpose(R),A)                      ! U = R^T*A

      end subroutine polarDecomp3
************************************************************************

      SUBROUTINE voigtAugment(vect2D, vect3D, ntens)
      ! this subroutine augemnts a 3x1 (plane) or 4x1 (axisymmetric)
      ! array to a 6x1 Voigt array of 3D dimensional case

      USE PARAMETERS

      IMPLICIT NONE

      integer:: ntens
      real*8 :: vect2D(ntens,1), vect3D(nSymm,1)

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

      RETURN
      END SUBROUTINE voigtAugment

************************************************************************

      SUBROUTINE voigtTruncate(vect3D, vect2D, ntens)
      ! this subroutine truncates a 6x1 Voigt array
      ! to a 3x1 (plane) or 4x1 (axisymmetry) Voigt array

      USE PARAMETERS

      IMPLICIT NONE

      integer:: ntens
      real*8 :: vect3D(nSymm,1), vect2D(ntens,1)

      vect2D = zero

      vect2D(1,1) = vect3D(1,1)
      vect2D(2,1) = vect3D(2,1)

      if (ntens.eq. 3) then
        vect2D(3,1) = vect3D(6,1)
      elseif (ntens.eq.4) then
        vect2D(3,1) = vect3D(3,1)
        vect2D(4,1) = vect3D(6,1)
      endif

      RETURN

      END SUBROUTINE voigtTruncate

************************************************************************
      SUBROUTINE symtensor2vector2(ATens,AVect)
      ! this subroutine maps a symmetric tensor to a vector
      ! for unSymmmetric tensor you can use "reshape" function

      USE PARAMETERS

      IMPLICIT NONE

      integer:: i
      real*8 :: ATens(2,2), AVect(3,1)

      do i = 1, 2
        AVect(i,1) = ATens(i,i)
      enddo
      AVect(3,1) = ATens(1,2)

      RETURN
      END SUBROUTINE symtensor2vector2

************************************************************************

      SUBROUTINE symtensor2vector3(ATens,AVect)
      ! this subroutine maps a symmetric tensor to a vector
      ! for unSymmmetric tensor you can use "reshape" function

      USE PARAMETERS

      IMPLICIT NONE

      integer:: i
      real*8 :: ATens(3,3), AVect(nSymm,1)

      do i = 1, 3
        AVect(i,1) = ATens(i,i)
      enddo

      AVect(4,1) = ATens(2,3)
      AVect(5,1) = ATens(1,3)
      AVect(6,1) = ATens(1,2)

      RETURN
      END SUBROUTINE symtensor2vector3

************************************************************************

      SUBROUTINE vector2symtensor2(Avect,Atens)
      ! this subroutine transforms a 3x1 Voigt vector to 2x2 symmetric tensor

      USE PARAMETERS

      IMPLICIT NONE

      integer:: i
      real*8 :: ATens(2,2), AVect(3,1)

      do i = 1, 2
        ATens(i,i) = AVect(i,1)
      enddo

      ATens(1,2) = AVect(3,1)
      ATens(2,1) = ATens(1,2)


      RETURN
      END SUBROUTINE vector2symtensor2

************************************************************************
      SUBROUTINE vector2symtensor3(Avect,Atens)
      ! this subroutine transforms a 6x1 Voigt vector to 3x3 symmetric tensor

      USE PARAMETERS

      IMPLICIT NONE

      integer:: i
      real*8 :: AVect(6,1), ATens(3,3)

      do i = 1, 3
        ATens(i,i) = AVect(i,1)
      enddo

      ATens(2,3) = AVect(4,1)
      ATens(1,3) = AVect(5,1)
      ATens(1,2) = AVect(6,1)
      ATens(2,1) = ATens(1,2)
      ATens(3,1) = ATens(1,3)
      ATens(3,2) = ATens(2,3)

      RETURN
      END SUBROUTINE vector2symtensor3

************************************************************************
      SUBROUTINE tangent2matrix(C,D)

      ! this subroutine maps the fourth order material/spatial tangent
      ! tensor (3x3x3x3) to a 2nd order stiffness tensor (6x6) using
      ! voigt notation: 11> 1, 22> 2, 33> 3, 23/32> 4, 13/31> 5, 12/21> 6

      USE PARAMETERS

      IMPLICIT NONE

      integer:: i, j, k, l, rw, cl
      integer:: Voigt(nSymm,2)
      real*8 :: C(3,3,3,3), D(nSymm,nSymm)

      ! Voigt convetion: (1,1) (2,2) (3,3) (2,3) (1,3) (1,2)
      Voigt = reshape((/ 1, 2, 3, 2, 1, 1,  1, 2, 3, 3, 3, 2 /),
     &        shape(Voigt))

      do rw = 1, nSymm
        do cl = 1, nSymm
          i = Voigt(rw,1)
          j = Voigt(rw,2)
          k = Voigt(cl,1)
          l = Voigt(cl,2)

          D(rw,cl) = C(i,j,k,l)
        enddo
      enddo

      RETURN
      END SUBROUTINE tangent2matrix

************************************************************************
************************************************************************

!     ADDITIONAL SUBROUTINES AVAILABEL THROUGH ABAQUS
!     CONSULT ABAQUS MANUAL FOR THE DETAILS

!     CALL SINV(STRESS,SINV1,SINV2,NDI,NSHR)
!     CALL SPRINC(S,PS,LSTR,NDI,NSHR)
!     CALL SPRIND(S,PS,AN,LSTR,NDI,NSHR)
!     CALL ROTSIG(S,R,SPRIME,LSTR,NDI,NSHR)

************************************************************************
************************************************************************