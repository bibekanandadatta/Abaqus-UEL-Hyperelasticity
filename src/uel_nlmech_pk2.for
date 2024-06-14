! **********************************************************************
! ********* ABAQUS/ STANDARD USER ELEMENT SUBROUTINE (UEL) *************
! **********************************************************************
!  large strain displacement element + Neo-Hookean & Arruda-Boyce model 
! **********************************************************************
!                     BIBEKANANDA DATTA (C) MAY 2024
!                 JOHNS HOPKINS UNIVERSITY, BALTIMORE, MD
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
!     Kappa       = props(2)        Bulk modulus
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
! **********************************************************************

      module user_element

      ! This module contains subroutines related to element formulation 
      ! and constitutive calculation. Abaqus user subroutines can not 
      ! be included in a module. Instead we extended the list of arguments
      ! of the Abaqus UEL subroutine and wrote another subroutine of
      ! similar kind which is included in the user_element module.
      ! Compilers can perform additional checks on the arguments when 
      ! any modularized subroutines are called. The first subroutine is
      ! called by UEL subroutine of Abaqus with an extended set of
      ! input arguments. The first subroutine calls other subroutines.

      contains

! **********************************************************************
! **********************************************************************

      subroutine uelNLMECH(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL)

      ! This subroutine contains the standard displacement-based 
      ! element formulation for static/ quasi-static large deformation
      ! of solids in total Lagrangian framework. This uses PK-2 stress 
      ! tangent formulation instead of PK-I and Cauchy stress. The later
      ! is known as updated Lagrangian framework. It calls the material
      ! model subroutine at each integration point to obtain the stress
      ! vector and tangent matrix used in the formulation. Currently
      ! available elements are 2D and 3D continuum elements of different
      ! shapes (TRI, QUAD, TET, HEX) and polynomial order (linear and
      ! qaudratic) with full and reduced integration. No specialzed
      ! numerical technique was employed to alleviate volumetric locking.

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

      ! input arguments to the subroutine
      integer, intent(in)   :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD
      integer, intent(in)   :: NNODE, JTYPE, KSTEP, KINC, JELEM
      integer, intent(in)   :: NDLOAD, JDLTYP, NPREDF, LFLAGS
      integer, intent(in)   :: MLVARX, MDLOAD, JPROPS, NJPROPS

      real(wp), intent(in)  :: PROPS, COORDS, DUall, Uall, Vel, Accn
      real(wp), intent(in)  :: TIME, DTIME, PARAMS, ADLMAG, PREDEF
      real(wp), intent(in)  :: DDLMAG, PERIOD

      character(len=2), intent(in)    :: analysis
      integer, intent(in)             :: nDim, nStress
      integer, intent(in)             :: uDOF, uDOFEL, nInt

      ! output of the suboutine
      real(wp), intent(out)           :: RHS, AMATRX
      real(wp), intent(out), optional :: SVARS, ENERGY, PNEWDT

      ! variables local to the subroutine
      real(wp)              :: ID(nDim,nDim)

      ! reshaped nodal degrees of freedom 
      real(wp)              :: uNode(nDim,nNode), duNode(nDim,nNode)

      ! Gauss quadrature weights and coordinates
      real(wp)              :: w(nInt), xi(nInt,nDim)
      
      ! interpolation function and their derivatives
      real(wp)              :: Nxi(nNode), dNdxi(nNode,nDim)

      ! element operator matrices
      real(wp)              :: dXdxi(nDim,nDim), dxidX(nDim,nDim)
      real(wp)              :: dNdX(nNode,nDim), detJ
      real(wp)              :: Na(nDim,nDim), Nmat(nDim,uDOFEl)
      real(wp)              :: Ba(nStress,nDim), Bmat(nStress,uDOFEl)
      real(wp)              :: Ga(nDim**2,nDim), Gmat(nDim**2,uDOFEl)
      real(wp)              :: stressTensorPK2(nDim,nDim)
      real(wp)              :: SIGMA_S(nDim**2,nDim**2)
      real(wp)              :: SIGMA_F(nDim*nNode,nDim*nNode)

      ! material point quantities (variables)
      real(wp)              :: F(3,3)
      real(wp)              :: strainLagrange(nStress,1)
      real(wp)              :: strainEuler(nStress,1)
      real(wp)              :: stressCauchy(nStress,1)
      real(wp)              :: stressPK1(nDim*nDim,1)
      real(wp)              :: stressPK2(nStress,1)
      real(wp)              :: Dmat(nStress,nStress)

      ! additional field variables (at nodes and int pt)
      real(wp)              :: fieldNode(npredf,nNode)
      real(wp)              :: dfieldNode(npredf,nNode)
      real(wp)              :: fieldVar(npredf), dfieldVar(npredf)

      ! element stiffness matrix and residual vector
      real(wp)              :: Kuu(uDOFEl,uDOFEl), Ru(uDOFEl,1)

      ! loop counter variables
      integer               :: i, j, intPt
      integer               :: matID

      ! element type
      type(element)         :: solidFiniteStrain
      ! error message logger
      type(logger)          :: msg

      ! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT

      ! set the element parameters
      solidFiniteStrain = element(nDim=nDim,analysis=analysis,
     &                            nNode=nNode,nInt=nInt)

      ! initialize the matrices and vectors
      F       = zero
      Na      = zero
      Ba      = zero
      Ga      = zero
      Nmat    = zero
      Bmat    = zero
      Gmat    = zero
      SIGMA_F = zero
      SIGMA_S = zero
      Kuu     = zero
      Ru      = zero

      ENERGY                    = zero
      AMATRX(1:NDOFEL,1:NDOFEL) = zero
      RHS(1:MLVARX,1)           = zero

      matID                     = jprops(2)

     !!!!!!!!!!! END VARIABLE DECLARATION AND INITIALIZATION !!!!!!!!!!!

      ! reshape the displacement vectors into matrix forms
      uNode  = reshape(UAll,[nDim,nNode])
      duNode = reshape(DUAll(:,1),[nDim,nNode])

      ! ! if applicable gather the prescribed field variables in a matrix
      ! ! such as temperature/ something (as shown below - not yet tested)
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

        call calcInterpFunc(solidFiniteStrain, xi(intPt,:), Nxi, dNdxi)

       ! calculate element jacobian and global shape func gradient
        dXdxi = matmul(coords,dNdxi)        ! calculate the jacobian (dXdxi)
        detJ  = det(dXdxi)                  ! calculate jacobian determinant

        if (detJ .le. zero) then
          call msg%ferror( flag=warn, src='uelNLMech',
     &     msg='Negative element jacobian.', ivec=[jelem, intpt])
        end if


        dxidX = inv(dXdxi)                  ! calculate jacobian inverse
        dNdX  = matmul(dNdxi,dxidX)         ! calculate dNdX


        ! loop over all the nodes (internal loop)
        do i=1,nNode

          ! form the nodal-level matrices: [Na], [Ga], [Ba]
          do j = 1, nDim
            Na(j,j) = Nxi(i)
            Ga(nDim*(j-1)+1:nDim*j,1:nDim) = dNdX(i,j)*ID
          end do

          ! form [Ba] matrix: plane stress/ plane strain case
          if (analysis .eq. 'PE') then
            Ba(1,1)       = dNdx(i,1)
            Ba(2,2)       = dNdx(i,2)
            Ba(3,1:nDim)  = [dNdx(i,2), dNdx(i,1)]

          ! form [Ba] matrix: 3D case
          else if (analysis .eq. '3D') then
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

        if (analysis .eq. 'PE') F(3,3) = one

    !     ! interpolate the field variables at the integration point
    !     ! (this not yet tested or used)
    !     do k = 1, npredf
    !       fieldVar(k)   = dot_product( Nxi, 
    !  &                    reshape( fieldNode(k,1:nNode), [nNode] ) )
    !       dfieldVar(k)  = dot_product( Nxi, 
    !  &                    reshape( dfieldNode(k,1:nNode), [nNode] ) )
    !     end do

        ! call material point subroutine (UMAT) for specific material
        if (matID .eq. 1) then
          call umatNeoHookean(kstep,kinc,time,dtime,nDim,analysis,
     &            nstress,nNode,jelem,coords,intpt,props,nprops,
     &            jprops,njprops,F,svars,nsvars,fieldVar,dfieldVar,
     &            npredf,
     &            stressPK2,Dmat,
     &            stressCauchy,stressPK1,strainLagrange,strainEuler)

        else if (matID .eq. 2) then
          call umatArrudaBoyce(kstep,kinc,time,dtime,nDim,analysis,
     &            nstress,nNode,jelem,coords,intpt,props,nprops,
     &            jprops,njprops,F,svars,nsvars,fieldVar,dfieldVar,
     &            npredf,
     &            stressPK2,Dmat,
     &            stressCauchy,stressPK1,strainLagrange,strainEuler)
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
          enddo
        enddo

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
        Kuu = Kuu + w(intpt) * detJ *
     &      ( matmul( transpose(matmul(Bmat,SIGMA_F)),
     &        matmul (Dmat, matmul(Bmat,SIGMA_F)) ) +
     &      matmul( transpose(Gmat), matmul(SIGMA_S,Gmat)) )

        Ru  = Ru - w(intpt) * detJ *
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

      subroutine umatNeoHookean(kstep,kinc,time,dtime,nDim,analysis,
     &            nstress,nNode,jelem,coords,intpt,props,nprops,
     &            jprops,njprops,F,svars,nsvars,fieldVar,dfieldVar,
     &            npredf,
     &            stressPK2,Dmat,
     &            stressCauchy,stressPK1,strainLagrange,strainEuler)

      ! This material point subroutine calculates constitutive response
      ! of a Neo-Hookean type material and returns the PK-II stress and 
      ! the material tangent as outputs. Optionally, it can return some 
      ! other strain and stress quantities in vector form (needs to be 
      ! programmed). All the constitutive calculations are initially done
      ! in 3D and later the corresponding matrices are reshaped based on
      ! the type of analysis is being performed by the user.
      ! This material subroutine also stores the user-defined element
      ! output in a global array for post=processing in Abaqus/Viewer.

      use global_parameters
      use error_logging
      use linear_algebra
      use lagrange_element
      use solid_mechanics
      use post_processing

      implicit none

      ! input arguments to the subroutine
      character(len=2), intent(in)  :: analysis

      integer, intent(in)   :: kstep, kinc, nDim, nstress
      integer, intent(in)   :: nNode, jelem, intpt, nprops
      integer, intent(in)   :: njprops, nsvars, npredf

      real(wp), intent(in)  :: time(2), dtime
      real(wp), intent(in)  :: coords(nDim,nNode)
      real(wp), intent(in)  :: props(nprops)
      integer,  intent(in)  :: jprops(njprops)

      real(wp), intent(in)  :: F(3,3)
      real(wp), intent(in)  :: fieldVar(npredf)
      real(wp), intent(in)  :: dfieldVar(npredf)

      ! output from the subroutine
      real(wp), intent(out) :: stressPK2(nStress,1)
      real(wp), intent(out) :: Dmat(nStress,nStress)

      ! optional output from the subroutine
      real(wp), intent(inout), optional :: svars(nsvars)
      real(wp), intent(out), optional   :: stressCauchy(nStress,1)
      real(wp), intent(out), optional   :: stressPK1(nDim*nDim,1)
      real(wp), intent(out), optional   :: strainLagrange(nstress,1)
      real(wp), intent(out), optional   :: strainEuler(nstress,1)  
      
      ! variables local to the subroutine
      ! material properties
      real(wp)              :: Gshear, Kappa, lam_L

      ! calculated kinematic quantities (variables)
      real(wp)              :: detF
      real(wp)              :: C(3,3), Cinv(3,3), trC, detC
      real(wp)              :: B(3,3), Binv(3,3), detB
      real(wp)              :: strainTensorEuler(3,3)
      real(wp)              :: strainTensorLagrange(3,3)

      ! stress tensors and material tangent
      real(wp)              :: stressTensorPK1(3,3)
      real(wp)              :: stressTensorPK2(3,3)
      real(wp)              :: stressTensorCauchy(3,3)
      real(wp)              :: Cmat(3,3,3,3)
      real(wp)              :: strainVoigtEuler(nSymm,1)
      real(wp)              :: stressVoigtPK1(nUnsymmm,1)
      real(wp)              :: stressVoigtPK2(nSymm,1)
      real(wp)              :: stressVoigtCauchy(nSymm,1)
      real(wp)              :: VoigtMat(nSymm,nSymm)

      integer               :: i, j, k, l             ! loop counters
      type(logger)          :: msg

      ! initialize matrial stiffness tensors
      Cmat      = zero
      Dmat      = zero

      ! assign material properties to variables
      Gshear    = props(1)        ! Shear modulus
      Kappa     = props(2)        ! Bulk modulus
      lam_L     = props(3)        ! Locking stretch


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

      B     = matmul(F,transpose(F))
      C     = matmul(transpose(F),F)

      Binv  = inv(B)
      Cinv  = inv(C)

      trC   = trace(C)

      ! calculate strain tensors
      strainTensorEuler = half*(ID3-Binv)

      ! calculate material tangent, C_ijkl
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              Cmat(i,j,k,l) = Cmat(i,j,k,l)
     &            + Kappa * Cinv(i,j) * Cinv(k,l)
     &            + ( Gshear-Kappa*log(detF) ) 
     &            * ( Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k) )
            end do
          end do
        end do
      end do

      ! calculate stress tensors
      stressTensorPK2     = Gshear*(ID3-Cinv) + Kappa*log(detF)*Cinv

      stressTensorCauchy  = (one/detF)*( Gshear*(B-ID3) +
     & 											Kappa*log(detF)*ID3 )


      ! transforms the stiffness tensor 3x3x3x3 to a 6x6 matrix
      call symtangent2matrix(Cmat,VoigtMat)

      ! transform the stress tensor (3x3) to Voigt vector form (6x1)
      call symtensor2vector(strainTensorEuler, strainVoigtEuler)
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

      subroutine umatArrudaBoyce(kstep,kinc,time,dtime,nDim,analysis,
     &            nstress,nNode,jelem,coords,intpt,props,nprops,
     &            jprops,njprops,F,svars,nsvars,fieldVar,dfieldVar,
     &            npredf,
     &            stressPK2,Dmat,
     &            stressCauchy,stressPK1,strainLagrange,strainEuler)

      ! This material point subroutine calculates constitutive response
      ! of a Arruda-Boyce type material and returns the PK-II stress and 
      ! the material tangent as outputs. Optionally, it can return some 
      ! other strain and stress quantities in vector form (needs to be 
      ! programmed). All the constitutive calculations are initially done
      ! in 3D and later the corresponding matrices are reshaped based on
      ! the type of analysis is being performed by the user.
      ! This material subroutine also stores the user-defined element
      ! output in a global array for post=processing in Abaqus/Viewer.

      use global_parameters
      use linear_algebra
      use lagrange_element
      use solid_mechanics
      use post_processing
      use error_logging

      implicit none

      ! input arguments to the subroutine
      character(len=2), intent(in)  :: analysis

      integer, intent(in)   :: kstep, kinc, nDim, nstress
      integer, intent(in)   :: nNode, jelem, intpt, nprops
      integer, intent(in)   :: njprops, nsvars, npredf

      real(wp), intent(in)  :: time(2), dtime
      real(wp), intent(in)  :: coords(nDim,nNode)
      real(wp), intent(in)  :: props(nprops)
      integer,  intent(in)  :: jprops(njprops)

      real(wp), intent(in)  :: F(3,3)
      real(wp), intent(in)  :: fieldVar(npredf)
      real(wp), intent(in)  :: dfieldVar(npredf)

      ! output from the subroutine
      real(wp), intent(out) :: stressPK2(nStress,1)
      real(wp), intent(out) :: Dmat(nStress,nStress)

      ! optional output from the subroutine
      real(wp), intent(inout), optional :: svars(nsvars)
      real(wp), intent(out), optional   :: stressCauchy(nStress,1)
      real(wp), intent(out), optional   :: stressPK1(nDim*nDim,1)
      real(wp), intent(out), optional   :: strainLagrange(nstress,1)
      real(wp), intent(out), optional   :: strainEuler(nstress,1)      
      
      ! variables local to the subroutine
      ! material properties
      real(wp)              :: Gshear, Kappa, lam_L
      real(wp)              :: lam_c, lam_r, beta_c, dBeta_c

      ! calculated kinematic quantities (variables)
      real(wp)              :: detF
      real(wp)              :: C(3,3), Cinv(3,3), trC, detC
      real(wp)              :: B(3,3), Binv(3,3), detB
      real(wp)              :: strainTensorEuler(3,3)
      real(wp)              :: strainTensorLagrange(3,3)

      ! stress tensors and material tangent
      real(wp)              :: stressTensorPK1(3,3)
      real(wp)              :: stressTensorPK2(3,3)
      real(wp)              :: stressTensorCauchy(3,3)
      real(wp)              :: Cmat(3,3,3,3)
      real(wp)              :: strainVoigtEuler(nSymm,1)
      real(wp)              :: stressVoigtPK1(nUnsymmm,1)
      real(wp)              :: stressVoigtPK2(nSymm,1)
      real(wp)              :: stressVoigtCauchy(nSymm,1)
      real(wp)              :: VoigtMat(nSymm,nSymm)
     
      integer               :: i, j, k, l
      type(logger)          :: msg

      ! initialize matrial stiffness tensors
      Cmat    = zero
      Dmat    = zero

      ! assign material properties to variables
      Gshear  = props(1)        ! Shear modulus
      Kappa   = props(2)        ! Bulk modulus
      lam_L   = props(3) 				! Locking stretch


      if (lam_L .le. zero) then
        call msg%ferror(flag=error, src='umatArrudaBoyce',
     &       msg='Incorrect material parameter (lam_L).', ra=lam_l)
        call xit
      end if

      ! perform all the constitutitve relations in 3D
      call detMat(F,detF)

      if (detF .le. zero) then
        call msg%ferror(flag = error, src = 'umatNeoHookean',
     &          msg = 'Issue with volume change (detF)', 
     &          ivec = [jelem, intpt], ra = detF)
        call xit
      end if

      B     = matmul(F,transpose(F))
      C     = matmul(transpose(F),F)

      Binv  = inv(B)
      Cinv  = inv(C)

      trC   = trace(C)

      ! calculate Euler-Almansi strain tensor
      strainTensorEuler = half*(ID3-Binv)

      lam_c     = sqrt(trC/three)
      lam_r     = lam_c/lam_L
      beta_c    = InvLangevin(lam_r)
      dBeta_c   = DInvLangevin(lam_r)

      ! form material tangent, C_ijkl
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              Cmat(i,j,k,l) = Cmat(i,j,k,l) 
     &            + Gshear/(nine*lam_c**two)
     &            * ( dBeta_c- lam_r*beta_c ) * ID3(i,j)*ID3(k,l)
     &            + Kappa * Cinv(i,j)*Cinv(k,l)
     &            + ( (Gshear/three) * lam_r*beta_c - Kappa*log(detF) )
     &            * ( Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k) )
            end do
          end do
        end do
      end do

      stressTensorCauchy  = (1/detF) * ( (Gshear/three)*lam_r*beta_c*B
     &      - ( (Gshear*lam_L)/three - Kappa*log(detF) ) * ID3 )

      stressTensorPK2     = (Gshear/three) * lam_r * beta_c * ID3 -
     &      ( Gshear*lam_L/three - Kappa*log(detF) ) * Cinv


      ! transforms the stiffness tensor 3x3x3x3 to a 6x6 matrix
      call symtangent2matrix(Cmat,VoigtMat)

      ! transform the stress tensor (3x3) to Voigt vector form (6x1)
      call symtensor2vector(strainTensorEuler, strainVoigtEuler)
      call symtensor2vector(stressTensorCauchy, stressVoigtCauchy)
      call symtensor2vector(stressTensorCauchy, stressVoigtPK2)

      !!!!!!!!!!!!!!!! END OF CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!


      ! reshape the Voigt matrix and tensor based on the analysis
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
        strainEuler     = strainVoigtEuler
        stressPK2       = stressVoigtPK2
        stressCauchy    = stressVoigtCauchy

      else
        call msg%ferror(flag=error, src='umatArrudaBoyce',
     &            msg='Wrong analysis.', ch=analysis)
        call xit
      end if

      ! save the variables to be post-processed in globalPostVars
      globalPostVars(jelem,intpt,1:nStress) = stressCauchy(1:nStress,1)
      globalPostVars(jelem,intpt,nStress+1:2*nStress) =
     & 								  strainEuler(1:nStress,1)

! **********************************************************************

      contains

      function InvLangevin(x)

      ! calculates an approximation of the inverse Langevin function
      ! reference: Bergstorm (PhD thesis, MIT, 1999)

      implicit none

      real(wp),intent(in)  :: x
      real(wp)             :: InvLangevin

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

! **********************************************************************

      function DInvLangevin(x)

      ! calculates an approximation of derivative of
      ! the inverse Langevin function
      ! reference: Bergstorm (PhD thesis, MIT, 1999)

      implicit none

      real(wp), intent(in)   :: x
      real(wp)               :: DInvLangevin, sec

      if (abs(x) .lt. 0.84136_wp) then
        DInvLangevin = 2.0898073756_wp*(tan(1.58986_wp*x))**two
     &                + 3.0018973756_wp
      else if ((abs(x) .ge. 0.84136_wp) .and. (abs(x) .lt. one)) then
        DInvLangevin = one/( (sign(one,x)-x)**two )
      else
       call msg%ferror(flag=error, src='umatArrudaBoyce:DInvLangevin',
     &                  msg='Unbound argument.', ra = x)
        call xit
      end if

      return

      end function DInvLangevin

      end subroutine umatArrudaBoyce

      end module user_element

! **********************************************************************
! ****************** ABAQUS USER ELEMENT SUBROUTINE ********************
! **********************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)

      ! This subroutine is called by Abaqus with above arguments
      ! for each user elements defined in an Abaqus model. Users are
      ! responsible for programming the element tangent/ stiffness
      ! matrix and residual vectors which will be then assembled and
      ! solved by Abaqus after applying the boundary conditions.
      
      use global_parameters
      use error_logging
      use user_element
      use post_processing

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     & DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      integer, intent(in)   :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD
      integer, intent(in)   :: NNODE, JTYPE, KSTEP, KINC, JELEM
      integer, intent(in)   :: NDLOAD, JDLTYP, NPREDF, LFLAGS
      integer, intent(in)   :: MLVARX, MDLOAD, JPROPS, NJPROPS

      real(wp), intent(in)  :: PROPS, COORDS, DUall, Uall, Vel, Accn
      real(wp), intent(in)  :: TIME, DTIME, PARAMS, ADLMAG, PREDEF
      real(wp), intent(in)  :: DDLMAG, PERIOD

      real(wp), intent(out)           :: RHS, AMATRX
      real(wp), intent(out), optional :: SVARS, ENERGY, PNEWDT


      integer             :: nInt, matID, nPostVars
      integer             :: nDim, nStress, uDOF, uDOFEL
      character(len=2)    :: analysis
      character(len=8)    :: abqProcedure
      logical             :: nlgeom
     

      integer             :: lenJobName,lenOutDir
      character(len=256)  :: outDir
      character(len=256)  :: jobName
      character(len=512)  :: errFile, dbgFile
      type(logger)        :: msg


      ! ! open a debug file for the current job
      ! call getJobName(jobName,lenJobName)
      ! call getOutDir(outDir,lenOutDir)
      ! errFile = trim(outDir)//'\aaERR_'//trim(jobName)//'.dat'
      ! dbgFile = trim(outDir)//'\aaDBG_'//trim(jobName)//'.dat'
      ! call msg%fopen( errfile=errFile, dbgfile=dbgFile )


      ! change the LFLAGS criteria as needed (check abaqus UEL manual)
      if((lflags(1) .eq. 1).or.(lflags(1) .eq. 2)) then
        abqProcedure = 'STATIC'
      else
       call msg%ferror(flag=error, src='UEL',
     &                msg='Incorrect Abaqus procedure.', ia=lflags(1))
        call xit
      end if

      ! check if the procedure is linear or nonlinear
      if (lflags(2) .eq. 0) then
        nlgeom = .false.
      else if (lflags(2) .eq. 1) then
        nlgeom = .true.
      end if

      ! check to see if it's a general step or a linear purturbation step
      if(lflags(4) .eq. 1) then
         call msg%ferror(flag=error, src='UEL',
     &      msg='The step should be a GENERAL step.', ia=lflags(4))
        call xit
      end if


       ! set parameters specific to analysis and element types
      if ((jtype .ge. 1).and.(jtype .le. 4)) then
        nDim      = 3
        analysis  = '3D'            ! three-dimensional analysis
        nStress   = 6
        uDOF      = nDim            ! displacement degrees of freedom of a node
        uDOFEL    = nNode*uDOF      ! total displacement degrees of freedom in element
      else if ((jtype .ge. 5).and.(jtype .le. 8)) then
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
        call msg%finfo('------- Abaqus FINITE STRAIN UEL ------')
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

       ! call the element subroutine with extended input arguments
       call uelNLMECH(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL)

      END SUBROUTINE UEL

! **********************************************************************
! ************** ABAQUS USER OUTPUT VARIABLES SUBROUTINE ***************
! **********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     & NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     & JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is called by Abaqus at each material point (int pt)
      ! to obtain the user defined output variables for standard Abaqus 
      ! elements. We used an additional layer of standard Abaqus elements
      ! with same topology (same number of nodes and int pts) on top of 
      ! the user elements with an offset in the numbering between the user 
      ! elements and standard elements. This number is defined in the
      ! post_processing module and should match with Abaqus input file.

      use global_parameters
      use post_processing

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      ! the dimensions of the variables FLGRAY, ARRAY and JARRAY
      ! must be set equal to or greater than 15.

      ! explicityly define the type for uvar to avoid issues
      real(wp)        :: uvar

      ! assign the stored global variables to the UVAR for Abaqus to process
      uvar(1:nuvarm)  = globalPostVars(noel-elemOffset,npt,1:nuvarm)


      END SUBROUTINE UVARM

! **********************************************************************
! **********************************************************************