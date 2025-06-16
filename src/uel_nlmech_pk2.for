! **********************************************************************
! ********* ABAQUS/ STANDARD USER ELEMENT SUBROUTINE (UEL) *************
! **********************************************************************
!  large strain displacement element + Neo-Hookean & Arruda-Boyce model
!   linear quad and hex element formulation use F-bar method to avoid
!   volumetric locking at near-incompressibility limit (de Souza Neto)
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
! ***********************************************************************
!          VOIGT NOTATION CONVENTION FOR STRESS/ STRAIN TENSORS
!
!       In this subroutine we adopted the following convention for
!       symmetric stress and strain tensor following Voigt notation
!       This is different than what is followed by Abaqus/ Standard
!
!          sigma11, sigma22, sigma33, sigma23, sigma13, sigma12
!       strain11, strain22, strain33, strain23, strain13, strain12
! **********************************************************************
!                        LIST OF ELEMENT PROPERTIES
!
!     jprops(1)   = nInt            no of integration points in element
!     jprops(2)   = fbarFlag        flag to use F-bar formulation
!     jprops(3)   = matID           constitutive relation for material
!     jprops(4)   = nPostVars       no of local (int pt) post-processing variables
! **********************************************************************
!                        POST-PROCESSED VARIABLES
!                     (follows the convention above)
!
!     uvar(1:nStress)               Cauchy stress tensor components
!     uvar(nStress+1:2*nStress)     Euler strain tensor components
! **********************************************************************
!
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
!     MDLOAD                        Total number of distributed loads and/or fluxes defined on this element
!     PERIOD                        Time period of the current step
!
! **********************************************************************

      ! make sure to have the correct directory
      include 'global_parameters.for'     ! parameter module
      include 'error_logging.for'         ! error/ debugging module
      include 'linear_algebra.for'        ! linear algebra module
      include 'nonlinear_solver.for'      ! Newton-Raphson solver module
      include 'lagrange_element.for'      ! module for Lagrange elements
      include 'gauss_quadrature.for'      ! Guassian quadrature module
      include 'solid_mechanics.for'       ! solid mechanics module
      include 'post_processing.for'       ! post-processing module

! **********************************************************************
! **********************************************************************

      module hyperelastic_material

      contains

      subroutine mat_NeoHookean(kstep,kinc,time,dtime,nDim,analysis,
     &            nstress,nNode,jelem,coords,intpt,props,nprops,
     &            jprops,njprops,F,svars,nsvars,fieldVar,dfieldVar,
     &            npredf,stressPK2,Dmat,Cmat)

! **********************************************************************
!     This material point subroutine calculates constitutive response
!     of a Neo-Hookean type material and returns the PK-II stress and
!     the material tangent as outputs. Optionally, it can return some
!     other strain and stress quantities in vector form (needs to be
!     programmed). All the constitutive calculations are initially done
!     in 3D and later the corresponding matrices are reshaped based on
!     the type of analysis is being performed by the user.
!     This material subroutine also stores the user-defined element
!     output in a global array for post=processing in Abaqus/Viewer.
!
!                       LIST OF MATERIAL PROPERTIES
!
!     G           = props(1)        Shear modulus
!     Kappa       = props(2)        Bulk modulus
!     lambda_L    = props(3)        Locking stretch for AB model
!                                   This is zero for the NH model
! **********************************************************************

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

      ! output arguments from the subroutine
      real(wp), intent(out) :: stressPK2(nStress,1)
      real(wp), intent(out) :: Dmat(nStress,nStress)
      real(wp), intent(out) :: Cmat(3,3,3,3)

      ! optional output from the subroutine
      real(wp), intent(inout), optional :: svars(nsvars)


      ! material properties
      real(wp)              :: Gshear, Kappa, lam_L

      ! calculated kinematic quantities (variables)
      real(wp)              :: detF
      real(wp)              :: C(3,3), Cinv(3,3), trC, detC
      real(wp)              :: B(3,3), Binv(3,3), detB

      ! intermediate variables for stress tensors and material tangent output
      real(wp)              :: stressTensorPK2(3,3)
      real(wp)              :: stressVectPK2(nSymm,1)
      real(wp)              :: VoigtMat(nSymm,nSymm)

      ! output variables (3x3 stress and strain tensors)
      real(wp)              :: strainTensorLagrange(3,3)
      real(wp)              :: strainTensorEuler(3,3)
      real(wp)              :: stressTensorPK1(3,3)
      real(wp)              :: stressTensorCauchy(3,3)

      ! vector form (6x1 or 9x1) of stress and strain tensors
      real(wp)              :: strainVectLagrange(nSymm,1)
      real(wp)              :: strainVectEuler(nSymm,1)
      real(wp)              :: stressVectPK1(nUnSymm,1)
      real(wp)              :: stressVectCauchy(nSymm,1)

      ! final vector form of stress and strain tensors based on analysis
      real(wp)              :: stressPK1(nDim*nDim,1)
      real(wp)              :: stressCauchy(nStress,1)
      real(wp)              :: strainLagrange(nStress,1)
      real(wp)              :: strainEuler(nStress,1)

      integer               :: i, j, k, l     ! loop counters
      type(logger)          :: msg            ! error message logger

      ! initialize matrial stiffness tensors
      Cmat      = zero
      Dmat      = zero

      ! assign material properties to variables
      Gshear    = props(1)        ! Shear modulus
      Kappa     = props(2)        ! Bulk modulus
      lam_L     = props(3)        ! Locking stretch


      ! locking stretch should be infinity (0 as input) for NH model
      if (lam_L .ne. zero) then
        call msg%ferror(flag=error, src='mat_Neohookean',
     &       msg='Incorrect material parameter (lam_L).', ra=lam_l)
        call xit
      end if


      ! perform all the constitutitve relations in 3D
      detF  = det(F)

      if (detF .le. zero) then
        call msg%ferror(flag=error, src='mat_NeoHookean',
     &          msg='Issue with volume change (detF)',
     &          ivec=[jelem, intpt], ra= detF)
        call xit
      end if

      C     = matmul(transpose(F),F)
      Cinv  = inv(C)
      trC   = trace(C)
      B     = matmul(F,transpose(F))
      Binv  = inv(B)


      ! calculate strain tensors
      strainTensorLagrange  = half*(C-ID3)
      strainTensorEuler     = half*(ID3-Binv)


      ! calculate the stress tensors
      stressTensorPK2     = Gshear*(ID3-Cinv) + Kappa*log(detF)*Cinv
      stressTensorCauchy  = (one/detF)
     &                    * ( Gshear*(B-ID3) + Kappa*log(detF)*ID3 )
      stressTensorPK1     = matmul(F,stressTensorPK2)


      ! calculate the material tangent: C_ijkl = dS_ij/dC_kl
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

      !!!!!!!!!!!!!!! END OF CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!


      ! reshape the stiffness tensor 3x3x3x3 to a 6x6 matrix
      call voigtMatrix(Cmat,VoigtMat)

      ! reshape the strain and stress tensors into vectors
      ! dimension: (SYMM 6x1) or (UNSYMM 9x1)
      call voigtVector(strainTensorLagrange, strainVectLagrange)
      call voigtVector(strainTensorEuler, strainVectEuler)
      call unsymmVector(stressTensorPK1,stressVectPK1)
      call voigtVector(stressTensorPK2, stressVectPK2)
      call voigtVector(stressTensorCauchy, stressVectCauchy)


      ! reshape the Voigt matrix and tensor based on the analysis
      if ((analysis .eq. '3D') .or. (analysis .eq. 'PE')
     &      .or. (analysis .eq. 'AX')) then
        call voigtVectorTruncate(stressVectPK2,stressPK2)
        call voigtMatrixTruncate(VoigtMat,Dmat)

      else
        ! plane stress is not available
        call msg%ferror(flag=error, src='mat_NeoHookean',
     &            msg='Wrong analysis.', ch=analysis)
        call xit
      end if

      ! additional variable for post-processing
      call voigtVectorTruncate(strainVectLagrange,strainLagrange)
      call voigtVectorTruncate(strainVectEuler,strainEuler)
      call voigtVectorTruncate(stressVectCauchy,stressCauchy)


      ! save the variables to be post-processed in globalPostVars
      globalPostVars(jelem,intpt,1:nStress) = stressCauchy(1:nStress,1)
      globalPostVars(jelem,intpt,nStress+1:2*nStress)
     &                                      = strainEuler(1:nStress,1)

      end subroutine mat_NeoHookean

! **********************************************************************
! **********************************************************************

      subroutine mat_ArrudaBoyce(kstep,kinc,time,dtime,nDim,analysis,
     &            nstress,nNode,jelem,coords,intpt,props,nprops,
     &            jprops,njprops,F,svars,nsvars,fieldVar,dfieldVar,
     &            npredf,stressPK2,Dmat,Cmat)

! **********************************************************************
!     This material point subroutine calculates constitutive response
!     of a Arruda-Boyce type material and returns the PK-II stress and
!     the material tangent as outputs. Optionally, it can return some
!     other strain and stress quantities in vector form (needs to be
!     programmed). All the constitutive calculations are initially done
!     in 3D and later the corresponding matrices are reshaped based on
!     the type of analysis is being performed by the user.
!     This material subroutine also stores the user-defined element
!     output in a global array for post=processing in Abaqus/Viewer.
!
!                       LIST OF MATERIAL PROPERTIES
!
!     G           = props(1)        Shear modulus
!     Kappa       = props(2)        Bulk modulus
!     lambda_L    = props(3)        Locking stretch for AB model
! **********************************************************************

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
      real(wp), intent(out) :: Cmat(3,3,3,3)

      ! optional output from the subroutine
      real(wp), intent(inout), optional :: svars(nsvars)


      ! material properties
      real(wp)              :: Gshear, Kappa, lam_L

      ! calculated kinematic quantities (variables)
      real(wp)              :: detF
      real(wp)              :: C(3,3), Cinv(3,3), trC, detC
      real(wp)              :: B(3,3), Binv(3,3), detB
      real(wp)              :: lam_c, lam_r, beta_c, dBeta_c

      ! intermediate variables for stress tensors and material tangent output
      real(wp)              :: stressTensorPK2(3,3)
      real(wp)              :: stressVectPK2(nSymm,1)
      real(wp)              :: VoigtMat(nSymm,nSymm)

      ! output variables (3x3 stress and strain tensors)
      real(wp)              :: strainTensorLagrange(3,3)
      real(wp)              :: strainTensorEuler(3,3)
      real(wp)              :: stressTensorPK1(3,3)
      real(wp)              :: stressTensorCauchy(3,3)

      ! vector form (6x1 or 9x1) of stress and strain tensors
      real(wp)              :: strainVectLagrange(nSymm,1)
      real(wp)              :: strainVectEuler(nSymm,1)
      real(wp)              :: stressVectPK1(nUnSymm,1)
      real(wp)              :: stressVectCauchy(nSymm,1)

      ! final vector form of stress and strain tensors based on analysis
      real(wp)              :: stressPK1(nDim*nDim,1)
      real(wp)              :: stressCauchy(nStress,1)
      real(wp)              :: strainLagrange(nStress,1)
      real(wp)              :: strainEuler(nStress,1)

      integer               :: i, j, k, l     ! loop counters
      type(logger)          :: msg            ! error message logger

      ! initialize matrial stiffness tensors
      Cmat    = zero
      Dmat    = zero

      ! assign material properties to variables
      Gshear  = props(1)        ! Shear modulus
      Kappa   = props(2)        ! Bulk modulus
      lam_L   = props(3) 				! Locking stretch


      if (lam_L .le. zero) then
        call msg%ferror(flag=error, src='mat_ArrudaBoyce',
     &       msg='Incorrect material parameter (lam_L).', ra=lam_l)
        call xit
      end if

      ! perform all the constitutitve relations in 3D
      detF  = det(F)

      if (detF .le. zero) then
        call msg%ferror(flag = error, src = 'mat_NeoHookean',
     &          msg = 'Issue with volume change (detF)',
     &          ivec = [jelem, intpt], ra = detF)
        call xit
      end if

      ! calculate the kinematic quantities
      C         = matmul(transpose(F),F)
      Cinv      = inv(C)
      trC       = trace(C)
      B         = matmul(F,transpose(F))
      Binv      = inv(B)
      lam_c     = sqrt(trC/three)
      lam_r     = lam_c/lam_L
      beta_c    = InvLangevin(lam_r)
      dBeta_c   = DInvLangevin(lam_r)



      ! calculate strain tensors
      strainTensorLagrange  = half*(C-ID3)
      strainTensorEuler     = half*(ID3-Binv)


      ! calculate stress tensors
      stressTensorPK2     = (Gshear/three) * lam_r * beta_c * ID3
     &      - ( (Gshear*lam_L)/three - Kappa*log(detF) ) * Cinv

      stressTensorCauchy  = (1/detF) * ( (Gshear/three)*lam_r*beta_c*B
     &      - ( (Gshear*lam_L)/three - Kappa*log(detF) ) * ID3 )

      stressTensorPK1     = matmul(F,stressTensorPK2)


      ! form material tangent, C_ijkl = dS_ij/dC_kl
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              Cmat(i,j,k,l) = Cmat(i,j,k,l)
     &            + Gshear/(nine*lam_c**two)
     &            * ( dBeta_c- lam_r*beta_c ) * ID3(i,j)*ID3(k,l)
     &            + Kappa * Cinv(i,j)*Cinv(k,l)
     &            + ( (Gshear*lam_L)/three - Kappa*log(detF) )
     &            * ( Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k) )
            end do
          end do
        end do
      end do

      !!!!!!!!!!!!!!!! END OF CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!


       ! reshape the stiffness tensor 3x3x3x3 to a 6x6 matrix
      call voigtMatrix(Cmat,VoigtMat)

      ! reshape the strain and stress tensors into vectors
      ! dimension: (SYMM 6x1) or (UNSYMM 9x1)
      call voigtVector(strainTensorLagrange, strainVectLagrange)
      call voigtVector(strainTensorEuler, strainVectEuler)
      call unsymmVector(stressTensorPK1,stressVectPK1)
      call voigtVector(stressTensorPK2, stressVectPK2)
      call voigtVector(stressTensorCauchy, stressVectCauchy)


      ! reshape the Voigt matrix and tensor based on the analysis
      if ((analysis .eq. '3D') .or. (analysis .eq. 'PE')
     &      .or. (analysis .eq. 'AX')) then
        call voigtVectorTruncate(stressVectPK2,stressPK2)
        call voigtMatrixTruncate(VoigtMat,Dmat)

      else
        ! plane stress is not available
        call msg%ferror(flag=error, src='mat_ArrudaBoyce',
     &            msg='Wrong analysis.', ch=analysis)
        call xit
      end if

      ! additional variable for post-processing
      call voigtVectorTruncate(strainVectLagrange,strainLagrange)
      call voigtVectorTruncate(strainVectEuler,strainEuler)
      call voigtVectorTruncate(stressVectPK2,stressPK2)
      call voigtVectorTruncate(stressVectCauchy,stressCauchy)


      ! save the variables to be post-processed in globalPostVars
      globalPostVars(jelem,intpt,1:nStress) = stressCauchy(1:nStress,1)
      globalPostVars(jelem,intpt,nStress+1:2*nStress)
     & 								                      = strainEuler(1:nStress,1)

! **********************************************************************

      contains

      function InvLangevin(x)

      ! calculates an approximation of the inverse Langevin function
      ! reference: Bergstorm (PhD thesis, MIT, 1999)

      implicit none

      real(wp), intent(in)  :: x
      real(wp)              :: InvLangevin

      if (abs(x) .lt. 0.84136_wp) then
        InvLangevin = 1.31446_wp*tan(1.58986_wp*x) + 0.91209_wp*x

      else if ((abs(x) .ge. 0.84136_wp) .and. (abs(x) .lt. one)) then
        InvLangevin = one/(sign(one,x)-x)

      else
        call msg%ferror(flag=error, src='mat_ArrudaBoyce:InvLangevin',
     &                  msg='Unbound argument.', ra = x)
        call xit
      end if

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
       call msg%ferror(flag=error, src='mat_ArrudaBoyce:DInvLangevin',
     &                  msg='Unbound argument.', ra = x)
        call xit
      end if

      end function DInvLangevin

      end subroutine mat_ArrudaBoyce


      end module hyperelastic_material

! **********************************************************************
! **********************************************************************
! **********************************************************************

      module finite_strain_element

! **********************************************************************
!     This module contains subroutines related to element formulation
!     and constitutive calculation. Abaqus user subroutines can not
!     be included in a module. Instead we extended the list of arguments
!     of the Abaqus UEL subroutine and wrote another subroutine of
!     similar kind which is included in the finite_strain_element module.
!     Compilers can perform additional checks on the arguments when
!     any modularized subroutines are called. The first subroutine is
!     called by UEL subroutine of Abaqus with an extended set of
!     input arguments. The first subroutine calls other subroutines.
! **********************************************************************

      contains

      subroutine elem_nlmech(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & NDIM,ANALYSIS,NSTRESS,NINT)

! **********************************************************************
!     This subroutine contains the standard displacement-based
!     element formulation for static/ quasi-static large deformation
!     of solids in total Lagrangian framework. This uses PK-I stress
!     tangent formulation instead of PK-II and Cauchy stress. The later
!     is known as updated Lagrangian framework. It calls the material
!     model subroutine at each integration point to obtain the stress
!     vector and tangent matrix used in the formulation. Currently
!     available elements are 2D and 3D continuum elements of different
!     shapes (TRI, QUAD, TET, HEX) and polynomial order (linear and
!     qaudratic) with full and reduced integration.
!     4-node linear QUAD and 8-node linear HEX elements use F-bar
!     technique proposed by de Souza Neto (IJSS, 1996) to alleviate
!     volumetric locking near incompressibility limit.
! **********************************************************************

      use global_parameters
      use error_logging
      use linear_algebra
      use lagrange_element
      use gauss_quadrature
      use solid_mechanics
      use post_processing
      use hyperelastic_material

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
      integer, intent(in)             :: nDim, nStress, nInt

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

      ! finite element matrix operators
      real(wp)              :: Na(nDim,nDim)
      real(wp)              :: Ba(nStress,nDim)
      real(wp)              :: Ga(nDim*nDim,nDim)
      real(wp)              :: Nmat(nDim,nDOFEl)
      real(wp)              :: NmatT(nDOFEL,nDim)
      real(wp)              :: Bmat(nStress,nDOFEl)
      real(wp)              :: BmatT(nDOFEL,nStress)
      real(wp)              :: Gmat(nDim*nDim,nDim*nNode)
      real(wp)              :: GmatT(nDim*nNode,nDim*nDim)

      ! additional reshaped matrices for element formulation
      real(wp)              :: stressTensorPK2(nDim,nDim)
      real(wp)              :: SIGMA_S(nDim**2,nDim**2)
      real(wp)              :: SIGMA_F(nDim*nNode,nDim*nNode)
      real(wp)              :: BNLmat(nStress,nDim*nNode)
      real(wp)              :: BNLmatT(nDim*nNode,nStress)

      ! material point quantities (variables)
      real(wp)              :: F(3,3), detF, Fbar(3,3)
      real(wp)              :: FInv(3,3), FInvT(3,3)
      real(wp)              :: strainLagrange(nStress,1)
      real(wp)              :: strainEuler(nStress,1)
      real(wp)              :: stressCauchy(nStress,1)
      real(wp)              :: stressPK1(nDim*nDim,1)
      real(wp)              :: stressPK2(nStress,1)
      real(wp)              :: Dmat(nStress,nStress)
      real(wp)              :: Cmat(3,3,3,3)


      ! additional variables for F-bar method (element and material)
      logical               :: fbarFlag
      real(wp)              :: centroid(nDim)
      real(wp)              :: Nxi0(nNode), dNdxi0(nNode,nDim)
      real(wp)              :: dXdxi0(nDim,nDim), dxidX0(nDim,nDim)
      real(wp)              :: dNdX0(nNode,nDim), detJ0
      real(wp)              :: Ga0(nDim**2,nDim), Gmat0(nDim**2,nDOFEL)
      real(wp)              :: F0(3,3), detF0
      real(wp)              :: F0Inv(3,3), F0InvT(3,3)
      real(wp)              :: QR0Tensor(nDim,nDim,nDim,nDim)
      real(wp)              :: QRTensor(nDim,nDim,nDim,nDim)
      real(wp)              :: QR0mat(nDim*nDim,nDim*nDim)
      real(wp)              :: QRmat(nDim*nDim,nDim*nDim)
      real(wp)              :: tanFac1, tanFac2, resFac


      ! additional field variables (at nodes and int pt)
      real(wp)              :: fieldNode(npredf,nNode)
      real(wp)              :: dfieldNode(npredf,nNode)
      real(wp)              :: fieldVar(npredf), dfieldVar(npredf)

      integer               :: i, j, k, l, m, n, p, q, intPt
      integer               :: matID                    ! constitutive mocel

      type(element)         :: solidFiniteStrain        ! element type
      type(logger)          :: msg                      ! error message logger

      ! element stiffness matrix and residual vector
      real(wp)              :: Kuu(nDOFEL,nDOFEL), Ru(nDOFEL,1)


      ! set the element parameters
      solidFiniteStrain = element(nDim=nDim,analysis=analysis,
     &                            nNode=nNode,nInt=nInt)

      ! initialize the matrices and vectors
      F0      = zero
      F       = zero
      Fbar    = zero
      Ga0     = zero
      Gmat0   = zero
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

      fbarFlag  = jprops(2)
      matID     = jprops(3)

      !!!!!!!!!! END VARIABLE DECLARATION AND INITIALIZATION !!!!!!!!!!!

      ! reshape the displacement vectors into matrix forms
      uNode  = reshape(UAll,[nDim,nNode])
      duNode = reshape(DUAll(:,1),[nDim,nNode])

      ! ! if applicable gather the prescribed field variables in a matrix
      ! ! such as temperature/ something (as shown below - not yet tested)
      ! do k = 1 , npredf
      !   fieldNode(k,1:nNode) = predef(1,k,1:nNode)
      !   dfieldNode(k,1:nNode) = predef(2,k,1:nNode)
      ! end do

      !!!!!!!!!!!!!!!!! ELEMENT RELATED OPERATIONS !!!!!!!!!!!!!!!!!!!!!

      call eyeMat(ID)

      ! For fully-integrated linear quad and hex elements, calculate Gmat0.
      ! These calculations are done to evaluate volumetric deformation
      ! gradient at centroid which will be used in to define F-bar later.
      if (fbarFlag .eq. .true.) then

        if ( ((jtype .eq. 3) .and. (nInt .eq. 8))
     &    .or. ((jtype .eq. 7) .and. (nInt .eq. 4)) ) then

          centroid = zero

          ! evaluate the interpolation functions and derivates at centroid
          call calcInterpFunc(solidFiniteStrain, centroid, Nxi0, dNdxi0)

          ! calculate element jacobian and global shape func gradient at centroid
          dXdxi0  = matmul(coords,dNdxi0)       ! calculate the jacobian (dXdxi) at centroid
          detJ0   = det(dXdxi0)                 ! calculate jacobian determinant at centroid

          if (detJ0 .le. zero) then
            call msg%ferror( flag=warn, src='elem_nlmech',
     &        msg='Negative element jacobian at centroid.', ia=jelem)
          end if

          dxidX0 = inv(dXdxi0)                  ! calculate jacobian inverse
          dNdX0  = matmul(dNdxi0,dxidX0)        ! calculate dNdX at centroid

          do i=1,nNode

            ! form the nodal-level matrix: [Ga0] at the centroid
            do j = 1, nDim
              Ga0(nDim*(j-1)+1:nDim*j,1:nDim) = dNdX0(i,j)*ID
            end do

            ! form the [G0] matrix at the centroid
            Gmat0(1:nDim**2,nDim*(i-1)+1:nDim*i) = Ga0(1:nDim**2,1:nDim)
          end do                             ! end of nodal point loop

          F0(1:nDim,1:nDim) = ID + matmul(uNode,dNdX0)

          if (analysis .eq. 'PE') F0(3,3) = one

          detF0   = det(F0)
          F0Inv   = inv(F0)
          F0InvT  = transpose(F0Inv)

        else
          call msg%ferror( flag=warn, src='elem_nlmech',
     &        msg='F-bar is not available: ', ivec=[jtype, nInt])
          call xit
        end if

      end if

      ! obtain gauss quadrature points and weights
      call getGaussQuadrtr(solidFiniteStrain,w,xi)

      ! loop through all the integration points (main/ external loop)
      do intPt = 1, nInt

        call calcInterpFunc(solidFiniteStrain, xi(intPt,:), Nxi, dNdxi)

       ! calculate element jacobian and global shape func gradient
        dXdxi = matmul(coords,dNdxi)        ! calculate the jacobian (dXdxi)
        detJ  = det(dXdxi)                  ! calculate jacobian determinant

        if (detJ .le. zero) then
          call msg%ferror( flag=warn, src='elem_nlmech',
     &         msg='Negative element jacobian.', ivec=[jelem, intpt])
        end if


        dxidX = inv(dXdxi)                  ! calculate jacobian inverse
        dNdX  = matmul(dNdxi,dxidX)         ! calculate dNdX


        ! loop over all the nodes (internal loop)
        do i = 1, nNode

          ! form the nodal-level matrices: [Na], [Ga], [Ba]
          do j = 1, nDim
            Na(j,j) = Nxi(i)
            Ga(nDim*(j-1)+1:nDim*j,1:nDim) = dNdX(i,j)*ID
          end do

          ! form [Ba] matrix: plane stress/ plane strain case
          if (analysis .eq. 'PE') then
            Ba(1,1)       = dNdX(i,1)
            Ba(2,2)       = dNdX(i,2)
            Ba(3,1:nDim)  = [dNdX(i,2), dNdX(i,1)]

          ! form [Ba] matrix: 3D case
          else if (analysis .eq. '3D') then
            Ba(1,1)       = dNdX(i,1)
            Ba(2,2)       = dNdX(i,2)
            Ba(3,3)       = dNdX(i,3)
            Ba(4,1:nDim)  = [  zero,      dNdX(i,3),  dNdX(i,2)]
            Ba(5,1:nDim)  = [dNdX(i,3),     zero,     dNdX(i,1)]
            Ba(6,1:nDim)  = [dNdX(i,2),   dNdX(i,1),    zero   ]

          else
            call msg%ferror( flag=error, src='elem_nlmech',
     &                  msg='Wrong analysis.', ch=analysis )
            call xit
          end if

          ! form the [N], [B], and [G] matrix
          Nmat(1:nDim,nDim*(i-1)+1:nDim*i)    = Na(1:nDim,1:nDim)
          Bmat(1:nStress,nDim*(i-1)+1:nDim*i) = Ba(1:nStress,1:nDim)
          Gmat(1:nDim**2,nDim*(i-1)+1:nDim*i) = Ga(1:nDim**2,1:nDim)
        end do                             ! end of nodal point loop


        ! transpose the vector field matrix operators
        NmatT       = transpose(Nmat)
        BmatT       = transpose(Bmat)
        GmatT       = transpose(Gmat)

        !!!!!!!!!!!!!! COMPLETE ELEMENT RELATED OPERATIONS !!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!! CONSTITUTIVE MODEL !!!!!!!!!!!!!!!!!!!!!!

        ! calculate deformation gradient and deformation tensors
        F(1:nDim,1:nDim) = ID + matmul(uNode,dNdX)

        if (analysis .eq. 'PE')  F(3,3) = one

        ! calculate material point jacobian (volume change)
        detF    = det(F)
        FInv    = inv(F)
        FInvT   = transpose(FInv)

        !! definition of modified deformation gradient, F-bar
        if (fbarFlag .eq. .true.) then
          if ( (jtype .eq. 3) .and. (nInt .eq. 8) ) then
            ! plane strain linear quad element
            Fbar    = (detF0/detF)**(third) * F
            resFac  = (detF0/detF)**(-two/three)
            tanFac1 = (detF0/detF)**(-one/three)
            tanFac2 = (detF0/detF)**(-two/three)

          else if ( (jtype .eq. 7) .and. (nInt .eq. 4) ) then
            ! three-dimensional linear hex element
            Fbar(1:nDim,1:nDim) = (detF0/detF)**(half)*F(1:nDim,1:nDim)
            Fbar(3,3)           = one
            resFac              = (detF0/detF)**(-half)
            tanFac1             = one
            tanFac2             = (detF0/detF)**(-half)

          else
            ! standard F for all other available elements
            Fbar    = F
            resFac  = one
            tanFac1 = one
            tanFac2 = one

            call msg%ferror( flag=warn, src='elem_nlmech',
     &          msg='F-bar is not available: ', ivec=[jtype, nInt])
            call xit
          end if
        else
          ! set F-bar = F if fbarFlag is .false. for all element
          Fbar    = F
          resFac  = one
          tanFac1 = one
          tanFac2 = one
        end if


    !     ! interpolate the field variables at the integration point
    !     ! (CAUTION: this not yet tested or used)
    !     do k = 1, npredf
    !       fieldVar(k)   = dot_product( Nxi,
    !  &                    reshape( fieldNode(k,1:nNode), [nNode] ) )
    !       dfieldVar(k)  = dot_product( Nxi,
    !  &                    reshape( dfieldNode(k,1:nNode), [nNode] ) )
    !     end do

        ! call material point subroutine (UMAT) for specific material
        if (matID .eq. 1) then
          call mat_NeoHookean(kstep,kinc,time,dtime,nDim,analysis,
     &            nstress,nNode,jelem,coords,intpt,props,nprops,
     &            jprops,njprops,Fbar,svars,nsvars,fieldVar,dfieldVar,
     &            npredf,stressPK2,Dmat,Cmat)

        else if (matID .eq. 2) then
          call mat_ArrudaBoyce(kstep,kinc,time,dtime,nDim,analysis,
     &            nstress,nNode,jelem,coords,intpt,props,nprops,
     &            jprops,njprops,Fbar,svars,nsvars,fieldVar,dfieldVar,
     &            npredf,stressPK2,Dmat,Cmat)
        else
          call msg%ferror( flag=error, src='elem_nlmech',
     &                    msg='Wrong material ID.', ia=matID )
          call xit
        end if


        !!!!!!!!!!!!!!!!!!!! END CONSTITUTIVE MODEL !!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!! FORM ADDITIONAL ELEMENT OPERATORS !!!!!!!!!!!!!!!!

        ! form the [SIGMA_S] matrix for geometric stiffness
        call voigtVectorScatter(stressPK2, stressTensorPK2)

        do i = 1, nDim
          do j = 1, nDim
              SIGMA_S(nDim*(i-1)+1:nDim*i,nDim*(j-1)+1:nDim*j)
     &                    = stressTensorPK2(i,j)*ID
          end do
        end do

        ! form [SIGMA_F] matrix for material stiffness
        do i = 1, nNode
          do j = 1, nNode
              if (i .eq. j) then                    ! banded diagonal
                  SIGMA_F(nDim*(i-1)+1:nDim*i,nDim*(j-1)+1:nDim*j)
     &                              = Fbar(1:nDim,1:nDim)
              end if
          end do
        end do

        BNLmat  = matmul(Bmat,transpose(SIGMA_F))
        BNLmatT = transpose(BNLmat)

        !!!!!!!!!! END FORMING ADDITIONAL ELEMENT OPERATORS !!!!!!!!!!!!


        ! form the stiffness matrix and residual vector
        ! tanFac1 and resFac1 will perform modification on the tangent
        ! and residual depending on the type of element being used


        !!!!!!!!!!!!!! TANGENT MATRIX AND RESIDUAL VECTOR !!!!!!!!!!!!!!
        Ru  = Ru - w(intpt) * detJ * resFac *
     &      matmul( transpose(matmul(Bmat,SIGMA_F)), stressPK2 )


        Kuu = Kuu + w(intpt) * detJ * tanFac1 *
     &      (
     &        matmul( matmul( GmatT, SIGMA_S ), Gmat )
     &        + matmul( BNLmatT, matmul(Dmat,BNLmat) )
     &      )


        ! tangent modification for F-bar
        if (fbarFlag .eq. .true.) then

          ! form fourth-order QR0 and QR tensor
          QR0Tensor = zero
          QRTensor  = zero

        ! perform F-bar modification on linear hex element
          if ((jtype .eq. 3) .and. (nInt .eq. 8)) then
            do i = 1,nDim
              do j = 1,nDim
                do k = 1,nDim
                  do l = 1,nDim
                    do m = 1,nDim
                      do n = 1,nDim
                        do p = 1,nDim
                          do q = 1,nDim

                            QR0Tensor(i,j,k,l) = QR0Tensor(i,j,k,l)
     &                        + third * F0InvT(k,l) *
     &                          (
     &                            Fbar(i,p) * Cmat(p,j,m,n)
     &                            * Fbar(q,m) * Fbar(q,n)
     &                            - Fbar(i,q) * stressTensorPK2(q,j)
     &                          )

                            QRTensor(i,j,k,l) = QRTensor(i,j,k,l)
     &                        + third * FInvT(k,l) *
     &                          (
     &                            Fbar(i,p) * Cmat(p,j,m,n)
     &                            * Fbar(q,m) * Fbar(q,n)
     &                            - Fbar(i,q) * stressTensorPK2(q,j)
     &                          )
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do

          ! do F-bar modification on linear quad element
          else if ( (jtype .eq. 7) .and. (nInt .eq. 4) ) then
            do i = 1,nDim
              do j = 1,nDim
                do k = 1,nDim
                  do l = 1,nDim
                    do m = 1,nDim
                      do n = 1,nDim
                        do p = 1,nDim
                          do q = 1,nDim
                            QR0Tensor(i,j,k,l) = QR0Tensor(i,j,k,l)
     &                        + half * Fbar(i,p) * Cmat(p,j,q,n)
     &                        * Fbar(m,n) * Fbar(m,q) * F0InvT(k,l)

                            QRTensor(i,j,k,l) = QRTensor(i,j,k,l)
     &                        + half * Fbar(i,p) * Cmat(p,j,q,n)
     &                        * Fbar(m,n) * Fbar(m,q) * FInvT(k,l)
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do

          end if

          ! reshape QR and QR0 tensor into matrix form
          call unsymmMatrix(QR0Tensor,QR0mat)
          call unsymmMatrix(QRTensor,QRmat)

          ! modify the element tangent matrix
          Kuu = Kuu + w(intPt) * detJ * tanFac2  *
     &              (
     &              matmul(transpose(Gmat), matmul(QR0mat,Gmat0))
     &              - matmul(transpose(Gmat), matmul(QRmat,Gmat))
     &              )
        end if

      !!!!!!!!!!!!! END TANGENT MATRIX AND RESIDUAL VECTOR !!!!!!!!!!!!!

      end do                           ! end of integration point loop

      ! body force and surface load can be added using dummy elements

      ! assign the element stiffness matrix to abaqus-defined variable
      amatrx(1:NDOFEL,1:NDOFEL) = Kuu(1:nDOFEL,1:nDOFEL)
      rhs(1:NDOFEL,1)           = Ru(1:nDOFEL,1)

      end subroutine elem_nlmech

      end module finite_strain_element

! **********************************************************************
! ****************** ABAQUS USER ELEMENT SUBROUTINE ********************
! **********************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)

! **********************************************************************
!     This subroutine is called by Abaqus with above arguments
!     for each user elements defined in an Abaqus model. Users are
!     responsible for programming the element tangent/ stiffness
!     matrix and residual vectors which will be then assembled and
!     solved by Abaqus after applying the boundary conditions.
! **********************************************************************

      use global_parameters
      use error_logging
      use post_processing
      use finite_strain_element

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
      integer             :: nDim, nStress
      character(len=2)    :: analysis
      character(len=8)    :: abqProcedure
      logical             :: nlgeom


      logical, parameter  :: dbgMode = .false.
      integer             :: lenJobName,lenOutDir
      character(len=256)  :: outDir
      character(len=256)  :: jobName
      character(len=512)  :: errFile, dbgFile
      type(logger)        :: msg


      ! initialize the output variables to be defined for Abaqus subroutine
      energy        = zero
      amatrx        = zero
      rhs(:,nrhs)   = zero


      ! open logs files for the current job from Abaqus run
      if (dbgMode .eq. .false.) then
        call getJobName(jobName,lenJobName)
        call getOutDir(outDir,lenOutDir)
        errFile = trim(outDir)//'\aaERR_'//trim(jobName)//'.dat'
        dbgFile = trim(outDir)//'\aaDBG_'//trim(jobName)//'.dat'
        call msg%fopen( errfile=errFile, dbgfile=dbgFile )
      end if


      ! change the LFLAGS criteria as needed (check abaqus UEL manual)
      if ((lflags(1) .eq. 1) .or. (lflags(1) .eq. 2)) then
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

      else if ((jtype .ge. 5).and.(jtype .le. 8)) then
        nDim      = 2
        analysis  = 'PE'            ! plane strain analysis
        nStress   = 3
      else
        call msg%ferror( flag=error, src='UEL',
     &            msg='Element is unavailable.', ia=jtype )
        call xit
      end if


      nInt      = jprops(1)
      matID     = jprops(3)
      nPostVars = jprops(4)


      ! array containing variables for post-processing
      if (.not. allocated(globalPostVars)) then
        allocate(globalPostVars(numElem,nInt,nPostVars))

        ! print job-related information the first time
        call msg%finfo('---------------------------------------')
        call msg%finfo('------- Abaqus FINITE STRAIN UEL ------')
        call msg%finfo('-------- PK-II TOTAL LAGRANGIAN -------')
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
       call elem_nlmech(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & NDIM,ANALYSIS,NSTRESS,NINT)

      END SUBROUTINE UEL

! **********************************************************************
! ************** ABAQUS USER OUTPUT VARIABLES SUBROUTINE ***************
! **********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     & NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     & JMAC,JMATYP,MATLAYO,LACCFLA)

! **********************************************************************
!     This subroutine is called by Abaqus at each material point (int pt)
!     to obtain the user defined output variables for standard Abaqus
!     elements. We used an additional layer of standard Abaqus elements
!     with same topology (same number of nodes and int pts) on top of
!     the user elements with an offset in the numbering between the user
!     elements and standard elements. This number is defined in the
!     post_processing module and should match with Abaqus input file.
! **********************************************************************

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