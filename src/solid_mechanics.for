! **********************************************************************
! solid mechanics module with subroutines to perform tensor manipulation
! **********************************************************************
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2024. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************

      module solid_mechanics

      ! no of symmetric and unSymmetric stress tensor components (3D)
      integer, parameter    :: nSymm = 6, nUnsymm = 9

      contains

! **********************************************************************

      subroutine unsymmMatrix(tensor,matrix)
      ! This subroutine reshapes a general the fourth order tensor to a
      ! rank-2 matrix form. For example, this operation is necessary to 
      ! reshape 4th-order first elasticity tensor dP/dF (derivative of 
      ! PK-I stress wrt deformation gradient) to a matrix form.
      ! use this for 3D and 2D plane strain/stress cases (not axisymmetry)

        use global_parameters, only: wp, zero

        implicit none

        real(wp), intent(in)  :: tensor(:,:,:,:)
        real(wp), intent(out) :: matrix(size(tensor,1)*size(tensor,2),
     &                                  size(tensor,3)*size(tensor,4))
        integer               :: i, j, k, l, m, n, nDim

        nDim    = size(tensor,1)
        matrix  = zero

        do i = 1,nDim
          do j = 1,nDim
            do k = 1,nDim
              do l = 1,nDim
                m           = (i-1)*nDim + j
                n           = (k-1)*nDim + l
                matrix(m,n) = tensor(j,i,l,k)
              end do
            end do
          end do
        end do

      end subroutine unsymmMatrix

! **********************************************************************

      subroutine unsymmVector(tensor,vector)
      ! This subroutine reshapes a general 2nd order tensor (unsymmetric)
      ! to a column-ordered vector form. For example, this can be used
      ! to reshape PK-I stress tensor to a vector form.
      ! use this for 3D and 2D plane strain/stress cases (not axisymmetry)

        use global_parameters, only: wp, zero

        implicit none

        real(wp), intent(in)  :: tensor(:,:)      ! dimension: nDim x nDim
        real(wp), intent(out) :: vector(:,:)      ! dimension: nDim*nDim x 1
        integer               :: nDim

        nDim      = size(tensor,1)

        vector    = reshape(tensor(1:nDim,1:nDim), shape=[nDim*nDim, 1])

      end subroutine unsymmVector

! **********************************************************************

      subroutine unsymmVectorTruncate(vect3D,vect2D)
      ! This subroutine reshapes a 9x1 unsymmetric vector (deformation gradient 
      ! or PK-I stress) to a 4x1 or 5x1 vector based on the type of analysis
      ! output order: P11, P21, P12, P22, P33

        use global_parameters, only: wp, zero

        implicit none

        real(wp), intent(in)    :: vect3D(nUnsymm,1)
        real(wp), intent(out)   :: vect2D(:,:)
        integer                 :: nStress

        nStress     = size(vect2D,1)
        vect2D      = zero

        if (nStress .eq. 9) then          ! 3D analysis
          vect2D = vect3D

        else if (nStress .eq. 5) then     ! axisymmetric analysis
          vect2D(1:2,1) = vect3D(1:2,1)
          vect2D(3:4,1) = vect3D(4:5,1)
          vect2D(5,1)   = vect3D(9,1)

        else if (nStress .eq. 4) then     ! plane strain/ plane stress analysis
          vect2D(1:2,1) = vect3D(1:2,1)
          vect2D(3:4,1) = vect3D(4:5,1)
        end if

      end subroutine unsymmVectorTruncate

! **********************************************************************

      subroutine vector2tensor(vector,tensor)
      ! This subroutine reshapes a general a vector to a tensor using
      ! column-ordering convention. For example, this can be used to
      ! reshape a PK-I stress vector to a tensor form.

        use global_parameters, only: wp, zero

        implicit none

        real(wp), intent(in)  :: vector(:,:)    ! dimension: nDim*nDim x 1
        real(wp), intent(out) :: tensor(:,:)    ! dimension: nDim x nDim
        integer               :: nDim

        nDim    = size(tensor,1)

        tensor  = reshape(vector(1:nDim*nDim,1), shape=[nDim, nDim])

      end subroutine vector2tensor

! **********************************************************************
      
      subroutine voigtMatrix(tensor,voigtMat)
      ! this subroutine maps a symmetric fourth order material/spatial 
      ! tangent tensor (3x3x3x3) to a 2nd order stiffness matrix (6x6) using
      ! voigt notation: 11> 1, 22> 2, 33> 3, 23/32> 4, 13/31> 5, 12/21> 6

        use global_parameters, only: wp, zero

        implicit none

        real(wp), intent(in)    :: tensor(3,3,3,3)
        real(wp), intent(out)   :: voigtMat(nSymm,nSymm)
        integer                 :: Voigt(nSymm,2)
        integer                 :: i, j, k, l, rw, cl

        ! Voigt convention notation: (1,1) (2,2) (3,3) (2,3) (1,3) (1,2)
        Voigt = reshape(  [ 1, 2, 3, 2, 1, 1,
     &                      1, 2, 3, 3, 3, 2 ], shape(Voigt) )

        voigtMat = zero

        do rw = 1, nSymm
          do cl = 1, nSymm
            i = Voigt(rw,1)
            j = Voigt(rw,2)
            k = Voigt(cl,1)
            l = Voigt(cl,2)

            voigtMat(rw,cl) = tensor(i,j,k,l)
          end do
        end do

      end subroutine voigtMatrix

! **********************************************************************

      subroutine voigtMatrixTruncate(voigtMat,Dmat)
      ! This subroutine truncates a 6 x 6 Voigt matrix to 3 x 3 (plane 
      ! strain and plane stress) or 4 x 4 (axisymmetry) D-matrix for 
      ! finite element operation based on the type of analysis.

        use global_parameters, only: wp, zero

        real(wp), intent(in)    :: voigtMat(:,:)    ! dimension: 6x6
        real(wp), intent(out)   :: Dmat(:,:)        ! dimension: 3x3 or 4x4 or 6x6

        integer                 :: nstress
        integer                 :: nsize


        nsize       = size(voigtMat,1)
        nStress     = size(Dmat,1)
        

        ! truncating a full voigt matrix (6x6) for 2D analysis
        if (nsize .eq. 6) then 
          if (nStress .eq. 6) then              ! 3D analysis
            Dmat            = voigtMat

          else if (nStress .eq. 4) then         ! axisymmetry
            Dmat(1:3,1:3)   = voigtMat(1:3,1:3)
            Dmat(1:3,4)     = voigtMat(1:3,6)
            Dmat(4,1:3)     = voigtMat(6,1:3)
            Dmat(4,4)       = voigtMat(6,6)

          else if (nStress .eq. 3) then         ! plane strain
            Dmat(1:2,1:2)   = voigtMat(1:2,1:2)
            Dmat(1:2,3)     = voigtMat(1:2,6)
            Dmat(3,1:2)     = voigtMat(6,1:2)
            Dmat(3,3)       = voigtMat(6,6)
          end if

        ! truncating a general 2D voigt matrix (4x4) to a 3x3 one
        ! this operation is performed for plane stress case
        else if (nsize .eq. 4) then
          if (nStress .eq. 3) then
            Dmat(1:2,1:2)   = voigtMat(1:2,1:2)
            Dmat(1:2,3)     = voigtMat(1:2,4)
            Dmat(3,1:2)     = voigtMat(4,1:2)
            Dmat(3,3)       = voigtMat(4,4)
          end if
        end if
      end subroutine voigtMatrixTruncate

! **********************************************************************

      subroutine voigtVector(tensor,vector)
      ! this subroutine reshapes a symmetric 2nd order tensor to a vector.
      ! For example, it can reshape, Lagrange strain, Euler strain, PK-II
      ! stress, Cauchy stress tensors to vector form using the Voigt notation.
      ! voigt notation: 11> 1, 22> 2, 33> 3, 23/32> 4, 13/31> 5, 12/21> 6
      
      ! for 2D (plane)  : tensor(2,2) >> vector(3,1)
      ! for 3D          : tensor(3,3) >> vector(6,1)

        use global_parameters, only: wp

        implicit none

        real(wp), intent(in)    :: tensor(:,:)      ! dimension: nDim x nDim
        real(wp), intent(out)   :: vector(:,:)      ! dimension: nstress x 1
        integer                 :: nDim, nStress, i

        nDim    = size(tensor,1)
        nStress = size(vector,1)

        do i = 1, nDim                      ! axial/ direct components
          vector(i,1) = tensor(i,i)
        end do

        if (nDim .eq. 2) then               ! shear strain components
          vector(3,1) = tensor(1,2)

        else if (nDim .eq. 3) then
          vector(4,1) = tensor(2,3)
          vector(5,1) = tensor(1,3)
          vector(6,1) = tensor(1,2)
        end if

      end subroutine voigtVector

! **********************************************************************

      subroutine voigtVectorAugment(vect2D,vect3D)
      ! this subroutine augemnts a 3x1 (plane) or 4x1 (axisymmetric)
      ! vector to a 6x1 vector using Voigt notation convention.
      ! voigt notation: 11> 1, 22> 2, 33> 3, 23/32> 4, 13/31> 5, 12/21> 6.

        use global_parameters, only: wp, zero
        use error_logging

        implicit none

        real(wp), intent(in)    :: vect2D(:,:)
        real(wp), intent(out)   :: vect3D(:,:)
        integer                 :: nStress
        type(logger)            :: msg

        if (size(vect3D(:,1)) .ne. nSymm) then
          call msg%ferror(flag=error, src='voigtVectorTruncate',
     &         msg='Size of the output argument should be 6x1 vector.')
          return 
        end if

        nStress = size(vect2D,1)
        vect3D  = zero

        vect3D(1,1) = vect2D(1,1)
        vect3D(2,1) = vect2D(2,1)

        if (nStress .eq. 3) then            ! plane strain/stress
          vect3D(6,1) = vect2D(3,1)

        else if (nStress .eq. 4) then       ! axisymmetry
          vect3D(3,1) = vect2D(3,1)
          vect3D(6,1) = vect2D(4,1)

        else if (nStress .eq. nSymm) then   ! 3D case
          vect3D = vect2D
        end if

      end subroutine voigtVectorAugment

! **********************************************************************

      subroutine voigtVectorScatter(vector,tensor)
      ! this subroutine reshapes a Voigt vector to symmetric tensor.
      ! This operation is useful to reshape symmetric strain or stress
      ! vector to a matrix form for any consequent operation.
      ! voigt notation: 11> 1, 22> 2, 33> 3, 23/32> 4, 13/31> 5, 12/21> 6

      ! for 2D (plane)  : vector(3x1) >> tensor(2x2)
      ! for 3D          : vector(6x1) >> tensor(3x3)

        use global_parameters, only: wp

        implicit none

        real(wp), intent(in)    :: vector(:,:)
        real(wp), intent(out)   :: tensor(:,:)
        integer                 :: nDim, nStress, i

        nDim    = size(tensor,1)
        nStress = size(vector,1)

        ! direct / axial stresses
        do i = 1, nDim              
          tensor(i,i) = vector(i,1)
        end do

        ! shear stress mapping
        if (nDim .eq. 2) then
          tensor(1,2) = vector(3,1)
          ! fill out the symmetric components
          tensor(2,1) = tensor(1,2)

        else if (nDim .eq. 3) then
          tensor(2,3) = vector(4,1)
          tensor(1,3) = vector(5,1)
          tensor(1,2) = vector(6,1)
          ! fill out the symmetric components
          tensor(3,2) = tensor(2,3)
          tensor(3,1) = tensor(1,3)
          tensor(2,1) = tensor(1,2)
        end if

      end subroutine voigtVectorScatter

! **********************************************************************

      subroutine voigtVectorTruncate(vect3D, vect2D)
      ! this subroutine truncates a 6x1 Voigt array to a 3x1 (plane)
      ! or 4x1 (axisymmetry) Voigt array for finite element operations
      ! based on the type of analysis being performed (such as,
      ! plane strain/ stress or axisymmetric). If the analysis is 3D
      ! it just returns the same vector.
      ! voigt notation: 11> 1, 22> 2, 33> 3, 23/32> 4, 13/31> 5, 12/21> 6

        use global_parameters, only: wp, zero
        use error_logging

        implicit none

        real(wp), intent(in)    :: vect3D(:,:)
        real(wp), intent(out)   :: vect2D(:,:)
        integer                 :: nStress
        type(logger)            :: msg

        if (size(vect3D(:,1)) .ne. nSymm) then
          call msg%ferror(flag=error, src='voigtVectorTruncate',
     &         msg='Size of the input argument should be 6x1 vector.')
          return 
        end if

        nStress   = size(vect2D,1)
        vect2D    = zero

        if (nStress .eq. 6) then          ! 3D analysis
          vect2D  = vect3D

        else if (nStress .eq. 3) then     ! plane stress/ plane strain
          vect2D(1,1) = vect3D(1,1)
          vect2D(2,1) = vect3D(2,1)
          vect2D(3,1) = vect3D(6,1)

        else if (nStress .eq. 4) then     ! axisymmetric analysis
          vect2D(1,1) = vect3D(1,1)
          vect2D(2,1) = vect3D(2,1)
          vect2D(3,1) = vect3D(3,1)
          vect2D(4,1) = vect3D(6,1)
        end if

      end subroutine voigtVectorTruncate

      end module solid_mechanics

! **********************************************************************
! **********************************************************************