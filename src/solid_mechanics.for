! **********************************************************************
! solid mechanics module with subroutines to perform tensor manipulation
! **********************************************************************
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2024. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************

      module solid_mechanics

      !! no of symmetric and unSymmmetric stress tensor components (3D)
      integer, parameter    :: nSymm = 6, nUnsymmm = 9

      contains

! **********************************************************************

      subroutine symtangent2matrix(C,Dmat)
      ! this subroutine maps the fourth order material/spatial tangent
      ! tensor (3x3x3x3) to a 2nd order stiffness tensor (6x6) using
      ! voigt notation: 11> 1, 22> 2, 33> 3, 23/32> 4, 13/31> 5, 12/21> 6

        use global_parameters, only: wp, zero

        implicit none

        real(wp), intent(in)    :: C(3,3,3,3)
        real(wp), intent(out)   :: Dmat(nSymm,nSymm)
        integer                 :: Voigt(nSymm,2)
        integer                 :: i, j, k, l, rw, cl

        ! Voigt convetion: (1,1) (2,2) (3,3) (2,3) (1,3) (1,2)
        Voigt = reshape(  [ 1, 2, 3, 2, 1, 1,
     &                      1, 2, 3, 3, 3, 2 ], shape(Voigt) )

        Dmat = zero

        do rw = 1, nSymm
          do cl = 1, nSymm
            i = Voigt(rw,1)
            j = Voigt(rw,2)
            k = Voigt(cl,1)
            l = Voigt(cl,2)

            Dmat(rw,cl) = C(i,j,k,l)
          end do
        end do

      end subroutine symtangent2matrix

! **********************************************************************

      subroutine symtensor2vector(ATens, AVect)
      ! this subroutine maps a symmetric tensor to a vector
      ! for unSymmmetric tensor you can use "reshape" function
      ! for 2D (plane)  : ATens(2,2) >> AVect(3,1)
      ! for 3D          : ATens(3,3) >> AVect(6,1)

        use global_parameters, only: wp

        implicit none

        real(wp), intent(in)    :: ATens(:,:)
        real(wp), intent(out)   :: AVect(:,:)
        integer                 :: nDim, nStress, i

        nDim    = size(ATens,1)
        nStress = size(AVect,1)

        do i = 1, nDim
          AVect(i,1) = ATens(i,i)
        end do

        if (nDim .eq. 2) then
          AVect(3,1) = ATens(1,2)

        else if (nDim .eq. 3) then
          AVect(4,1) = ATens(2,3)
          AVect(5,1) = ATens(1,3)
          AVect(6,1) = ATens(1,2)
        end if

      end subroutine symtensor2vector

! **********************************************************************
      
      subroutine unsymtensor2matrix(tensor,matrix)
      ! this subroutine maps the fourth order tensor to a rank-2 matrix

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

      end subroutine unsymtensor2matrix

! **********************************************************************

      subroutine vector2symtensor(AVect, ATens)
      ! this subroutine transforms a Voigt vector to symmetric tensor
      ! for 2D (plane)  : AVect(3,1) >> ATens(2,2)
      ! for 3D          : AVect(6,1) >> ATens(3,3)

        use global_parameters, only: wp

        implicit none

        real(wp), intent(in)    :: AVect(:,:)
        real(wp), intent(out)   :: ATens(:,:)
        integer                 :: nDim, nStress, i

        nDim    = size(ATens,1)
        nStress = size(AVect,1)

        do i = 1, nDim              ! direct / axial stresses
          ATens(i,i) = AVect(i,1)
        end do

        ! shear stress mapping
        if (nDim .eq. 2) then
          ATens(1,2) = AVect(3,1)
          ATens(2,1) = ATens(1,2)

        else if (nDim .eq. 3) then
          ATens(2,3) = AVect(4,1)
          ATens(1,3) = AVect(5,1)
          ATens(1,2) = AVect(6,1)
          ATens(2,1) = ATens(1,2)
          ATens(3,1) = ATens(1,3)
          ATens(3,2) = ATens(2,3)
        end if

      end subroutine vector2symtensor

! **********************************************************************

      subroutine voigtAugment(vect2D, vect3D)
      ! this subroutine augemnts a 3x1 (plane) or 4x1 (axisymmetric)
      ! array to a 6x1 Voigt array of 3D dimensional case

        use global_parameters, only: wp, zero

        implicit none

        real(wp), intent(in)    :: vect2D(:,:)
        real(wp), intent(out)   :: vect3D(nSymm,1)
        integer                 :: nStress

        nStress = size(vect2D,1)
        vect3D  = zero

        vect3D(1,1) = vect2D(1,1)
        vect3D(2,1) = vect2D(2,1)

        if (nStress .eq. 3) then      ! plane strain/stress
          vect3D(6,1) = vect2D(3,1)
        else if (nStress .eq. 4) then  ! axisymmetry
          vect3D(3,1) = vect2D(3,1)
          vect3D(6,1) = vect2D(4,1)
        else if (nStress .eq. nSymm) then
          vect3D = vect2D
        end if

      end subroutine voigtAugment

! **********************************************************************

      subroutine voigtTruncate(vect3D, vect2D)
      ! this subroutine truncates a 6x1 Voigt array
      ! to a 3x1 (plane) or 4x1 (axisymmetry) Voigt array

        use global_parameters, only: wp, zero

        implicit none

        real(wp), intent(in)    :: vect3D(nSymm,1)
        real(wp), intent(out)   :: vect2D(:,:)
        integer                 :: nStress

        nStress = size(vect2D,1)
        vect2D  = zero

        vect2D(1,1) = vect3D(1,1)
        vect2D(2,1) = vect3D(2,1)

        if (nStress .eq. 3) then
          vect2D(3,1) = vect3D(6,1)
        else if (nStress .eq. 4) then
          vect2D(3,1) = vect3D(3,1)
          vect2D(4,1) = vect3D(6,1)
        end if

      end subroutine voigtTruncate

      end module solid_mechanics

! **********************************************************************
! **********************************************************************