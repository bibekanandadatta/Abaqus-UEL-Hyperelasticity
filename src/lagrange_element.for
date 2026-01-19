! **********************************************************************
! ********************** LAGRANGE ELEMENT MODULE ***********************
! **********************************************************************
!   collection of subroutines to perform element realted calculations
!            subroutines are arranged in alphabetical order.
! **********************************************************************
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2024. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************

      module lagrange_element

      use global_parameters

      type, public  :: element
        integer                     :: nDim
        character(len=8)            :: analysis
        integer                     :: nNode
        integer                     :: nInt
      end type element

! **********************************************************************
!     interface to evaluate shape functions at an interior coordinate
! **********************************************************************

      interface

        module subroutine calcInterpFunc(elem,xiCoord,Nxi,dNdxi)
          use global_parameters, only: wp
          implicit none
          type(element), intent(in)   :: elem
          real(wp), intent(in)        :: xiCoord(:)
          real(wp), intent(out)       :: Nxi(:), dNdxi(:,:)
        end subroutine calcInterpFunc

        module subroutine faceNodes(elem,face,nFaceNodes,list)
          implicit none
          type(element), intent(in)   :: elem
          integer, intent(in)         :: face
          integer, intent(out)        :: nFaceNodes
          integer, intent(out)        :: list(*)
        end subroutine faceNodes

      end interface

      end module lagrange_element

! **********************************************************************
! ****************** INTERPOLATION FUNCTION SUBMODULE ******************
! **********************************************************************
!  available elements:    (a) 1D bar/truss element (2 and 3 nodes)
!                         (b) 2D tri elements (3 and 6 nodes)
!                         (c) 2D quad elements (4 and 8 nodes)
!                         (d) 3D tet elements (4 and 10 nodes)
!                         (e) 3D hex elements (8 and 20 nodes)
! **********************************************************************

      submodule (lagrange_element) interpolation

      contains

      module subroutine calcInterpFunc(elem,xiCoord,Nxi,dNdxi)
      ! driver subroutine to calculate shape functions and derivatives

        use global_parameters, only: wp

        implicit none

        type(element), intent(in)   :: elem
        real(wp), intent(in)        :: xiCoord(:)
        real(wp), intent(out)       :: Nxi(:), dNdxi(:,:)

        if (elem%nDim  .eq.  1) then
          call interpFunc1(elem%nNode, xiCoord, Nxi, dNdxi)
        else if (elem%nDim  .eq.  2) then
          call interpFunc2(elem%nNode, xiCoord, Nxi, dNdxi)
        else if (elem%nDim  .eq.  3) then
          call interpFunc3(elem%nNode, xiCoord, Nxi, dNdxi)
        end if

      end subroutine calcInterpFunc

! **********************************************************************

      subroutine interpFunc1(nNode,xiCoord,Nxi,dNdxi)
      ! this subroutine calculates shape function of 1D elements
      ! available 1D elements are: 2 and 3 node bar/truss

      ! Nxi(i)          = shape function of node i at the intpt.
      ! dNdxi(i,j)      = derivative wrt j direction of shape fn of node i

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)         :: nNode
      real(wp), intent(in)        :: xiCoord(:)
      real(wp), intent(out)       :: Nxi(:), dNdxi(:,:)
      real(wp)                    :: xi

      ! location in the master element
      xi    = xiCoord(1)

      if (nNode  .eq.  2) then      ! 2 node linear bar element
        ! shape functions
        Nxi(1) = half*(one - xi)
        Nxi(2) = half*(one + xi)

        ! the first derivatives of the shape functions dN/dxi (2x1)
        dNdxi(1,1)  = -half
        dNdxi(2,1)  = half

      else if (nNode  .eq.  3) then  ! 3 node quadratic bar element
        ! shape functions
        Nxi(1)  = -half*xi*(one - xi)
        Nxi(2)  = one-xi**two
        Nxi(3)  = half*xi*(one + xi)

        ! the first derivatives of the shape functions dN/dxi (3x1)
        dNdxi(1,1)  = -half+xi
        dNdxi(2,1)  = -two*xi
        dNdxi(3,1)  = half+xi

      else
        call msg%ferror(flag=error, src='interpFunc1',
     &                  msg='Element is unavailable.', ia=nNode)
        return
      end if

      end subroutine interpFunc1

! **********************************************************************

      subroutine interpFunc2(nNode,xiCoord,Nxi,dNdxi)
      ! this subroutine calculates shape function of 2D elements
      ! available 2D elements are: 3 node tri, 6 node tri, 4 node quad, 8 node quad

      ! Nxi(i)          = shape function of node i at the intpt.
      ! dNdxi(i,j)      = derivative wrt j direction of shape fn of node i

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)         :: nNode
      real(wp), intent(in)        :: xiCoord(:)
      real(wp), intent(out)       :: Nxi(:), dNdxi(:,:)
      real(wp)                    :: xi, eta, lam

      ! location in the master element
      xi    = xiCoord(1)
      eta   = xiCoord(2)

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

      if (nNode .eq. 3) then        ! 3-noded tri3 linear element
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

      else if (nNode .eq. 6) then    ! 6-noded quadratic tri6 element
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

      else if (nNode .eq. 4) then    ! 4-noded bilinear quad4 element
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

      else if (nNode .eq. 8) then    ! 8-noded serendipity quad8 element
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
        call msg%ferror(flag=error, src='interpFunc2',
     &                  msg='Element is unavailable.', ia=nNode)
        return
      end if

      end subroutine interpFunc2

! **********************************************************************

      subroutine interpFunc3(nNode,xiCoord,Nxi,dNdxi)
      ! this subroutine calculates shape function of 3D elements
      ! available 3D elements are: 4 node tet, 10 node ter, 8 node hex.

      ! Nxi(i)          = shape function of node i at the intpt.
      ! dNdxi(i,j)      = derivative wrt j direction of shape fn of node i

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)         :: nNode
      real(wp), intent(in)        :: xiCoord(:)
      real(wp), intent(out)       :: Nxi(:), dNdxi(:,:)
      real(wp)                    :: xi, eta, zeta, lam

      ! Nxi(i)          = shape function of node i at the intpt.
      ! dNdxi(i,j)      = derivative wrt j direction of shape fn of node i

      ! location in the master element
      xi    = xiCoord(1)
      eta   = xiCoord(2)
      zeta  = xiCoord(3)

      Nxi   = zero
      dNdxi = zero

      if (nNode .eq. 4) then      ! 4-noded linear tet4 element
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

      else if (nNode .eq. 10) then  ! 10-noded quadratic tet10 element

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

      else if(nNode .eq. 8) then   ! 8-noded trilinear hex8 element

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

      else if (nNode .eq. 20) then   ! 20-noded serendipity hex20 element
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
        call msg%ferror(flag=error, src='interpFunc3',
     &                  msg='Element is unavailable.', ia=nNode)
        return
      end if

      end subroutine interpFunc3

      end submodule interpolation

! **********************************************************************
! **********************************************************************