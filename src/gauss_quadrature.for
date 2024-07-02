! **********************************************************************
! ******************** GAUSSIAN QUADRATURE MODULE **********************
! **********************************************************************
! integration schemes:  (a) reduced and full integration: bar elements
!                       (b) full integration: tri and tet elements
!                       (c) reduced and full integration: quad and hex elements
! **********************************************************************
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2024. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************

      module gauss_quadrature

      private :: gaussQuadrtr1, gaussQuadrtr2, gaussQuadrtr3
      public  :: getGaussQuadrtr

      contains

      subroutine getGaussQuadrtr(elem,w,xi)
      ! driver subroutine for gaussian quadrature tabluation

        use global_parameters, only: wp
        use lagrange_element, only: element

        implicit none

        type(element)           :: elem
        real(wp), intent(out)   :: w(:), xi(:,:)

        if (elem%nDim  .eq.  1) then
          call gaussQuadrtr1(elem%nNode,elem%nInt,w,xi)
        else if (elem%nDim  .eq.  2) then
          call gaussQuadrtr2(elem%nNode,elem%nInt,w,xi)
        else if (elem%nDim  .eq.  3) then
          call gaussQuadrtr3(elem%nNode,elem%nInt,w,xi)
        end if

      end subroutine getGaussQuadrtr

! **********************************************************************

      subroutine gaussQuadrtr1(nNode,nInt,w,xi)

      ! this subroutines returns the weights and
      ! gauss point coordinate of 2D Lagrangian elements
      ! currently supports: nInt = 1, 2, 3  (for bar2 and bar3 elements)

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)     :: nNode, nInt
      real(wp), intent(out)   :: w(:), xi(:,:)
      type(logger)            :: msg

      w = zero
      xi = zero

      ! 1D bar2 and bar3 elements
      if ((nNode .eq. 2) .or. (nNode .eq. 3)) then
        if (nInt .eq. 1) then
          w = zero
          xi = two

        else if (nInt  .eq.  2) then
          w(1:2)  = one
          xi(1,1) = -sqrt(third)
          xi(2,1) = sqrt(third)

        else if (nInt .eq. 3) then
          w(1) = five/nine
          w(2) = eight/nine
          w(3) = five/nine
          xi(1,1) = -sqrt(three/five)
          xi(2,1) = zero
          xi(3,1) = sqrt(three/five)
        else
          call msg%ferror(flag=error,src='gaussQuadrtr1',
     &         msg='Invalid Gauss points for BAR element.', ia=nInt)
          return
        end if

      else
        call msg%ferror(flag=error,src='gaussQuadrtr1',
     &       msg='Invalid number of nodes for 1D element.', ia=nNode)
        return
      end if

      end subroutine gaussQuadrtr1

! **********************************************************************

      subroutine gaussQuadrtr2(nNode,nInt,w,xi)

      ! this subroutines returns the weights and
      ! gauss point coordinate of 2D Lagrangian elements
      ! currently supports: nInt = 1, 3, 4  (for tri3 and tri4 elements)
      !                     nInt = 1, 4, 9 (for quad8 and quad8 elements)

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)     :: nNode, nInt
      real(wp), intent(out)   :: w(:), xi(:,:)
      real(wp)                :: x1D(4), w1D(4)
      type(logger)            :: msg

      w  = zero
      xi = zero

      ! 2D-PLANE: tri3 and tri6 elements
      if ( (nNode .eq. 3) .or. (nNode .eq. 6) ) then

        if (nInt .eq. 1) then
          w(1)    = half
          xi(1,1) = third
          xi(2,1) = third

        else if (nInt .eq. 3) then
          w(1:3) = sixth

          xi(1,1) = half
          xi(1,2) = half
          xi(2,1) = zero
          xi(2,2) = half
          xi(3,1) = half
          xi(3,2) = zero

        else if (nInt .eq. 4) then

          w(1) = -27.0_wp/96.0_wp
          w(2) = 25.0_wp/96.0_wp
          w(3) = w(2)
          w(4) = w(2)

          xi(1,1) = one/three
          xi(1,2) = xi(1,1)

          xi(2,1) = three/five
          xi(2,2) = one/five

          xi(3,1) = one/five
          xi(3,2) = three/five

          xi(4,1) = one/five
          xi(4,2) = one/five

        else if (nInt .eq. 6) then

          w1D(1)    = 0.0549758718227661_wp
          w1D(2)    = 0.11169079483905_wp
          x1D(1)    = 0.445948490915965_wp
          x1D(2)    = 0.091576213509771_wp

          w(1:3)    = w1D(1)
          w(4:6)    = w1D(2)
          xi(1,1)   = x1D(2)
          xi(1,2)   = x1D(2)
          xi(2,1)   = one-two*x1D(2)
          xi(2,2)   = x1D(2)
          xi(3,1)   = x1D(2)
          xi(3,2)   = one-two*x1D(2)
          xi(4,1)   = x1D(1)
          xi(4,2)   = one - two*x1D(1)
          xi(5,1)   = x1D(1)
          xi(5,2)   = x1D(1)
          xi(6,1)   = one - two*x1D(1)
          xi(6,2)   = x1D(1)

        else
         call msg%ferror(flag=error,src='gaussQuadrtr2',
     &        msg='Invalid Gauss points for TRI element.', ia=nInt)
          return
        end if

      ! 2D-PLANE: quad4 element and quad8 elements
      else if((nNode .eq. 4) .or. (nNode .eq. 8)) then

        if (nInt .eq. 1) then
          w = four

          xi(1,1) = zero
          xi(1,2) = zero

        else if (nInt .eq. 4) then

          w(1:4) = one

          x1D(1)  = sqrt(third)
          xi(1,1) = -x1D(1)
          xi(1,2) = -x1D(1)
          xi(2,1) = x1D(1)
          xi(2,2) = -x1D(1)
          xi(3,1) = -x1D(1)
          xi(3,2) = x1D(1)
          xi(4,1) = x1D(1)
          xi(4,2) = x1D(1)

        else if(nInt .eq. 9) then

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

          x1D(1)  = sqrt(three/five)
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
          call msg%ferror(flag=error,src='gaussQuadrtr2',
     &         msg='Invalid Gauss points for QUAD element.', ia=nInt)
          return
        end if

      else
        call msg%ferror(flag=error,src='gaussQuadrtr2',
     &       msg='Invalid number of nodes for 2D element.', ia=nNode)
        return
      end if

      end subroutine gaussQuadrtr2

! **********************************************************************

      subroutine gaussQuadrtr3(nNode,nInt,w,xi)

      ! this subroutines returns the weights and
      ! gauss point coordinate of 3D Lagrangian elements
      ! currently supports: nInt = 1, 4, 5  (for tet4 and tet10 elements)
      !                     nInt = 1, 8, 27 (for hex8 and hex20 elements)

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)     :: nNode, nInt
      real(wp), intent(out)   :: w(:), xi(:,:)
      real(wp)                :: x1D(4), w1D(4)
      integer                 :: i, j, k, n
      type(logger)            :: msg

      w  = zero
      xi = zero

      ! 3D-SOLID: tet4 and tet10 elements.
      if ( (nNode .eq. 4) .or. (nNode .eq. 10) )then
        if(nInt .eq. 1) then
          w(1) = sixth
          xi(1,1:3) = fourth

        else if (nInt .eq. 4) then
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

        else if (nInt .eq. 5) then
          w(1)    = -four/30.0_wp
          w(2:5)  = three/40.0_wp

          xi(1,1) = one/four
          xi(1,2) = one/four
          xi(1,3) = one/four
          xi(2,1) = half
          xi(2,2) = one/six
          xi(2,3) = one/six
          xi(3,1) = one/six
          xi(3,2) = half
          xi(3,3) = one/six
          xi(4,1) = one/six
          xi(4,2) = one/six
          xi(4,3) = half
          xi(5,1) = one/six
          xi(4,2) = one/six
          xi(5,3) = one/six

        else
          call msg%ferror(flag=error,src='gaussQuadrtr3',
     &         msg='Invalid Gauss points for TET element.', ia=nInt)
          return
        end if

      ! 3D-SOLID: hex8 and hex20 elements
      else if ((nNode .eq. 8) .or. (nNode .eq. 20)) then

        if(nInt .eq. 1) then
          w(1) = eight
          xi(1,1:3) = zero

        else if(nInt .eq. 8) then
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

        else if(nInt .eq. 27) then
          w1D(1) = five/nine
          w1D(2) = eight/nine
          w1D(3) = w1D(1)

          x1D(1) = -sqrt(0.6_wp)
          x1D(2) = zero
          x1D(3) = sqrt(0.6_wp)

          do k = 1,3
            do j = 1,3
              do i = 1,3
                n       = 9*(k-1) + 3*(j-1) + i
                w(n)    = w1D(i)*w1D(j)*w1D(k)

                xi(n,1) = x1D(i)
                xi(n,2) = x1D(j)
                xi(n,3) = x1D(k)
              end do
            end do
          end do

        else if (nInt .eq. 64) then
          w1D(1) = .3478548451374538_wp
          w1D(2) = .6521451548625461_wp
          w1D(3) = .6521451548625461_wp
          w1D(4) = .3478548451374538_wp

          x1D(1) = .8611363115940526_wp
          x1D(2) = .3399810435848563_wp
          x1D(3) = -.3399810435848563_wp
          x1D(4) = -.8611363115940526_wp

          do k = 1,4
            do j = 1,4
              do i = 1,4
                n       = 16*(k-1) + 4*(j-1) + i
                w(n)    = w1D(i)*w1D(j)*w1D(k)

                xi(n,1) = x1D(i)
                xi(n,2) = x1D(j)
                xi(n,3) = x1D(k)
              end do
            end do
          end do

        else
          call msg%ferror(flag=error,src='gaussQuadrtr3',
     &         msg='Invalid Gauss points for HEX element.', ia=nInt)
          return
        end if

      else
       call msg%ferror(flag=error,src='gaussQuadrtr3',
     &      msg='Invalid number of nodes for 3D element.', ia=nNode)
        return
      end if

      end subroutine gaussQuadrtr3

      end module gauss_quadrature

! **********************************************************************
! **********************************************************************