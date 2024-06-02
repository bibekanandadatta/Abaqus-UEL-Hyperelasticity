! **********************************************************************
! ******************** GAUSSIAN QUADRATURE MODULE **********************
! **********************************************************************
! integration schemes:  (a) reduced and full integration: bar elements
!                       (b) full integration: tri and tet elements
!                       (c) reduced and full integration: quad and hex elements
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

        type(element)               :: elem
        real(kind=wp), intent(out)  :: w(:), xi(:,:)

        if (elem%nDim .eq. 1) then
          call gaussQuadrtr1(elem%nNode,elem%nInt,w,xi)
        else if (elem%nDim .eq. 2) then
          call gaussQuadrtr2(elem%nNode,elem%nInt,w,xi)
        else if (elem%nDim .eq. 3) then
          call gaussQuadrtr3(elem%nNode,elem%nInt,w,xi)
        end if

      end subroutine getGaussQuadrtr

! **********************************************************************

      subroutine gaussQuadrtr1(nNode,nInt,w,xi)

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)         :: nNode, nInt
      real(kind=wp), intent(out)  :: w(:), xi(:,:)
      type(logger)                :: msg

      w = zero
      xi = zero

      if (nNode.eq.2) then
        if (nInt.eq.1) then       ! full integration for bar2
          w = zero
          xi = two
        else
          call msg%ferror(flag=error,src='gaussQuadrtr1',
     &                    msg='Wrong Gauss points.', ia=nInt)
          return
        end if           ! end int for 2-noded bar

      else if (nNode.eq.3) then
        if (nInt .eq. 2) then     ! reduced integration for bar3
          w(1:2)  = one
          xi(1,1) = -sqrt(third)
          xi(2,1) = sqrt(third)
        else if (nInt.eq.3) then   ! full integration for bar3
          w(1) = five/nine
          w(2) = eight/nine
          w(3) = five/nine
          xi(1,1) = -sqrt(three/five)
          xi(2,1) = zero
          xi(3,1) = sqrt(three/five)
        else
          call msg%ferror(flag=error,src='gaussQuadrtr1',
     &                    msg='Wrong Gauss points.', ia=nInt)
          return
        end if           ! end int for 3-noded bar

      else
        call msg%ferror(flag=error,src='gaussQuadrtr1',
     &                  msg='Wrong Gauss points.', ia=nInt)
        return
      end if

      end subroutine gaussQuadrtr1

! **********************************************************************

      subroutine gaussQuadrtr2(nNode,nInt,w,xi)

      ! this subroutines returns the weights and
      ! gauss point coordinate of 2D Lagrangian elements
      ! currently supports: nInt = 1 (tri3) and 3 (tri6)
      !                     nInt = 1, 4 (quad4) and 4, 9 (quad8)

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)         :: nNode, nInt
      real(kind=wp), intent(out)  :: w(:), xi(:,:)
      real(kind=wp)               :: x1D(4), w1D(4)
      type(logger)                :: msg

      w  = zero
      xi = zero

      if (nNode.eq.3) then      ! plane tri3 elements (full integration)
        if (nInt.eq.1) then
          w(1)    = half
          xi(1,1) = third
          xi(2,1) = third
        else
          call msg%ferror(flag=error,src='gaussQuadrtr2',
     &                    msg='Wrong Gauss points.', ia=nInt)
          return
        end if

      else if (nNode.eq.6) then  ! plane tri6 elements (full integration)
        if (nInt.eq.3) then
          w(1:3) = sixth

          xi(1,1) = half
          xi(1,2) = half
          xi(2,1) = zero
          xi(2,2) = half
          xi(3,1) = half
          xi(3,2) = zero
        else
         call msg%ferror(flag=error,src='gaussQuadrtr2',
     &                    msg='Wrong Gauss points', ia=nInt)
          return
        end if

      else if((nNode.eq.4)) then ! plane quad4 element

        if (nInt.eq.1) then     ! reduced integration for quad4
          w = four

          xi(1,1) = zero
          xi(1,2) = zero

        else if (nInt.eq.4) then ! full integration for quad4

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
        else
          call msg%ferror(flag=error,src='gaussQuadrtr2',
     &                    msg='Wrong Gauss points', ia=nInt)
          return
        end if

      else if (nNode.eq.8) then  ! plane quad8 element

        if (nInt.eq.4) then     ! reduced integration for quad8

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

        else if(nInt.eq.9) then  ! full integration for quad8
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
     &                     msg='Wrong Gauss points', ia=nInt)
          return
        end if

      else
        call msg%ferror(flag=error,src='gaussQuadrtr2',
     &                    msg='Wrong Gauss points', ia=nInt)
        return
      end if

      end subroutine gaussQuadrtr2

! **********************************************************************

      subroutine gaussQuadrtr3(nNode,nInt,w,xi)

      ! this subroutines returns the weights and
      ! gauss point coordinate of 3D Lagrangian elements
      ! currently supports: nInt = 1 (tet4) and 4 (tet10)
      !                     nInt = 1, 8 (hex8) and 8, 27 (hex20)

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)         :: nNode, nInt
      real(kind=wp), intent(out)  :: w(:), xi(:,:)
      real(kind=wp)               :: x1D(4), w1D(4)
      integer                     :: i, j, k, n
      type(logger)                :: msg

      w  = zero
      xi = zero

      if(nNode.eq.4) then       ! 3D tet4 element (full integration)
        if(nInt.eq.1) then
          w(1) = sixth
          xi(1,1:3) = fourth

        else
          call msg%ferror(flag=error,src='gaussQuadrtr3',
     &                    msg='Wrong Gauss points', ia=nInt)
          return
        end if

      else if(nNode.eq.10) then  ! 3D tet10 element (full integration)

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
          call msg%ferror(flag=error,src='gaussQuadrtr3',
     &                    msg='Wrong Gauss points', ia=nInt)
          return
        end if

      else if(nNode.eq.8) then   ! 3D hex8 element

        if(nInt.eq.1) then      ! reduced integration for hex8
          w(1) = eight
          xi(1,1:3) = zero

        else if(nInt.eq.8) then  ! full-integration for hex8
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
         call msg%ferror(flag=error,src='gaussQuadrtr3',
     &                   msg='Wrong Gauss points', ia=nInt)
          return
        end if

      else if(nNode.eq.20) then  ! 3D hex20 element

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

        else if(nInt.eq.27) then ! full integration for hex20
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
          call msg%ferror(flag=error,src='gaussQuadrtr3',
     &                    msg='Wrong Gauss points', ia=nInt)
          return
        end if

      else
       call msg%ferror(flag=error,src='gaussQuadrtr3',
     &                 msg='Wrong Gauss points', ia=nInt)
        return
      end if

      end subroutine gaussQuadrtr3

      end module gauss_quadrature

! **********************************************************************
! **********************************************************************