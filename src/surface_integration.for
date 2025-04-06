! **********************************************************************
! ************* LAGRANGE ELEMENT SURFACE INTEGRATION MODULE ************
! **********************************************************************
!   collection of subroutines to perform surface integration on 2D QUAD4
!        and 3D HEX8 subroutines are arranged in alphabetical order
! **********************************************************************
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2024. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************
      module surface_integration 

      private   :: gaussQuadrtrSurf2, gaussQuadrtrSurf3
      private   :: computeSurf2, computeSurf3

      public    :: getGaussQuadrtrSurf
      public    :: computeSurfArea

      contains
      
      subroutine getSurfGaussQuadrtr(face,w,xiSurf)

        use global_parameters,  only: wp

        implicit none

        integer, intent(in)     :: face
        real(wp), intent(out)   :: w(:), xiSurf(:,:)
        real(wp)                :: xi( size(xiSurf, dim=1) )
        real(wp)                :: eta( size(xiSurf, dim=1) )
        real(wp)                :: zeta( size(xiSurf, dim=1) )


        if (size(xiSurf, dim=2) .eq. 2) then
          call gaussQuadrtrSurf2(face,w,xi,eta)
          xiSurf(:,1)   = xi
          xiSurf(:,2)   = eta

        else if (size(xiSurf, dim=2) .eq. 3) then
          call gaussQuadrtrSurf3(face,w,xi,eta,zeta)
          xiSurf(:,1)   = xi
          xiSurf(:,2)   = eta
          xiSurf(:,3)   = zeta
        end if

      end subroutine getSurfGaussQuadrtr

! **********************************************************************

      subroutine computeSurfArea(xiIntS,face,coords,NxiS,dA)


      use global_parameters, only: wp
      
      implicit none 

      integer, intent(in)   :: face
      real(wp), intent(in)  :: xiIntS(:), coords(:,:)
      real(wp), intent(out) :: NxiS(:), dA

      if (size(xiIntS) .eq. 2) then
        call computeSurf2(xiIntS(1), xiIntS(2), face, coords, NxiS, dA)

      else if (size(xiIntS) .eq. 3) then
        call computeSurf3(xiIntS(1), xiIntS(2), xiIntS(3), 
     &                    face, coords, NxiS, dA)
      end if 


      end subroutine computeSurfArea

! **********************************************************************      
      subroutine gaussQuadrtrSurf2(face,w,xLocal,yLocal)

      ! This subroutine will get the integration point locations
      ! and corresponding gauss quadrature weights for 2D elements
      ! using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

        
      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)   :: face
      real(wp), intent(out) :: xLocal(2), yLocal(2), w(2)
      type(logger)          :: msg

      ! Gauss weights
      !
      w(1) = one
      w(2) = one
      
      ! Gauss pt locations in master element
      if(face .eq. 1) then
        xLocal(1) = -sqrt(one/three)
        yLocal(1) = -one
        xLocal(2) = sqrt(one/three)
        yLocal(2) = -one
      elseif(face  .eq.  2) then
        xLocal(1) = one
        yLocal(1) = -sqrt(one/three)
        xLocal(2) = one
        yLocal(2) = sqrt(one/three)
      elseif(face .eq. 3) then
        xLocal(1) = -sqrt(one/three)
        yLocal(1) = one
        xLocal(2) = sqrt(one/three)
        yLocal(2) = one
      elseif(face  .eq.  4) then
        xLocal(1) = -one
        yLocal(1) = sqrt(one/three)
        xLocal(2) = -one
        yLocal(2) = -sqrt(one/three)
      else
        call msg%ferror(flag=error,src='gaussQuadrtrSurf2',
     &         msg='Invalid face ID', ia=face)
      endif

      end subroutine gaussQuadrtrSurf2

!************************************************************************

      subroutine gaussQuadrtrSurf3(face,w,xLocal,yLocal,zLocal)

      ! This subroutine will get the integration point locations
      ! and corresponding gauss quadrature weights for 3D elements
      ! using 4 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  yLocal(nIntPt): z coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)   :: face
      real(wp), intent(out) :: xLocal(4), yLocal(4), zLocal(4), w(4)
      type(logger)          :: msg

      ! Gauss weights
      w(1) = one
      w(2) = one
      w(3) = one
      w(4) = one

      ! Gauss pt locations in master element
      if(face  .eq.  1) then
        xLocal(1) = -sqrt(one/three)
        yLocal(1) = -sqrt(one/three)
        zLocal(1) = -one
        xLocal(2) = sqrt(one/three)
        yLocal(2) = -sqrt(one/three)
        zLocal(2) = -one
        xLocal(3) = sqrt(one/three)
        yLocal(3) = sqrt(one/three)
        zLocal(3) = -one
        xLocal(4) = -sqrt(one/three)
        yLocal(4) = sqrt(one/three)
        zLocal(4) = -one
      elseif(face .eq. 2) then
        xLocal(1) = -sqrt(one/three)
        yLocal(1) = -sqrt(one/three)
        zLocal(1) = one
        xLocal(2) = sqrt(one/three)
        yLocal(2) = -sqrt(one/three)
        zLocal(2) = one
        xLocal(3) = sqrt(one/three)
        yLocal(3) = sqrt(one/three)
        zLocal(3) = one
        xLocal(4) = -sqrt(one/three)
        yLocal(4) = sqrt(one/three)
        zLocal(4) = one
      elseif(face .eq. 3) then
        xLocal(1) = -sqrt(one/three)
        yLocal(1) = -one
        zLocal(1) = -sqrt(one/three)
        xLocal(2) = sqrt(one/three)
        yLocal(2) = -one
        zLocal(2) = -sqrt(one/three)
        xLocal(3) = sqrt(one/three)
        yLocal(3) = -one
        zLocal(3) = sqrt(one/three)
        xLocal(4) = -sqrt(one/three)
        yLocal(4) = -one
        zLocal(4) = sqrt(one/three)
      elseif(face .eq. 4) then
        xLocal(1) = one
        yLocal(1) = -sqrt(one/three)
        zLocal(1) = -sqrt(one/three)
        xLocal(2) = one
        yLocal(2) = sqrt(one/three)
        zLocal(2) = -sqrt(one/three)
        xLocal(3) = one
        yLocal(3) = sqrt(one/three)
        zLocal(3) = sqrt(one/three)
        xLocal(4) = one
        yLocal(4) = -sqrt(one/three)
        zLocal(4) = sqrt(one/three)
      elseif(face .eq. 5) then
        xLocal(1) = -sqrt(one/three)
        yLocal(1) = one
        zLocal(1) = -sqrt(one/three)
        xLocal(2) = sqrt(one/three)
        yLocal(2) = one
        zLocal(2) = -sqrt(one/three)
        xLocal(3) = sqrt(one/three)
        yLocal(3) = one
        zLocal(3) = sqrt(one/three)
        xLocal(4) = -sqrt(one/three)
        yLocal(4) = one
        zLocal(4) = sqrt(one/three)
      elseif(face .eq. 6) then
        xLocal(1) = -one
        yLocal(1) = -sqrt(one/three)
        zLocal(1) = -sqrt(one/three)
        xLocal(2) = -one
        yLocal(2) = sqrt(one/three)
        zLocal(2) = -sqrt(one/three)
        xLocal(3) = -one
        yLocal(3) = sqrt(one/three)
        zLocal(3) = sqrt(one/three)
        xLocal(4) = -one
        yLocal(4) = -sqrt(one/three)
        zLocal(4) = sqrt(one/three)
      else
        call msg%ferror(flag=error,src='gaussQuadrtrSurf3',
     &         msg='Invalid face ID', ia=face)
          call xit
      endif

      end subroutine gaussQuadrtrSurf3


      !************************************************************************

      subroutine computeSurf2(xLocal,yLocal,face,coords,sh,ds)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the length ds, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 4-node quadrilateral elements

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)     :: face
      real(wp), intent(in)    :: xLocal, yLocal, coords(2,4)
      real(wp), intent(out)   :: ds,sh(4)
      real(wp)                :: dshxi(4,2),dXdXi,dXdEta,dYdXi
      real(wp)                :: dYdEta, normal(2,1)
      type(logger)            :: msg

      sh(1) = fourth*(one - xLocal)*(one - yLocal)
      sh(2) = fourth*(one + xLocal)*(one - yLocal)
      sh(3) = fourth*(one + xLocal)*(one + yLocal)
      sh(4) = fourth*(one - xLocal)*(one + yLocal)
      
      dshxi(1,1) = -fourth*(one - yLocal)
      dshxi(1,2) = -fourth*(one - xLocal)
      dshxi(2,1) = fourth*(one - yLocal)
      dshxi(2,2) = -fourth*(one + xLocal)
      dshxi(3,1) = fourth*(one + yLocal)
      dshxi(3,2) = fourth*(one + xLocal)
      dshxi(4,1) = -fourth*(one + yLocal)
      dshxi(4,2) = fourth*(one - xLocal)

      dXdXi = dshxi(1,1)*coords(1,1)+dshxi(2,1)*coords(1,2)
     &     + dshxi(3,1)*coords(1,3)+dshxi(4,1)*coords(1,4)
      dXdEta = dshxi(1,2)*coords(1,1)+dshxi(2,2)*coords(1,2)
     &     + dshxi(3,2)*coords(1,3)+dshxi(4,2)*coords(1,4)
      dYdXi = dshxi(1,1)*coords(2,1)+dshxi(2,1)*coords(2,2)
     &     + dshxi(3,1)*coords(2,3)+dshxi(4,1)*coords(2,4)
      dYdEta = dshxi(1,2)*coords(2,1)+dshxi(2,2)*coords(2,2)
     &     + dshxi(3,2)*coords(2,3)+dshxi(4,2)*coords(2,4)


      ! Jacobian of the mapping
      !
      if((face .eq. 2) .or. (face .eq. 4)) then
          ds = sqrt(dXdEta*dXdEta + dYdEta*dYdEta)
      elseif((face .eq. 1) .or. (face .eq. 3)) then
          ds = sqrt(dXdXi*dXdXi + dYdXi*dYdXi)
      else
          write(*,*) 'never should get here'
          call xit
      endif


      ! Surface normal, outward pointing in this case. Useful for
      ! ``follower'' type loads. The normal is referential or spatial
      ! depending on which coords were supplied to this subroutine
      ! (NOT fully tested)

      if((face .eq. 2) .or. (face .eq. 4)) then
        normal(1,1) = dYdEta/sqrt(dXdEta*dXdEta + dYdEta*dYdEta)
        normal(2,1) = -dXdEta/sqrt(dXdEta*dXdEta + dYdEta*dYdEta)
        if(face .eq. 4) normal = -normal
      elseif((face .eq. 1) .or. (face .eq. 3)) then
        normal(1,1) = dYdXi/sqrt(dXdXi*dXdXi + dYdXi*dYdXi)
        normal(2,1) = -dXdXi/sqrt(dXdXi*dXdXi + dYdXi*dYdXi)
        if(face .eq. 3) normal = -normal
      else
       call msg%ferror(flag=error,src='computeSurf2',
     &         msg='Invalid face ID', ia=face)
      endif

      return
      end subroutine computeSurf2

************************************************************************

      subroutine computeSurf3(xLocal,yLocal,zLocal,face,coords,sh,dA)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the area dA, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 8-node brick elements

      use global_parameters
      use error_logging

      implicit none

      integer, intent(in)   :: face

      real(wp), intent(in)  :: xLocal, yLocal ,zLocal, coords(3,8)
      real(wp), intent(out) :: dA, sh(8)
      real(wp)              :: dsh(8,3), dshxi(8,3), mapJ(3,3), mag
      real(wp)              :: dXdXi,dXdEta, dXdZeta, dYdXi, dYdEta
      real(wp)              :: dYdZeta, dZdXi, dZdZeta, dZdEta
      real(wp)              :: normal(3,1)
      integer               :: stat, i, j, k
      type(logger)          :: msg

      ! The shape functions
      !
      sh(1) = eighth*(one - xLocal)*(one - yLocal)*(one - zLocal)
      sh(2) = eighth*(one + xLocal)*(one - yLocal)*(one - zLocal)
      sh(3) = eighth*(one + xLocal)*(one + yLocal)*(one - zLocal)
      sh(4) = eighth*(one - xLocal)*(one + yLocal)*(one - zLocal)
      sh(5) = eighth*(one - xLocal)*(one - yLocal)*(one + zLocal)
      sh(6) = eighth*(one + xLocal)*(one - yLocal)*(one + zLocal)
      sh(7) = eighth*(one + xLocal)*(one + yLocal)*(one + zLocal)
      sh(8) = eighth*(one - xLocal)*(one + yLocal)*(one + zLocal)


      ! Shape function derivatives
      !
      dshxi(1,1) = -eighth*(one - yLocal)*(one - zLocal)
      dshxi(1,2) = -eighth*(one - xLocal)*(one - zLocal)
      dshxi(1,3) = -eighth*(one - xLocal)*(one - yLocal)
      dshxi(2,1) = eighth*(one - yLocal)*(one - zLocal)
      dshxi(2,2) = -eighth*(one + xLocal)*(one - zLocal)
      dshxi(2,3) = -eighth*(one + xLocal)*(one - yLocal)
      dshxi(3,1) = eighth*(one + yLocal)*(one - zLocal)
      dshxi(3,2) = eighth*(one + xLocal)*(one - zLocal)
      dshxi(3,3) = -eighth*(one + xLocal)*(one + yLocal)
      dshxi(4,1) = -eighth*(one + yLocal)*(one - zLocal)
      dshxi(4,2) = eighth*(one - xLocal)*(one - zLocal)
      dshxi(4,3) = -eighth*(one - xLocal)*(one + yLocal)
      dshxi(5,1) = -eighth*(one - yLocal)*(one + zLocal)
      dshxi(5,2) = -eighth*(one - xLocal)*(one + zLocal)
      dshxi(5,3) = eighth*(one - xLocal)*(one - yLocal)
      dshxi(6,1) = eighth*(one - yLocal)*(one + zLocal)
      dshxi(6,2) = -eighth*(one + xLocal)*(one + zLocal)
      dshxi(6,3) = eighth*(one + xLocal)*(one - yLocal)
      dshxi(7,1) = eighth*(one + yLocal)*(one + zLocal)
      dshxi(7,2) = eighth*(one + xLocal)*(one + zLocal)
      dshxi(7,3) = eighth*(one + xLocal)*(one + yLocal)
      dshxi(8,1) = -eighth*(one + yLocal)*(one + zLocal)
      dshxi(8,2) = eighth*(one - xLocal)*(one + zLocal)
      dshxi(8,3) = eighth*(one - xLocal)*(one + yLocal)


      dXdXi   = zero
      dXdEta  = zero
      dXdZeta = zero
      dYdXi   = zero
      dYdEta  = zero
      dYdZeta = zero
      dZdXi   = zero
      dZdEta  = zero
      dZdZeta = zero

      do k=1,8
        dXdXi = dXdXi + dshxi(k,1)*coords(1,k)
        dXdEta = dXdEta + dshxi(k,2)*coords(1,k)
        dXdZeta = dXdZeta + dshxi(k,3)*coords(1,k)
        dYdXi = dYdXi + dshxi(k,1)*coords(2,k)
        dYdEta = dYdEta + dshxi(k,2)*coords(2,k)
        dYdZeta = dYdZeta + dshxi(k,3)*coords(2,k)
        dZdXi = dZdXi + dshxi(k,1)*coords(3,k)
        dZdEta = dZdEta + dshxi(k,2)*coords(3,k)
        dZdZeta = dZdZeta + dshxi(k,3)*coords(3,k)
      enddo


      ! Jacobian of the mapping
      !
      if((face .eq. 1) .or. (face .eq. 2)) then
         ! zeta = constant on this face
         dA = sqrt(
     &          (dYdXi*dZdEta - dYdEta*dZdXi)**two
     &        + (dXdXi*dZdEta - dXdEta*dZdXi)**two
     &        + (dXdXi*dYdEta - dXdEta*dYdXi)**two
     &        )
      elseif((face .eq. 3) .or. (face .eq. 5)) then
         ! eta = constant on this face
         dA = sqrt(
     &          (dYdXi*dZdZeta - dYdZeta*dZdXi)**two
     &        + (dXdXi*dZdZeta - dXdZeta*dZdXi)**two
     &        + (dXdXi*dYdZeta - dXdZeta*dYdXi)**two
     &        )
      elseif((face .eq. 4) .or. (face .eq. 6)) then
         ! xi = constant on this face
         dA = sqrt(
     &          (dYdEta*dZdZeta - dYdZeta*dZdEta)**two
     &        + (dXdEta*dZdZeta - dXdZeta*dZdEta)**two
     &        + (dXdEta*dYdZeta - dXdZeta*dYdEta)**two
     &        )
         else
            write(*,*) 'never should get here'
            call xit
      endif


      ! Surface normal, outward pointing in this case. Useful for
      ! ``follower'' type loads. The normal is referential or spatial
      ! depending on which coords were supplied to this subroutine
      ! (NOT fully tested)

      if((face .eq. 1) .or. (face .eq. 2)) then
        ! zeta = constant on this face
        normal(1,1) = dYdXi*dZdEta - dYdEta*dZdXi
        normal(2,1) = dXdXi*dZdEta - dXdEta*dZdXi
        normal(3,1) = dXdXi*dYdEta - dXdEta*dYdXi

        if(face .eq. 1) normal = -normal

      elseif((face .eq. 3) .or. (face .eq. 5)) then
        ! eta = constant on this face
        normal(1,1) = dYdXi*dZdZeta - dYdZeta*dZdXi
        normal(2,1) = dXdXi*dZdZeta - dXdZeta*dZdXi
        normal(3,1) = dXdXi*dYdZeta - dXdZeta*dYdXi

        if(face .eq. 5) normal = -normal
      elseif((face .eq. 4) .or. (face .eq. 6)) then
        ! xi = constant on this face
        normal(1,1) = dYdEta*dZdZeta - dYdZeta*dZdEta
        normal(2,1) = dXdEta*dZdZeta - dXdZeta*dZdEta
        normal(3,1) = dXdEta*dYdZeta - dXdZeta*dYdEta
        if(face .eq. 6) normal = -normal
      else
        call msg%ferror(flag=error,src='computeSurf3',
     &         msg='Invalid face ID', ia=face)
      endif

      mag = sqrt(normal(1,1)**two+normal(2,1)**two+normal(3,1)**two)
      normal(1,1) = normal(1,1)/mag
      normal(2,1) = normal(2,1)/mag
      normal(3,1) = normal(3,1)/mag

      end subroutine computeSurf3

      end module surface_integration

************************************************************************