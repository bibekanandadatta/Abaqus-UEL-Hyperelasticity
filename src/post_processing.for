! **********************************************************************
! ******************* ABAQUS POST-PROCESSING MODULE ********************
! **********************************************************************
!         defines the mesh parameters and global variable for
!        post-processing element output of Abaqus user elements
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2024. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************
      module post_processing

      use global_parameters, only: wp

      !! post-processing related global variables for Abaqus.
      !! no of total user elements and offset for overlaying elements.
      !! numElem should be same (recommended) or higher than the
      !! total number of user elements in the overlaying element set.
      !1 elemOffset should match with the number used in input file.

      integer, parameter    :: numElem    = 50000
      integer, parameter    :: elemOffset = 100000
      real(wp), allocatable :: globalPostVars(:,:,:)

      end module post_processing
! **********************************************************************