! **********************************************************************
! ********************** GLOBAL PARAMETERS MODULE **********************
! **********************************************************************
!   * defines parameters used in various functions and subroutines  *
! **********************************************************************

      module global_parameters
        
        integer, parameter    :: wp = selected_real_kind(15,307)
        
        real(wp), parameter   :: zero = 0.0_wp, one =  1.0_wp
        real(wp), parameter   :: two = 2.0_wp, three = 3.0_wp
        real(wp), parameter   :: four = 4.0_wp, five = 5.0_wp
        real(wp), parameter   :: six = 6.0_wp, seven = 7.0_wp
        real(wp), parameter   :: eight = 8.0_wp, nine = 9.0_wp

        real(wp), parameter   :: half= 0.5_wp, third = one/three
        real(wp), parameter   :: fourth = 0.25_wp, fifth = 0.2_wp
        real(wp), parameter   :: sixth = 1.0_wp/6.0_wp
        real(wp), parameter   :: eighth = 0.125_wp
        
        real(wp), parameter   :: pi = 3.1415926535897932_wp
        real(wp), parameter   :: tolx = 1.0e-10_wp
        real(wp), parameter   :: eps = epsilon(one)
      

        ! identity matrices in 2- and 3-dimensions
        real(wp)              :: ID2(2,2)
        parameter(  ID2 = reshape([ one, zero,
     &                              zero, one ], shape(ID2) ) )

        real(wp)              :: ID3(3,3)
        parameter(  ID3 = reshape([ one, zero, zero,
     &                              zero, one, zero,
     &                              zero, zero, one ], shape(ID3)) )


        !! error handling variables: flags and file unit numbers
        integer, parameter  :: info = 0, error = 1, warn = 3, debug = 3
        integer, parameter  :: stderr = 15, stddbg = 16
        
      end module global_parameters

! **********************************************************************