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