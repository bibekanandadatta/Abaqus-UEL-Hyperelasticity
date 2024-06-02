! **********************************************************************
! ************************ ERROR LOGGING MODULE ************************
! **********************************************************************
!         Fortran class to print error and debugging messages
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2024. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************
      module error_logging

        !! global error handling variables: flags and file unit numbers
        integer, parameter  :: error  = 1,  warn      = 2
        integer, parameter  :: debug  = 3,  verbosity = debug
        integer, parameter  :: stderr = 15, stddbg    = 16
        
        type, public  :: logger
          integer, private            :: error      = error
          integer, private            :: warn       = warn
          integer, private            :: debug      = debug
          integer, private            :: verbosity  = verbosity
          integer, private            :: stderr     = stderr
          integer, private            :: stddbg     = stddbg
          character(len=512), private :: errfile    = 'error.dat'
          character(len=512), private :: dbgfile    = 'debug.dat'

          contains
          private
          procedure, public :: initialize, fopen, finfo, ferror, fdebug
        end type logger
        
! **********************************************************************

      contains

        subroutine initialize(self, verbosity, stderr, stddbg)

          implicit none

          class(logger), intent(inout)    :: self
          integer, intent(in), optional   :: verbosity
          integer, intent(in), optional   :: stderr, stddbg

          if (present(verbosity))   self%verbosity  = verbosity
          if (present(stderr))      self%stderr     = stderr
          if (present(stddbg))      self%stddbg     = stddbg

        end subroutine initialize

! **********************************************************************

        subroutine fopen(self, errfile, dbgfile)

          implicit none

          class(logger), intent(inout)            :: self
          character(len=*), intent(in), optional  :: errfile, dbgfile

          if (present(errfile))   self%errfile    = errfile
          if (present(dbgfile))   self%dbgfile    = dbgfile

          ! open the files (debug and error)
          if (self%verbosity .ge. self%debug) then
            open(unit=self%stddbg, file=self%dbgfile, status='unknown')
            open(unit=self%stderr, file=self%errfile, status='unknown')
          end if

          ! open error file only
          if ((self%verbosity .ge. self%error) .and.
     &        (self%verbosity .lt. self%debug) ) then
            open(unit=self%stderr, file=self%errfile, status='unknown')
          end if

        end subroutine fopen

! **********************************************************************

        subroutine finfo(self,msg,ch,la,ia,ra)

          use global_parameters, only: wp

          implicit none

          class(logger), intent(in)               :: self
          character(len=*), intent(in), optional  :: msg
          character(len=*), intent(in), optional  :: ch
          integer, intent(in), optional           :: ia
          real(wp), intent(in), optional          :: ra
          logical, intent(in), optional           :: la
          integer                                 :: funit

          funit = self%stderr

          if (present(msg)) then
            write(funit,'(A)',advance='no')  trim(msg)
            write(*,'(A)',advance='no')      trim(msg)
          end if

          if (present(ia)) then
            write(funit,'(1x,I0)',advance='no') ia
            write(*,'(1x,I0)',advance='no')     ia
          end if

          if (present(ra)) then
            write(funit,'(1x,g0)',advance='no')  ra
            write(*,'(1x,g0)',advance='no')      ra
          end if

          if (present(ch)) then
            write(funit,'(1x,A)',advance='no') ch
            write(*,'(1x,A)',advance='no')     ch
          end if

          if (present(la)) then
            write(funit,'(L)',advance='no')  la
            write(*,'(L)',advance='no')      la
          end if

          !! print a new line after each pass
          write(funit,'(A)',advance='no')  new_line(ch)
          write(*,'(A)',advance='no')      new_line(ch)

        end subroutine finfo

! **********************************************************************

        subroutine ferror(self,flag,src,msg,ch,la,ia,ra,ivec,rvec)

          use global_parameters, only: wp

          implicit none

          class(logger), intent(in)               :: self
          integer, intent(in)                     :: flag
          character(len=*), intent(in)            :: src
          character(len=*), intent(in), optional  :: msg
          character(len=*), intent(in), optional  :: ch
          logical, intent(in), optional           :: la
          integer, intent(in), optional           :: ia
          real(wp), intent(in), optional          :: ra
          integer, intent(in), optional           :: ivec(:)
          real(wp), intent(in), optional          :: rvec(:)
          integer                                 :: funit
          integer                                 :: values(8)
          integer                                 :: YR, MM, DD
          integer                                 :: HR, MIN, SEC

          call date_and_time(VALUES=values)
          YR    = values(1)
          MM    = values(2)
          DD    = values(3)
          HR    = values(5)
          MIN   = values(6)
          SEC   = values(7)

          if (self%verbosity .ge. self%error) then

            funit = self%stderr

            write(funit,'(6(A,I0),A)',advance='no')
     &        '[', YR,'-',MM,'-',DD, ' ', HR,':',MIN,':',SEC, '] '
            write(*,'(6(A,I0),A)',advance='no')
     &        '[', YR,'-',MM,'-',DD, ' ', HR,':',MIN,':',SEC, '] '

            if (flag .eq. self%error) then
              write(funit,'(A)', advance='no')  '(ERROR) '
              write(*,'(A)', advance='no')      '(ERROR) '
            elseif ( (flag .eq. self%warn) ) then
              write(funit,'(A)',advance='no')  '(WARNING) '
              write(*,'(A)',advance='no')      '(WARNING) '
            endif

            write(funit,'(A)',advance='no') '<'//trim(src)//'> '
            write(*,'(A)',advance='no')     '<'//trim(src)//'> '


            if (present(msg)) then
              write(funit,'(A)',advance='no')  trim(msg)
              write(*,'(A)',advance='no')      trim(msg)
            end if
  
            if (present(ia)) then
              write(funit,'(1x,I0)',advance='no') ia
              write(*,'(1x,I0)',advance='no')     ia
            end if
  
            if (present(ra)) then
              write(funit,'(1x,g0)',advance='no')  ra
              write(*,'(1x,g0)',advance='no')      ra
            end if
  
            if (present(ch)) then
              write(funit,'(1x,A)',advance='no') ch
              write(*,'(1x,A)',advance='no')     ch
            end if
  
            if (present(la)) then
              write(funit,'(L)',advance='no')  la
              write(*,'(L)',advance='no')      la
            end if

            if (present(ivec)) then
              write(funit,'(10000(2x,g0))')  ivec
              write(*,'(10000(2x,g0))')      ivec
            endif
  
            if (present(rvec)) then
              write(funit,'(10000(g0,2x))')  rvec
              write(*,'(10000(g0,2x))')      rvec
            endif
  
            !! print a new line after each pass
            write(funit,'(A)',advance='no')  new_line(ch)
            write(*,'(A)',advance='no')      new_line(ch)

          endif

        end subroutine ferror

! **********************************************************************

        subroutine fdebug(self,src,msg,ch,la,ia,ra,ivec,rvec,imat,rmat)

          use global_parameters, only: wp

          implicit none
          
          class(logger), intent(in)               :: self
          character(len=*), intent(in)            :: src
          character(len=*), intent(in), optional  :: msg
          character(len=*), intent(in), optional  :: ch
          logical, intent(in), optional           :: la
          integer, intent(in), optional           :: ia
          real(wp), intent(in), optional          :: ra
          integer, intent(in), optional           :: ivec(:)
          real(wp), intent(in), optional          :: rvec(:)
          integer, intent(in), optional           :: imat(:,:)
          real(wp), intent(in), optional          :: rmat(:,:)
          integer                                 :: funit, i


          if (self%verbosity .ge. self%debug) then

            funit = self%stddbg

            write(funit,'(A)',advance='no') '<'//trim(src)//'> '
            write(*,'(A)',advance='no')     '<'//trim(src)//'> '

            if (present(msg)) then
              write(funit,'(A)',advance='no')  trim(msg)
              write(*,'(A)',advance='no')      trim(msg)
            end if
  
            if (present(ia)) then
              write(funit,'(1x,I0)',advance='no') ia
              write(*,'(1x,I0)',advance='no')     ia
            end if
  
            if (present(ra)) then
              write(funit,'(1x,g0)',advance='no')  ra
              write(*,'(1x,g0)',advance='no')      ra
            end if
  
            if (present(ch)) then
              write(funit,'(1x,A)',advance='no') ch
              write(*,'(1x,A)',advance='no')     ch
            end if
  
            if (present(la)) then
              write(funit,'(L)',advance='no')  la
              write(*,'(L)',advance='no')      la
            end if

            if (present(ivec)) then
              write(funit,'(10000(2x,g0))')  ivec
              write(*,'(10000(2x,g0))')      ivec
            endif
  
            if (present(rvec)) then
              write(funit,'(10000(g0,2x))')  rvec
              write(*,'(10000(g0,2x))')      rvec
            endif
  
            if (present(imat)) then
              !! print a new line before matrix
              write(funit,'(A)', advance='no')  new_line(ch)
              write(*,'(A)', advance='no')      new_line(ch)
              
              do i = 1, size(imat,1)
                write(funit,'(10000(I0,2x))')  imat(i,:)
                write(*,'(10000(I0,2x))')      imat(i,:)
              enddo
            endif
  
            if (present(rmat)) then
              !! print a new line after each pass
              write(funit,'(A)', advance='no')  new_line(ch)
              write(*,'(A)', advance='no')      new_line(ch)

              do i = 1, size(rmat,1)
                write(funit,'(10000(g0,2x))')  rmat(i,:)
                write(*,'(10000(g0,2x))')      rmat(i,:)
              enddo
            endif

            !! print a new line after each pass
            write(funit,'(A)', advance='no')  new_line(ch)
            write(*,'(A)', advance='no')      new_line(ch)

          endif

        end subroutine fdebug

      end module error_logging
      
! **********************************************************************
! **********************************************************************