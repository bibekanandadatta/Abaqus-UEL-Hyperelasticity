! **********************************************************************
! ************************ ERROR LOGGING MODULE ************************
! **********************************************************************
!         Fortran class to print error and debugging messages
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2024. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************
      module error_logging

      ! global error handling variables: flags and file unit numbers
      integer, parameter  :: debug = 1, warn = 2, error = 3   ! log levels
      integer, parameter  :: verbosity = debug                ! default verbosity
      integer, parameter  :: stderr = 15, stddbg = 16         ! file units
      
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

      ! This subroutine allows to change the verbosity of the logging 
      ! as well as to set the file unit numbers for error and debugging
      ! files. The names of the file are set while opening the file.

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
      ! this subroutine opens the debugging and error logging files
      ! user can set the name of file at the time of opening or it
      ! it will use the default name. If the files are not open
      ! then the logs will only be printed on the terminal or screen

        implicit none

        class(logger), intent(inout)            :: self
        character(len=*), intent(in), optional  :: errfile, dbgfile

        if (present(errfile))   self%errfile    = errfile
        if (present(dbgfile))   self%dbgfile    = dbgfile

        ! open files for debugging and error logging
        if (self%verbosity .ge. self%debug) then
          open(unit=self%stddbg, file=self%dbgfile, status='unknown')
          open(unit=self%stderr, file=self%errfile, status='unknown')
        else if ((self%verbosity .ge. self%warn)) then
          open(unit=self%stderr, file=self%errfile, status='unknown')
        end if

      end subroutine fopen

! **********************************************************************

      subroutine finfo(self,msg,ch,la,ia,ra,time)

      ! This subroutine is used to print simple information messages on
      ! the terminal and the file. If requeted with a logical flag,
      ! the subroutine will also print a time-stamp to the message.

        use global_parameters, only: wp

        implicit none

        class(logger), intent(in)               :: self
        character(len=*), intent(in), optional  :: msg
        character(len=*), intent(in), optional  :: ch
        logical, intent(in), optional           :: la
        integer, intent(in), optional           :: ia
        real(wp), intent(in), optional          :: ra
        logical, intent(in), optional           :: time
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

        funit = self%stderr

        if ((present(time)) .and. (time .eq. .true.)) then
          write(funit,'(6(A,I0),A)',advance='no')
     &        '[', YR,'-',MM,'-',DD, ' ', HR,':',MIN,':',SEC, '] '
          write(*,'(6(A,I0),A)',advance='no')
     &        '[', YR,'-',MM,'-',DD, ' ', HR,':',MIN,':',SEC, '] '
        end if

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

        ! print a new line after each pass
        write(funit,'(A)',advance='no')  new_line(ch)
        write(*,'(A)',advance='no')      new_line(ch)

      end subroutine finfo

! **********************************************************************

      subroutine ferror(self,flag,src,msg,ch,la,ia,ra,ivec,rvec)
      
      ! This subroutine is used to print warning and error messages on
      ! screen and in the error file (if the files are opened).
      ! It can not print matrices, and a flag and sorce file name has to
      ! be passed to the function.

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

        if (self%verbosity .lt. self%error) then

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
            write(funit,'(1x,I0,1x)',advance='no') ia
            write(*,'(1x,I0,1x)',advance='no')     ia

            ! print a new line after each pass
            write(funit,'(A)', advance='no')  new_line(ch)
            write(*,'(A)', advance='no')      new_line(ch)
          end if

          if (present(ra)) then
            write(funit,'(1x,g0,1x)',advance='no')  ra
            write(*,'(1x,g0,1x)',advance='no')      ra

            ! print a new line after each pass
            write(funit,'(A)', advance='no')  new_line(ch)
            write(*,'(A)', advance='no')      new_line(ch)
          end if

          if (present(ch)) then
            write(funit,'(1x,A,1x)',advance='no') ch
            write(*,'(1x,A,1x)',advance='no')     ch

            ! print a new line after each pass
            write(funit,'(A)', advance='no')  new_line(ch)
            write(*,'(A)', advance='no')      new_line(ch)
          end if

          if (present(la)) then
            write(funit,'(L)',advance='no')  la
            write(*,'(L)',advance='no')      la

            ! print a new line after each pass
            write(funit,'(A)', advance='no')  new_line(ch)
            write(*,'(A)', advance='no')      new_line(ch)
          end if

          if (present(ivec)) then
            write(funit,'(10000(2x,g0))')  ivec
            write(*,'(10000(2x,g0))')      ivec
          endif

          if (present(rvec)) then
            write(funit,'(10000(2x,g0))')  rvec
            write(*,'(10000(2x,g0))')      rvec
          endif

          ! print a new line after each pass
          write(funit,'(A)',advance='no')  new_line(ch)
          write(*,'(A)',advance='no')      new_line(ch)

        endif

      end subroutine ferror

! **********************************************************************

      subroutine fdebug(self,src,msg,ch,la,ia,ra,ivec,rvec,imat,rmat)

      ! This subroutine is can be used to debugging messages on the screen
      ! as well as to the debugging file (if that was opened before).
      ! All the arguments of this subroutine is optional. User can print
      ! simple messages to matrix (as of now with default formatting only).

        use global_parameters, only: wp

        implicit none
        
        class(logger), intent(in)               :: self
        character(len=*), intent(in), optional  :: src
        character(len=*), intent(in), optional  :: msg
        character(len=*), intent(in), optional  :: ch
        logical,  intent(in), optional          :: la
        integer,  intent(in), optional          :: ia
        real(wp), intent(in), optional          :: ra
        integer,  intent(in), optional          :: ivec(:)
        real(wp), intent(in), optional          :: rvec(:)
        integer,  intent(in), optional          :: imat(:,:)
        real(wp), intent(in), optional          :: rmat(:,:)
        integer                                 :: funit, i


        if (self%verbosity .le. self%debug) then

          funit = self%stddbg

          if (present(src)) then
            write(funit,'(A)',advance='no') '<'//trim(src)//'> '
            write(*,'(A)',advance='no')     '<'//trim(src)//'> '
          end if

          if (present(msg)) then
            write(funit,'(A)',advance='no')  trim(msg)
            write(*,'(A)',advance='no')      trim(msg)
          end if

          if (present(ia)) then
            write(funit,'(1x,I0,1x)',advance='no') ia
            write(*,'(1x,I0,1x)',advance='no')     ia

            ! print a new line after each pass
            write(funit,'(A)', advance='no')  new_line(ch)
            write(*,'(A)', advance='no')      new_line(ch)
          end if

          if (present(ra)) then
            write(funit,'(1x,g0,1x)',advance='no')  ra
            write(*,'(1x,g0,1x)',advance='no')      ra

            ! print a new line after each pass
            write(funit,'(A)', advance='no')  new_line(ch)
            write(*,'(A)', advance='no')      new_line(ch)
          end if

          if (present(ch)) then
            write(funit,'(1x,A,1x)',advance='no') ch
            write(*,'(1x,A,1x)',advance='no')     ch

            ! print a new line after each pass
            write(funit,'(A)', advance='no')  new_line(ch)
            write(*,'(A)', advance='no')      new_line(ch)
          end if

          if (present(la)) then
            write(funit,'(L)',advance='no')  la
            write(*,'(L)',advance='no')      la

            ! print a new line after each pass
            write(funit,'(A)', advance='no')  new_line(ch)
            write(*,'(A)', advance='no')      new_line(ch)
          end if

          if (present(ivec)) then
            write(funit,'(10000(2x,g0))')  ivec
            write(*,'(10000(2x,g0))')      ivec
          endif

          if (present(rvec)) then
            write(funit,'(10000(2x,g0))')  rvec
            write(*,'(10000(2x,g0))')      rvec
          endif

          if (present(imat)) then
            ! print a new line before matrix
            write(funit,'(A)', advance='no')  new_line(ch)
            write(*,'(A)', advance='no')      new_line(ch)
            
            do i = 1, size(imat,1)
              write(funit,'(10000(I0,2x))')  imat(i,:)
              write(*,'(10000(I0,2x))')      imat(i,:)
            enddo
          endif

          if (present(rmat)) then
            ! print a new line after each pass
            write(funit,'(A)', advance='no')  new_line(ch)
            write(*,'(A)', advance='no')      new_line(ch)

            do i = 1, size(rmat,1)
              write(funit,'(10000(g0,2x))')  rmat(i,:)
              write(*,'(10000(g0,2x))')      rmat(i,:)
            enddo
          endif

        endif

      end subroutine fdebug

      end module error_logging
      
! **********************************************************************
! **********************************************************************