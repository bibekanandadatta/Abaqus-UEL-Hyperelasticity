 ! **********************************************************************
! ************************ ERROR LOGGING MODULE ************************
! **********************************************************************
!         Fortran class to print error and debugging messages
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2024. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************
      module error_logging
      implicit none

      private
      integer, parameter, public    :: info  =  1,  warn  = -1
      integer, parameter, public    :: error = -2,  fatal = -3

      type, public :: logger
        integer, private            :: stderr  = 15
        character(len=512), private :: errfile = 'error.dat'

      contains
        procedure, public           :: initialize, fopen, finfo, ferror

      end type logger

      type(logger), public, save    :: msg

! **********************************************************************

      contains

      subroutine initialize(self, stderr, errfile)

      class(logger), intent(inout)            :: self
      integer, intent(in), optional           :: stderr
      character(len=*), intent(in), optional  :: errfile

      if (present(stderr))  self%stderr  = stderr
      if (present(errfile)) self%errfile = errfile

      call self%fopen()

      end subroutine initialize


! **********************************************************************

      subroutine fopen(self, errfile)

      class(logger), intent(inout)            :: self
      character(len=*), intent(in), optional  :: errfile
      integer                                 :: ios

      if ( present(errfile) ) self%errfile = errfile

      open( unit=self%stderr, file=self%errfile,
     &      status='unknown', iostat=ios )

      if (ios .ne. 0) then
        write(*,'(A)') 'ERROR: cannot open log file.'
      end if

      end subroutine fopen

! **********************************************************************

      subroutine finfo(self, msg, ch, la, ia, ra)

      use global_parameters

      class(logger), intent(inout)            :: self
      character(len=*), intent(in), optional  :: msg, ch
      logical, intent(in), optional           :: la
      integer, intent(in), optional           :: ia
      real(wp), intent(in), optional          :: ra

      integer                                 :: funit

      funit   = self%stderr

      if (present(msg)) then
        write(*,'(A)',advance='no') trim(msg)
        write(funit,'(A)',advance='no') trim(msg)
      end if

      if (present(ia)) then
        write(*,'(1x,I0)',advance='no') ia
        write(funit,'(1x,I0)',advance='no') ia
      end if

      if (present(ra)) then
        write(*,'(1x,g0)',advance='no') ra
        write(funit,'(1x,g0)',advance='no') ra
      end if

      if (present(ch)) then
        write(*,'(1x,A)',advance='no') trim(ch)
        write(funit,'(1x,A)',advance='no') trim(ch)
      end if

      if (present(la)) then
        write(*,'(1x,L1)',advance='no') la
        write(funit,'(1x,L1)',advance='no') la
      end if

      write(*,'(A)') ''
      write(funit,'(A)') ''
      flush(funit)

      end subroutine finfo

! **********************************************************************

      subroutine ferror(self, flag, src, msg, ch, la, ia, ra,
     &                  ivec, rvec, imat, rmat)

      use global_parameters

      class(logger), intent(inout)            :: self
      integer, intent(in)                     :: flag
      character(len=*), intent(in)            :: src
      character(len=*), intent(in), optional  :: msg, ch
      logical, intent(in), optional           :: la
      integer, intent(in), optional           :: ia
      real(wp), intent(in), optional          :: ra
      integer, intent(in), optional           :: ivec(:)
      real(wp), intent(in), optional          :: rvec(:)
      integer, intent(in), optional           :: imat(:,:)
      real(wp), intent(in), optional          :: rmat(:,:)

      integer                                 :: funit, values(8)
      integer                                 :: YR, MM, DD, HR, MN, SC
      integer                                 :: i

      funit   = self%stderr

      ! timestamp is printed for warn/error/fatal only
      if (flag .ne. info) then
        call date_and_time(values=values)
        YR = values(1); MM = values(2); DD = values(3)
        HR = values(5); MN = values(6); SC = values(7)

        write(*,'(6(A,I0),A)',advance='no')
     &        '[',YR,'-',MM,'-',DD,' ',HR,':',MN,':',SC,'] '
        write(funit,'(6(A,I0),A)',advance='no')
     &        '[',YR,'-',MM,'-',DD,' ',HR,':',MN,':',SC,'] '
      end if

      ! prefix by flag level
      select case (flag)
      case (info)
        write(*,'(A)',advance='no') '(INFO) '
        write(funit,'(A)',advance='no') '(INFO) '
      case (warn)
        write(*,'(A)',advance='no') '(WARNING) '
        write(funit,'(A)',advance='no') '(WARNING) '
      case (error, fatal)
        write(*,'(A)',advance='no') '(ERROR) '
        write(funit,'(A)',advance='no') '(ERROR) '
      case default

      end select

      if (len_trim(src) .gt. 0) then
        write(*,'(A)',advance='no') '<'//trim(src)//'> '
        write(funit,'(A)',advance='no') '<'//trim(src)//'> '
      end if

      if (present(msg)) then
        write(*,'(A)',advance='no') trim(msg)
        write(funit,'(A)',advance='no') trim(msg)
      end if

      if (present(ia)) then
        write(*,'(1x,I0)',advance='no') ia
        write(funit,'(1x,I0)',advance='no') ia
      end if

      if (present(ra)) then
        write(*,'(1x,g0)',advance='no') ra
        write(funit,'(1x,g0)',advance='no') ra
      end if

      if (present(ch)) then
        write(*,'(1x,A)',advance='no') trim(ch)
        write(funit,'(1x,A)',advance='no') trim(ch)
      end if

      if (present(la)) then
        write(*,'(1x,L1)',advance='no') la
        write(funit,'(1x,L1)',advance='no') la
      end if

      write(*,'(A)') ''
      write(funit,'(A)') ''

      if (present(ivec)) then
        write(*,'(10000(2x,g0))') ivec
        write(funit,'(10000(2x,g0))') ivec
      end if

      if (present(rvec)) then
        write(*,'(10000(2x,g0))') rvec
        write(funit,'(10000(2x,g0))') rvec
      end if

      if (present(imat)) then
        do i = 1, size(imat,1)
          write(*,'(10000(I0,2x))') imat(i,:)
          write(funit,'(10000(I0,2x))') imat(i,:)
        end do
      end if

      if (present(rmat)) then
        do i = 1, size(rmat,1)
          write(*,'(10000(g0,2x))') rmat(i,:)
          write(funit,'(10000(g0,2x))') rmat(i,:)
        end do
      end if

      if ((flag .eq. error) .or. (flag .eq. fatal)) flush(funit)

      end subroutine ferror

! **********************************************************************

      end module error_logging
! **********************************************************************
! **********************************************************************