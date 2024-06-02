! **********************************************************************
! *************** NONLINEAR NEWTON-RAPHSON SOLVER MODULE ***************
! **********************************************************************
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2024. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************

      module nonlinear_solver

      use global_parameters, only: wp, eps

      ! solver options for newton-raphson algorithm
      ! last two attributes are not applicable to single variable solver
      type, public  :: options
        integer           :: maxIter    = 100
        real(wp)          :: tolfx      = 1.0e-13_wp
        real(wp)          :: tolx       = 1.0e-10_wp      ! unused criterion
        character(len=8)  :: fdScheme   = 'Forward'       ! 'Backward' or 'Central'
        real(wp)          :: fdStep     = sqrt(eps)
        character(len=16) :: algo       = 'Newton'        ! 'Linesearch`
        real(wp)          :: minAlpha   = 0.1_wp
        real(wp)          :: maxAlpha   = 1.0_wp          ! start w/ NR solver
        real(wp)          :: c          = 0.5_wp
        real(wp)          :: tau        = 0.5_wp
        character(len=8)  :: lib        = 'LAPACK'        ! highly recommended
        character(len=8)  :: method     = 'QR'            ! alternative: 'LU'
      end type options

      ! as of now first two criteria are in use for convergence or exit
      ! only first 5 criterion are used in single variable solver
      ! minAlpha ... tau => are only used if 'Linesearch' is chosen
      ! alternative option for lib = 'Standard' (only 'LU' method available)
      ! alternative option for method = 'LU' or 'QR'


      ! generic interface for single variable newton-raphson solver
      interface fzero
        module procedure newton
        module procedure newton_hybrid
      end interface fzero


      ! generic interface for gradient and jacobian calculation
      interface fdjac
        module procedure dfdx       ! calculates grad of single func
        module procedure dfdx_n     ! calculates grad of multi func
      end interface fdjac


      abstract interface

        ! abstract interface for a single nonlinear equation
        subroutine func_interface(x, fx, dfx, stateVars)
          use global_parameters, only: wp
          implicit none
          real(wp), intent(in)                  :: x
          real(wp), intent(out)                 :: fx
          real(wp), intent(out), optional       :: dfx
          real(wp), intent(in), optional        :: stateVars(:)
        end subroutine func_interface

        ! abstract interface for the system of "n" nonlinear equations
        subroutine func_interface_n(x, fvec, fjac, stateVars)
          use global_parameters, only: wp
          implicit none
          real(wp), intent(in)                  :: x(:)
          real(wp), intent(out)                 :: fvec(:)
          real(wp), intent(out), optional       :: fjac(:,:)
          real(wp), intent(in), optional        :: stateVars(:)
        end subroutine func_interface_n

      end interface

! **********************************************************************

      contains

! **********************************************************************

      subroutine newton(func, xOld, x, jac, stateVars, opts)
      ! standard Newton-Raphson solver for a single nonlinear equation

        use global_parameters, only: wp
        use linear_algebra
        use error_logging

        implicit none

        procedure(func_interface)               :: func
        real(wp), intent(in)                    :: xOld
        real(wp), intent(out)                   :: x
        logical, intent(in), optional           :: jac
        real(wp), intent(in), optional          :: stateVars(:)
        type(options), intent(in), optional     :: opts
        type(options)                           :: params
        real(wp)                                :: fx, dfx, dx
        real(wp)                                :: fx0
        integer                                 :: iter
        type(logger)                            :: msg

        if ( present(opts) ) params = opts

        x = xOld

        do iter = 1, params%maxIter

          ! evaluate the function and its derivative
          if ( present(jac) .and. (jac .eq. .true.) )  then
            if ( present(stateVars) ) then
              call func(x, fx, dfx, stateVars=stateVars)
            else
              call func(x, fx, dfx)
            end if

          else if ( ( present(jac) .and. (jac .eq. .false.) )
     &              .or. (.not. present(jac)) )   then
            if ( present(stateVars) ) then
              call func(x, fx, stateVars=stateVars)
              call fdjac(func, x, fx, dfx, stateVars=stateVars)
            else
              call func(x, fx)
              call fdjac(func, x, fx, dfx)
            end if
          else
            call msg%ferror(flag=error, src='newton',
     &          msg='Illegal argument.')
            return
          end if

          if (iter .eq. 1) fx0 = fx

          if ( ( abs(fx)/abs(fx0) .le. params%tolfx )
     &        .or. ( abs(fx) .le. params%tolfx ) ) then
            return
          else
            dx  = - fx/dfx
            x   = x + dx
          end if

        end do

        call msg%ferror(flag=error, src='newton',
     &          msg='Execeeded maximum iterations.')

      end subroutine newton

! **********************************************************************

      subroutine newton_hybrid(func, xOld, x, xMin, xMax,
     &                         jac, stateVars, opts)
      ! Implements a hybrid bisection-Newton-Raphson approach
      ! follows the exit criterion of 'fail-safe Newton-Raphson' algorithm
      ! from Numerical Recipes by Press et al. (second volume)

        use global_parameters, only: wp, zero, half, two
        use error_logging

        implicit none

        procedure(func_interface)               :: func
        real(wp), intent(inout)                 :: xOld
        real(wp), intent(in)                    :: xMin, xMax
        real(wp), intent(out)                   :: x
        logical, intent(in), optional           :: jac
        real(wp), intent(in), optional          :: stateVars(:)
        type(options), intent(in), optional     :: opts
        type(options)                           :: params
        real(wp)                                :: fx, dfx
        real(wp)                                :: xl, xh, temp
        real(wp)                                :: fxMin, fxMax
        real(wp)                                :: fx0
        real(wp)                                :: dx, dxOld
        integer                                 :: iter
        type(logger)                            :: msg

        if ( present(opts) ) params = opts

        call func(xMin, fxMin, stateVars=stateVars)
        call func(xMax, fxMax, stateVars=stateVars)

        if ( fxMin*fxMax .gt. zero ) then
          x = xOld
          call msg%ferror(flag=error, src='newton_bisect',
     &     msg='Roots are not bound within limit.', rvec=[xMin, xMax])
          return
        end if

        if ( fxMin .eq. zero ) then
          x = xMin
          return
        else if ( fxMax .eq. zero ) then
          x = xMax
          return
        end if

        if ( fxMin .lt. zero ) then
          xl = xMin
          xh = xMax
        else
          xl = xMax
          xh = xMin
        end if

        if ( xOld .lt. xMin ) xOld = xMin
        if ( xOld .gt. xMax ) xOld = xMax

        x = xOld
        dxOld = abs(xMax - xMin)
        dx = dxOld

        ! evaluate the function and its derivative
        if ( present(jac) .and. (jac .eq. .true.) )  then
          if ( present(stateVars) ) then
            call func(x, fx, dfx, stateVars=stateVars)
          else
            call func(x, fx, dfx)
          end if

        else if ( ( present(jac) .and. (jac .eq. .false.) )
     &            .or. (.not. present(jac)) )   then
          if ( present(stateVars) ) then
            call func(x, fx, stateVars=stateVars)
            call fdjac(func, x, fx, dfx, stateVars=stateVars)
          else
            call func(x, fx)
            call fdjac(func, x, fx, dfx)
          end if
        else
          call msg%ferror(flag=error, src='fsolve',
     &                msg='Illegal argument.')
          return
        end if

        fx0 = fx

        do iter = 1, params%maxIter

          if ( ( (  ((x-xh)*dfx - fx) * ((x-xl)*dfx - fx) ) .gt. zero )
     &       .or. ( abs(two*fx) .gt. abs(dxOld*dfx) ) ) then

            dxOld = dx
            dx = half*(xh-xl)
            x = xl + dx

            if (xl .eq. x) return

          else
            dxOld = dx
            dx = -fx/dfx
            temp = x
            x = x + dx

            if (temp .eq. x) return
          end if

          ! evaluate the function and its derivative
          if ( present(jac) .and. (jac .eq. .true.) )  then
            if ( present(stateVars) ) then
              call func(x, fx, dfx, stateVars=stateVars)
            else
              call func(x, fx, dfx)
            end if

          else if ( ( present(jac) .and. (jac .eq. .false.) )
     &        .or. (.not. present(jac)) )   then
            if ( present(stateVars) ) then
              call func(x, fx, stateVars=stateVars)
              call fdjac(func, x, fx, dfx, stateVars=stateVars)
            else
              call func(x, fx)
              call fdjac(func, x, fx, dfx)
            end if
          else
            call msg%ferror(flag=error, src='newton_hybrid',
     &                msg='Illegal argument.')
            return
          end if

          if ( ( abs(fx)/abs(fx0) .le. params%tolfx )
     &        .or. ( abs(fx) .le. params%tolfx ) ) then
            return
          end if

          if ( fx .lt. zero ) then
            xl = x
          else
            xh = x
          end if

        end do

        call msg%ferror(flag=error, src='newton_hybrid',
     &          msg='Execeeded maximum iterations.')

      end subroutine newton_hybrid

! **********************************************************************

      subroutine fsolve(func, xOld, x, jac, stateVars, opts)
      ! Newton-Raphson solver for a system of nonlinear equations
      ! A backtracking 'Linesearch' option can be used for solution

        use global_parameters, only: wp
        use linear_algebra
        use error_logging

        implicit none

        procedure(func_interface_n)             :: func
        real(wp), intent(in)                    :: xOld(:)
        real(wp), intent(out)                   :: x(:)
        logical, optional                       :: jac
        real(wp), intent(in), optional          :: stateVars(:)
        type(options), intent(in), optional     :: opts
        type(options)                           :: params
        real(wp)                                :: fvec(size(x))
        real(wp)                                :: rhs(size(x))
        real(wp)                                :: fjac(size(x),size(x))
        real(wp)                                :: dx(size(x))
        real(wp)                                :: fvec0(size(x))
        real(wp)                                :: fnorm
        integer                                 :: iter
        type(logger)                            :: msg

        if ( present(opts) ) params = opts
        x = xOld

        do iter = 1, params%maxIter

          ! evaluate the function (n) and its jacobian (nxn)
          if ( present(jac) .and. (jac .eq. .true.) )  then
            if ( present(stateVars) ) then
              call func(x, fvec, fjac, stateVars=stateVars)
            else
              call func(x, fvec, fjac)
            end if

          else if ( ( present(jac) .and. (jac .eq. .false.) )
     &        .or. (.not. present(jac)) )   then
            if ( present(stateVars) ) then
              call func(x, fvec, stateVars=stateVars)
              call fdjac(func, x, fvec, fjac, stateVars=stateVars)
            else
              call func(x, fvec)
              call fdjac(func, x, fvec, fjac)
            end if
          else
            call msg%ferror(flag=error, src='fsolve',
     &                msg='Illegal argument.')
            return
          end if

          if (iter .eq. 1) fvec0 = fvec

          if ( ( norm2(fvec)/norm2(fvec0) .le. params%tolfx )
     &        .or. ( norm2(fvec) .le. params%tolfx ) ) then
            return
          end if

          rhs = -fvec

          if (  (trim(params%lib) .eq. 'Standard') .and.
     &          (trim(params%method) .eq. 'LU') ) then
            call linSolve(fjac, rhs, dx)
          else if (  (trim(params%lib) .eq. 'LAPACK') .and.
     &           (trim(params%method) .eq. 'LU') ) then
            call linSolve(fjac, rhs, dx, params%lib)
          else if (  (trim(params%lib) .eq. 'LAPACK') .and.
     &            (trim(params%method)) .eq. 'QR' )  then
            call linSolve(fjac, rhs, dx, params%lib, params%method)
          else
            call msg%ferror(flag=error, src='fsolve', 
     &            msg='Illegal arguments for library and method.')
            return
          end if

          if ( trim(params%algo) .eq. 'Newton' ) then
            x = x + dx

          else if ( trim(params%algo) .eq. 'Linesearch' ) then
            call lineSearch(func, x, fjac, fvec,
     &                  dx, params, x, stateVars)
          else
            call msg%ferror(flag=error, src='fsolve',
     &                msg='Illegal argument.', ch=trim(params%algo))
            return
          end if

        end do

        call msg%ferror(flag=error, src='fsolve',
     &          msg='Execeeded maximum iterations.')

      end subroutine fsolve

! **********************************************************************

      subroutine lineSearch(func, xOld, fjac, fvec,
     &                      dx, params, x, stateVars)
      ! uses a backtracking linesearch algorithm. see details below:
      ! https://en.wikipedia.org/wiki/Backtracking_line_search

        use global_parameters, only: wp, two, half

        implicit none

        procedure(func_interface_n)         :: func
        real(wp), intent(in)                :: xOld(:)
        real(wp), intent(inout)             :: fvec(:)
        real(wp), intent(in)                :: fjac(:,:)
        real(wp), intent(in)                :: dx(:)
        type(options), intent(in)           :: params
        real(wp), intent(inout)             :: x(:)
        real(wp), intent(in), optional      :: stateVars(:)
        real(wp)                            :: gradf(size(x))
        real(wp)                            :: t, slope, alpha
        logical                             :: checkAlpha
        real(wp)                            :: xtmp(size(x))
        real(wp)                            :: fvectmp(size(x))
        real(wp)                            :: fnorm0, fnorm, ftmp
        integer                             :: i

        fnorm0  = norm2(fvec)
        gradf   = matmul(fvec, fjac)
        slope   = dot_product(dx, gradf)
        t       = - params%c * slope

        alpha       = params%maxAlpha
        checkAlpha  = .false.

        do
          xtmp = xOld + alpha * dx

          if ( present(stateVars) ) then
            call func(xtmp, fvectmp, stateVars=stateVars)
          else
            call func(xtmp, fvectmp)
          end if

          ftmp = norm2(fvectmp)

          if ( ((fnorm0-ftmp)/two .ge. alpha*t) .or. checkAlpha ) then
            x     = xtmp
            fvec  = fvectmp
            fnorm = ftmp
            exit
          end if

          alpha = alpha * params%tau    ! reduce step size

          if (alpha .le. params%minAlpha) then
            alpha       = params%minAlpha
            checkAlpha  = .true.        ! will stop on the next step
          end if

        end do

      end subroutine lineSearch

! **********************************************************************
!       subroutines to calculate numerical derivative or jacobian
! **********************************************************************

      subroutine dfdx(func,x,fx,dfx,stateVars,opts)
      ! subroutine to calculate numerical derivative of a single function

        use global_parameters, only: wp, zero, two, eps
        use error_logging

        implicit none

        procedure(func_interface)               :: func
        real(wp), intent(in)                    :: x, fx
        real(wp), intent(out)                   :: dfx
        real(wp), intent(in), optional          :: stateVars(:)
        type(options), intent(inout), optional  :: opts
        type(options)                           :: params
        real(wp)                                :: h
        real(wp)                                :: x_h, fx_h
        real(wp)                                :: fx_h1, fx_h2
        integer                                 :: i
        type(logger)                            :: msg

        if ( present(opts) ) params = opts

        x_h = x

        h = abs(x_h)*params%fdStep
        if (h .eq. zero) then
          h = params%fdStep
        end if

        if ( trim(params%fdScheme).eq. 'Forward') then

          x_h = x_h + h
          call func(x_h, fx_h, stateVars=stateVars)
          dfx = (fx_h - fx)/h

        else if( trim(params%fdScheme) .eq. 'Backward') then

          x_h = x_h - h
          call func(x_h, fx_h, stateVars=stateVars)
          dfx = (fx_h - fx)/h

        else if ( trim(params%fdScheme) .eq. 'Central' ) then

          x_h = x_h + h
          call func(x_h, fx_h1, stateVars=stateVars)
          x_h = x_h - two*h
          call func(x_h, fx_h2, stateVars=stateVars)
          dfx = (fx_h1 - fx_h2)/(two*h)

        else
          call msg%ferror(flag=error, src='dfdx_n',
     &          msg='Illegal argument.', ch=trim(params%fdScheme))
          return
        end if


      end subroutine dfdx

! **********************************************************************

      subroutine dfdx_n(func,x,fvec,fjac,stateVars,opts)
      ! subroutine to calculate numerical jacobian of a set of functions

        use global_parameters, only: wp, zero, two, eps
        use error_logging

        implicit none

        procedure(func_interface_n)             :: func
        real(wp), intent(in)                    :: x(:), fvec(:)
        real(wp), intent(out)                   :: fjac(:,:)
        real(wp), intent(in), optional          :: stateVars(:)
        type(options), intent(in), optional     :: opts
        type(options)                           :: params
        real(wp)                                :: h
        real(wp)                                :: x_h(size(x))
        real(wp)                                :: fvec_h(size(x))
        real(wp)                                :: fvec_h1(size(x))
        real(wp)                                :: fvec_h2(size(x))
        integer                                 :: i
        type(logger)                            :: msg

        if ( present(opts) ) params = opts

        do i = 1, size(x)

          x_h = x

          h = abs(x_h(i))*params%fdStep
          if (h .eq. zero) then
            h = params%fdStep
          end if

          if ( trim(params%fdScheme).eq. 'Forward' ) then

            x_h(i) = x_h(i) + h
            call func(x_h, fvec_h, stateVars=stateVars)
            fjac(:,i) = (fvec_h - fvec)/h

          else if( trim(params%fdScheme) .eq. 'Backward' ) then

            x_h(i) = x_h(i) - h
            call func(x_h, fvec_h, stateVars=stateVars)
            fjac(:,i) = (fvec_h - fvec)/h

          else if ( trim(params%fdScheme) .eq. 'Central' )  then

            x_h(i) = x_h(i) + h
            call func(x_h, fvec_h1, stateVars=stateVars)
            x_h(i) = x_h(i) - two*h
            call func(x_h, fvec_h2, stateVars=stateVars)
            fjac(:,i) = (fvec_h1 - fvec_h2)/(two*h)

          else
            call msg%ferror(flag=error, src='dfdx_n',
     &          msg='Illegal argument.', ch=trim(params%fdScheme))
            return
          end if

        end do

      end subroutine dfdx_n

      end module nonlinear_solver

! **********************************************************************
! **********************************************************************