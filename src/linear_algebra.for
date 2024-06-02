! **********************************************************************
! *********************** LINEAR ALGEBRA MODULE ************************
! **********************************************************************
!     collection of subroutines to perform standard linear algebraic
!     computations. subroutines are arranged in alphabetical order.
! **********************************************************************
!   this module has three submodules: utilities, computational, lapack
!     * requires Intel oneMKL libraries for compiling and linking *
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2024. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************

      module linear_algebra

! **********************************************************************
!   generic interfaces for the subroutines in the utilities submodule
! **********************************************************************


      interface diag

        module function diag_m(mat) result(vec)
          use global_parameters, only: wp
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: mat(:,:)
          real(kind=wp), allocatable  :: vec(:)
        end function diag_m

        module function diag_v(vec) result(mat)
          use global_parameters, only: wp
          implicit none
          real(kind=wp), intent(in)   :: vec(:)
          real(kind=wp), allocatable  :: mat(:,:)
        end function diag_v

      end interface diag


      interface isEqual

        module function isEqualMat(A,B) result(comp)
          use global_parameters, only: wp
          implicit none
          real(kind=wp), intent(in) :: A(:,:), B(:,:)
          logical                   :: comp
        end function isEqualMat

        module function isEqualVec(a,b) result(comp)
          use global_parameters, only: wp
          implicit none
          real(kind=wp), intent(in) :: a(:), b(:)
          logical                   :: comp
        end function isEqualVec

      end interface isEqual

! **********************************************************************
!   interfaces for the other subroutines from the utilities submodule
! **********************************************************************

      interface

        module function cross_product(a,b) result(c)
          use global_parameters, only: wp
          implicit none
          real(kind=wp), intent(in)   :: a(3), b(3)
          real(kind=wp)               :: c(3)
        end function cross_product

        module subroutine crossProduct(a,b,c)
          use global_parameters, only: wp
          implicit none
          real(kind=wp), intent(in)  :: a(3), b(3)
          real(kind=wp), intent(out) :: c(3)
        end subroutine crossProduct

        module subroutine diagonal(A,v,ask)
          use global_parameters, only: wp, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(inout)  :: A(:,:)
          real(kind=wp), intent(inout)  :: v(:)
          character(len=*), intent(in)  :: ask
        end subroutine diagonal

        module function eye(m, n) result(ID)
          use global_parameters, only: wp, zero, one
          integer, intent(in)           :: m
          integer, intent(in), optional :: n
          real(kind=wp), allocatable    :: ID(:,:)
        end function eye

        module subroutine eyeMat(A)
          use global_parameters, only: wp, zero, one, error, warn
          implicit none
          real(kind=wp), intent(out)  :: A(:,:)
          integer                     :: i
        end subroutine eyeMat

        module function isSkew(A) result(skw)
          use global_parameters, only: wp
          implicit none
          real(kind=wp), intent(in)     :: A(:,:)
          logical                       :: skw
        end function isSkew

        module function isSquare(A) result(sqr)
          use global_parameters, only: wp
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          logical                     :: sqr 
        end function isSquare

        module function isSymmetric(A) result(sym)
          use global_parameters, only: wp
          implicit none
          real(kind=wp), intent(in)     :: A(:,:)
          logical                       :: sym  
        end function isSymmetric

         module function skew(A) result(skewA)
          use global_parameters, only: wp, half, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp)               :: skewA(size(A,1),size(A,1))
        end function skew

        module subroutine skewMat(A,skewA)
          use global_parameters, only: wp, half, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp), intent(out)  :: skewA(:,:)
        end subroutine skewMat

        module function sym(A) result(symA)
          use global_parameters, only: wp, half, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp)               :: symA(size(A,1),size(A,2))
        end function sym

        module subroutine symMat(A,symA)
          use global_parameters, only: wp, half, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp), intent(out)  :: symA(:,:)
        end subroutine symMat

        module function trace(A) result(trA)
          use global_parameters, only: wp, zero, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp)               :: trA
        end function trace

        module subroutine traceMat(A,trA)
          use global_parameters, only: wp, zero, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp), intent(out)  :: trA
        end subroutine traceMat

        module function triL(A) result(L)
          use global_parameters, only: wp, zero, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp)               :: L(size(A,1),size(A,2))
        end function triL

        module function triU(A) result(U)
          use global_parameters, only: wp, zero, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp)               :: U(size(A,1),size(A,2))
        end function triU

      end interface

! **********************************************************************
!     generic interface for the computational and LAPACK submodules
! **********************************************************************

      interface detMat

        module subroutine detMat_std(A,detA)
          use global_parameters, only: wp, zero, one, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp), intent(out)  :: detA
        end subroutine detMat_std

       module subroutine detMat_lapack(A,detA,lib)
          use global_parameters, only: wp, zero, one, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp), intent(out)  :: detA
          character(len=*), intent(in):: lib
        end subroutine detMat_lapack

      end interface detMat

! **********************************************************************

      interface det

        module function det_std(A) result(detA)
          use global_parameters, only: wp
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp)               :: detA
        end function det_std
          
        module function det_lapack(A,lib) result(detA)
          use global_parameters, only: wp
          implicit none 
          real(kind=wp), intent(in)   :: A(:,:)
          character(len=*), intent(in):: lib
          real(kind=wp)               :: detA
        end function det_lapack
        
      end interface det

! **********************************************************************

      interface eigen

        module subroutine eigen_lapack_sym(A,eVal,eVec,lib)
          use global_parameters, only: wp, zero, one, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)     :: A(:,:)
          real(kind=wp), intent(out)    :: eVal(:), eVec(:,:)
          character(len=*), intent(in)  :: lib
        end subroutine eigen_lapack_sym

      end interface eigen

! **********************************************************************

      interface inv

        module function inv_std(A) result(Ainv)
          use global_parameters, only: wp
          implicit none
          real(kind=wp), intent(in)     :: A(:,:)
          real(kind=wp)                 :: Ainv(size(A,1),size(A,2))
        end function inv_std 

        module function inv_lapack(A,lib) result(Ainv)
          use global_parameters, only: wp
          implicit none
          real(kind=wp), intent(in)     :: A(:,:)
          character(len=*), intent(in)  :: lib
          real(kind=wp)                 :: Ainv(size(A,1),size(A,2))
        end function inv_lapack

      end interface inv 

! **********************************************************************
      interface inverse

        module subroutine inverseMat_std(A,Ainv)
          use global_parameters, only: wp, zero, one, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp), intent(out)  :: Ainv(:,:)
        end subroutine inverseMat_std

        module subroutine inverseMat_lapack(A,Ainv,lib)
          use global_parameters, only: wp, error, warn
          implicit none
          real(kind=wp), intent(in)     :: A(:,:)
          character(len=*), intent(in)  :: lib
          real(kind=wp), intent(out)    :: Ainv(:,:)
        end subroutine inverseMat_lapack

      end interface

! **********************************************************************

      interface linSolve

        module subroutine linSolve_std_LU(A,b,x)
          use global_parameters, only: wp, zero, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(inout)  :: A(:,:), b(:)
          real(kind=wp), intent(out)    :: x(:)
        end subroutine linSolve_std_LU

        module subroutine linSolve_lapack_LU(A,b,x,lib)
          use global_parameters, only: wp
          use error_logging
          implicit none
          real(kind=wp), intent(inout)    :: A(:,:), b(:)
          character(len=*), intent(in)    :: lib
          real(kind=wp), intent(out)      :: x(:)
        end subroutine linSolve_lapack_LU

        module subroutine linSolve_lapack_QR(A,b,x,lib,method)
          use global_parameters, only: wp, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(inout)    :: A(:,:), b(:)
          character(len=*), intent(in)    :: lib, method
          real(kind=wp), intent(out)      :: x(:)
        end subroutine linSolve_lapack_QR

      end interface linSolve

! **********************************************************************

      interface LUfact

        module subroutine luDecompose_std(A,L,U)
          use global_parameters, only: wp, zero, one, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp), intent(out)  :: L(:,:), U(:,:)
        end subroutine luDecompose_std

        module subroutine luDecompose_lapack(A,P,L,U,lib)
          use global_parameters, only: wp, zero, one, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          character(len=*), intent(in):: lib
          real(kind=wp), intent(out)  :: P(size(A,1),size(A,1))
          real(kind=wp), intent(out)  :: L(size(A,1),size(A,1))
          real(kind=wp), intent(out)  :: U(size(A,1),size(A,2))
        end subroutine luDecompose_lapack

      end interface LUfact

! **********************************************************************
!     interfaces for the subroutines from the LAPACK-only submodule
! **********************************************************************

      interface

        module function cond(A,norm) result(con)
          use global_parameters, only: wp
          implicit none 
          real(kind=wp), intent(in)       :: A(:,:)
          character(len=1), intent(inout) :: norm
          real(kind=wp)                   :: con
        end function cond

        module subroutine condition(A,cond,norm)
          use global_parameters, only: wp, one, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(in)       :: A(:,:)
          real(kind=wp), intent(out)      :: cond
          character(len=1), intent(inout) :: norm
        end subroutine condition

        module function norm(A,normType) result(anorm)
          use global_parameters, only: wp
          implicit none
          real(kind=wp), intent(in)       :: A(:,:)
          character(len=1), intent(inout) :: normType
          real(kind=wp)                   :: anorm
        end function norm

        module subroutine polarDecompose(A,V,U,R)
          use global_parameters, only: wp, zero, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent (in)    :: A(:,:)
          real(kind=wp), intent (out)   :: R(:,:), U(:,:), V(:,:)
        end subroutine polarDecompose

        module function solve(A,b,lib,method) result(x)
          use global_parameters, only: wp
          use error_logging
          implicit none
          real(kind=wp), intent(in)               :: A(:,:)
          real(kind=wp), intent(in)               :: b(:)
          character(len=*), intent(in), optional  :: lib
          character(len=*), intent(in), optional  :: method
          real(kind=wp)                           :: x(size(b))
        end function solve

        module function sqrtm(A) result(sqrtA)
          use global_parameters, only: wp
          implicit none
          real(kind=wp), intent(in)   :: A(:,:)
          real(kind=wp)               :: sqrtA(size(A,1),size(A,2))  
        end function sqrtm

        module subroutine sqrtMat(A,B)
          use global_parameters, only: wp, zero, error, warn
          use error_logging
          implicit none
          real(kind=wp), intent(inout)  :: A(:,:)
          real(kind=wp), intent(out)    :: B(:,:)
        end subroutine sqrtMat

      end interface

      end module linear_algebra

! **********************************************************************
! **********************************************************************
!       submodule utilities contain functions and subroutines for
!         linear algebraic manipulations of vectors and matrices
! **********************************************************************
! **********************************************************************

      submodule (linear_algebra) utilities

      contains

! **********************************************************************
      
      module function cross_product(a,b) result(c)

        use global_parameters, only: wp

        implicit none

        real(kind=wp), intent(in)   :: a(3), b(3)
        real(kind=wp)               :: c(3)

        c(1) = a(2)*b(3)-a(3)*b(2)
        c(2) = b(1)*a(3)-a(1)*b(3)
        c(3) = a(1)*b(2)-a(2)*b(1)

      end function cross_product

! **********************************************************************

      module subroutine crossProduct(a,b,c)
      ! this subroutine computes the cross product of two vectors
      ! the input/ output arguments are 1D array

        use global_parameters, only: wp

        implicit none

        real(kind=wp), intent(in)  :: a(3), b(3)
        real(kind=wp), intent(out) :: c(3)

        c(1) = a(2)*b(3)-a(3)*b(2)
        c(2) = b(1)*a(3)-a(1)*b(3)
        c(3) = a(1)*b(2)-a(2)*b(1)

      end subroutine crossProduct

! **********************************************************************

      module function diag_m(mat) result(vec)
      ! returns the diagonal of a square matrix

        use global_parameters, only: wp, error, warn
        use error_logging

        implicit none
        real(kind=wp), intent(in)   :: mat(:,:)
        real(kind=wp), allocatable  :: vec(:)
        integer                     :: m, n, i
        type(logger)                :: msg

        m = size(mat,1)
        n = size(mat,2)

        if(m .ne. n) then
          call msg%ferror(flag=error,src='diag_m',
     &          msg ='Matrix is not square.', ivec=[m, n])
          return
        end if

        allocate(vec(m))

        do i = 1, m
          vec(i) = mat(i,i)
        end do

      end function diag_m

! **********************************************************************

      module function diag_v(vec) result(mat)
      ! sets the diagonal of a square matrix 

      use global_parameters, only: wp
      implicit none

      real(kind=wp), intent(in)   :: vec(:)
      real(kind=wp), allocatable  :: mat(:,:)
      integer                     :: dim, i

      dim = size(vec)
      allocate( mat(dim,dim) )

      do i = 1, dim
        mat(i,i) = vec(i)
      end do

      end function diag_v

! **********************************************************************

      module subroutine diagonal(A,v,ask)
      ! this subroutine returns the diagonal of a square martrix

        use global_parameters, only: wp, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(inout)  :: A(:,:)
        real(kind=wp), intent(inout)  :: v(:)
        character(len=*), intent(in)  :: ask
        integer                       :: m, n, i
        type(logger)                  :: msg

        m = size(A,1)
        n = size(A,2)

        if (m .ne. n) then
          call msg%ferror(flag=error,src='diagonal',
     &          msg ='Matrix is not square.', ivec=[m, n])
          return
        end if

        if ( trim(ask) .eq. 'get' ) then
          do i = 1, m
            v(i) = A(i,i)
          end do

        else if ( trim(ask) .eq. 'put' ) then
          do i = 1, m
            A(i,i) = v(i)
          end do

        else
          call msg%ferror(flag=error, src='diagonal',
     &            msg='Illegal argument.', ch=trim(ask))
          return
        end if

      end subroutine diagonal

! **********************************************************************

      module function eye(m, n) result(ID)

        use global_parameters, only: wp, zero, one

        integer, intent(in)           :: m
        integer, intent(in), optional :: n
        real(kind=wp), allocatable    :: ID(:,:)
        integer                       :: i, dim

        if ( present(n) ) then
          allocate(ID(m,n))
          dim = min(m,n)
        else
          allocate(ID(m,m))
          dim = m
        end if

        ID = zero
        do i = 1, dim
          ID(i, i) = one
        end do

      end function eye

  ! **********************************************************************

      module subroutine eyeMat(A)
      ! returns an identity matrix of dimension n x n

        use global_parameters, only: wp, zero, one, error, warn

        implicit none

        real(kind=wp), intent(out)  :: A(:,:)
        integer                     :: m, n, i

        m = size(A,1)
        n = size(A,2)
        A = zero

        do i = 1, min(m,n)
          A(i,i) = one
        end do

      end subroutine eyeMat

! **********************************************************************
      
      module function isEqualMat(A,B) result(comp)

        use global_parameters, only: wp

        implicit none

        real(kind=wp), intent(in) :: A(:,:), B(:,:)
        logical                   :: comp
        integer                   :: m1, n1, m2, n2

        m1  = size(A,1)
        n1  = size(A,2)
        m2  = size(B,1)
        n2  = size(B,2)

        if ( (m1 .ne. m2) .and. (n1 .ne. n2) ) then
          comp = .false.
        else
          comp = all( abs(A-B) .lt. epsilon(A) )
        end if

      end function isEqualMat

! **********************************************************************

      module function isEqualVec(a,b) result(comp)

        use global_parameters, only: wp

        implicit none

        real(kind=wp), intent(in) :: a(:), b(:)
        logical                   :: comp
        integer                   :: n1, n2

        n1  = size(a)
        n2  = size(b)

        if (n1 .ne. n2) then
          comp = .false.
        else
          comp = all( abs(a-b) .lt. epsilon(a) )
        end if

      end function isEqualVec

! **********************************************************************
      
      module function isSkew(A) result(skw)

        use global_parameters, only: wp

        implicit none

        real(kind=wp), intent(in)     :: A(:,:)
        logical                       :: skw

        if ( all( abs(A+transpose(A)) .lt. epsilon(A) ) ) then
          skw = .true.
        else 
          skw = .false.
        end if

      end function isSkew

! **********************************************************************

      module function isSquare(A) result(sqr)

        use global_parameters, only: wp

        implicit none
        
        real(kind=wp), intent(in)   :: A(:,:)
        logical                     :: sqr 
        integer                     :: m, n

        m = size(A,1)
        n = size(A,2)

        if (m .eq. n) then 
          sqr = .true.
        else 
          sqr = .false.
        end if 

      end function isSquare

! **********************************************************************
      
      module function isSymmetric(A) result(sym)

        use global_parameters, only: wp

        implicit none

        real(kind=wp), intent(in)     :: A(:,:)
        logical                       :: sym

        if ( all( abs(A-transpose(A)) .lt. epsilon(A) ) ) then
          sym = .true.
        else 
          sym = .false.
        end if

      end function isSymmetric

! **********************************************************************

      module function skew(A) result(skewA)
      ! calculates the skew part of a square matrix

        use global_parameters, only: wp, half, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        real(kind=wp)               :: skewA(size(A,1),size(A,1))
        integer                     :: m, n
        type(logger)                :: msg

        m  = size(A,1)
        n  = size(A,2)

        if ( m .ne. n ) then
          call msg%ferror(flag=error,src='skew',
     &          msg ='Matrix is not square.', ivec=[m, n])
          return
        else
          skewA = half * ( A - transpose(A) )
        end if

      end function skew

! **********************************************************************

      module subroutine skewMat(A,skewA)
      ! calculates the skew part of a square matrix

        use global_parameters, only: wp, half, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        real(kind=wp), intent(out)  :: skewA(:,:)
        integer                     :: m1, n1, m2, n2, i
        type(logger)                :: msg

        m1  = size(A,1)
        n1  = size(A,2)
        m2  = size(skewA,1)
        n2  = size(skewA,2)

        if ( (m1 .ne. m2) .and. (n1 .ne. n2) ) then
          call msg%ferror(flag=error,src='diag_m',
     &          msg ='Matrix is not square.', ivec=[m1, n1, m2, n2])
          return
        else
          skewA = half * (A - transpose(A))
        end if

      end subroutine skewMat

! **********************************************************************

      module function sym(A) result(symA)

        use global_parameters, only: wp, half, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        real(kind=wp)               :: symA(size(A,1),size(A,2))
        integer                     :: m, n
        type(logger)                :: msg

        m  = size(A,1)
        n  = size(A,2)

        if ( m .ne. n ) then
          call msg%ferror(flag=error,src='sym',
     &          msg ='Matrix is not square.', ivec=[m, n])
          return
        else
          symA = half * (A + transpose(A))
        end if

      end function sym

! **********************************************************************

      module subroutine symMat(A,symA)
      ! calculates the symmetric part of a square matix

        use global_parameters, only: wp, half, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        real(kind=wp), intent(out)  :: symA(:,:)
        integer                     :: m1, n1, m2, n2, i
        type(logger)                :: msg

        m1  = size(A,1)
        n1  = size(A,2)
        m2  = size(symA,1)
        n2  = size(symA,2)

        if ( (m1 .ne. m2) .and. (n1 .ne. n2) ) then
          call msg%ferror(flag=error,src='symMat',
     &          msg ='Matrix is not square.', ivec=[m1, n1, m2, n2])
          return
        else
          symA = half * (A + transpose(A))
        end if

      end subroutine symMat

! **********************************************************************

      module function trace(A) result(trA)

        use global_parameters, only: wp, zero, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        real(kind=wp)               :: trA
        integer                     :: i, m, n
        type(logger)                :: msg

        m = size(A,1)
        n = size(A,2)

        if (m .ne. n) then
          call msg%ferror(flag=error,src='trace',
     &          msg ='Matrix is not square.', ivec=[m, n])
          return
        end if

        trA = zero
        do i = 1, m
          trA = trA + A(i,i)
        end do

      end function trace

! **********************************************************************

      module subroutine traceMat(A,trA)

        use global_parameters, only: wp, zero, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        real(kind=wp), intent(out)  :: trA
        integer                     :: i, m, n
        type(logger)                :: msg

        m = size(A,1)
        n = size(A,2)

        if (m .ne. n) then
          call msg%ferror(flag=error,src='traceMat',
     &          msg ='Matrix is not square.', ivec=[m, n])
          return
        end if

        trA = zero
        do i = 1, m
          trA = trA + A(i,i)
        end do

      end subroutine traceMat

! **********************************************************************

      module function triL(A) result (L)
      ! this subroutine returns the lower triangular part of a square matrix

        use global_parameters, only: wp, zero, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        real(kind=wp)               :: L(size(A,1),size(A,2))
        integer                     :: i, m, n
        type(logger)                :: msg

        m = size(A,1)
        n = size(A,2)

        if (m .ne. n) then
          call msg%ferror(flag=error,src='triL',
     &          msg ='Matrix is not square.', ivec=[m , n])
          return
        end if

        L = zero
        do i = 1, m
          L(i:m,i) = A(i:m,i)
        end do

      end function triL

! **********************************************************************

      module function triU(A) result(U)
      ! this subroutine returns the upper triangular part of a square matrix

        use global_parameters, only: wp, zero, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        real(kind=wp)               :: U(size(A,1),size(A,2))
        integer                     :: i, m, n
        type(logger)                :: msg

        m = size(A,1)
        n = size(A,2)

        if (m .ne. n) then
         call msg%ferror(flag=error,src='triU',
     &          msg ='Matrix is not square.', ivec=[m , n])
          return
        end if

        U = zero
        do i = 1, m
          U(1:i,i) = A(1:i,i)
        end do

      end function triU

      end submodule utilities

! **********************************************************************
! **********************************************************************
!         computational submodule contains some subroutines for
!         linear algebraic calculations using LU decomposition
! **********************************************************************
! **********************************************************************
!   TODO: add the following subroutines to the computational submodule
!   QRdecompose_std, linSolveQR_std, eigen_std_sym, norm, condition
! **********************************************************************

      submodule (linear_algebra) computational

      contains

! **********************************************************************

      module function det_std(A) result(detA)

        use global_parameters, only: wp

        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        real(kind=wp)               :: detA

        call detMat_std(A,detA)

      end function det_std

! **********************************************************************

      module subroutine detMat_std(A,detA)
      ! this subroutine calculates the determinant of a square matrix [A]

        use global_parameters, only: wp, zero, one, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        real(kind=wp), intent(out)  :: detA
        real(kind=wp)               :: L(size(A,1),size(A,2))
        real(kind=wp)               :: U(size(A,1),size(A,2))
        integer                     :: m, n, i
        type(logger)                :: msg

        m = size(A,1)
        n = size(A,2)

        if (m .ne. n) then
          call msg%ferror(flag=error,src='detMat_std',
     &          msg ='Matrix is not square.', ivec=[m , n])
          return
        end if

        if (m .eq. 2) then
          detA = A(1,1)*A(2,2) - A(1,2)*A(2,1)

        else if (m .eq. 3) then
          detA =  A(1,1)*A(2,2)*A(3,3)
     &          + A(1,2)*A(2,3)*A(3,1)
     &          + A(1,3)*A(2,1)*A(3,2)
     &          - A(3,1)*A(2,2)*A(1,3)
     &          - A(3,2)*A(2,3)*A(1,1)
     &          - A(3,3)*A(2,1)*A(1,2)

        else if ( n .gt. 3) then

          call LUfact(A,L,U)

          detA = one
          do i = 1, m
            detA = detA*L(i,i)*U(i,i)
          end do

        end if

        if (detA .le. zero) then
          call msg%ferror(flag=warn, src='detMat_std',
     &          msg='Negative determinant.', ra=detA)
        end if

      end subroutine detMat_std

! **********************************************************************

      module function inv_std(A) result(Ainv)

        use global_parameters, only: wp

        implicit none

        real(kind=wp), intent(in)     :: A(:,:)
        real(kind=wp)                 :: Ainv(size(A,1),size(A,2))

        call inverseMat_std(A,Ainv)

      end function inv_std 

! **********************************************************************

      module subroutine inverseMat_std(A,Ainv)
      ! this subroutine returns Ainv and detA for a square matrix A

      use global_parameters, only: wp, zero, one, error, warn
      use error_logging

      implicit none

      real(kind=wp), intent(in)   :: A(:,:)
      real(kind=wp), intent(out)  :: Ainv(:,:)
      real(kind=wp)               :: detA, rdetA
      real(kind=wp)               :: L(size(A,1),size(A,2))
      real(kind=wp)               :: U(size(A,1),size(A,2))
      real(kind=wp)               :: b(size(A,1)), d(size(A,1))
      real(kind=wp)               :: x(size(A,1))
      integer                     :: m1, n1, m2, n2, i, j, k
      type(logger)                :: msg


      m1 = size(A,1)
      n1 = size(A,2)
      m2 = size(Ainv,1)
      n2 = size(Ainv,2)
      Ainv = zero

      if ( (m1 .ne. n1) .and. (m2 .ne. n2) ) then
        call msg%ferror(flag=error,src='inverseMat_std',
     &          msg ='Matrix is not square.', ivec=[m1 , n1, m2, n2])
        return
      end if

      call detMat(A,detA)

      rdetA = one/detA

      if (m1 .eq. 2) then

        Ainv(1,1) =  rdetA*A(2,2)
        Ainv(1,2) = -rdetA*A(1,2)
        Ainv(2,1) = -rdetA*A(2,1)
        Ainv(2,2) =  rdetA*A(1,1)

      else if (m1 .eq. 3) then

        Ainv(1,1) = rdetA*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
        Ainv(1,2) = rdetA*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
        Ainv(1,3) = rdetA*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
        Ainv(2,1) = rdetA*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
        Ainv(2,2) = rdetA*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
        Ainv(2,3) = rdetA*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
        Ainv(3,1) = rdetA*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
        Ainv(3,2) = rdetA*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
        Ainv(3,3) = rdetA*(A(1,1)*A(2,2)-A(2,1)*A(1,2))

      else if (m1 .gt. 3) then

        call LUfact(A,L,U)

        do k = 1,m1
          b(k)=one
          d(1) = b(1)
          do i = 2, m1
            d(i) = b(i)
            d(i) = d(i) - dot_product(L(i,1:i-1),d(1:i-1))
          end do

          x(m1)=d(m1)/U(m1,m1)
          do i = m1-1,1,-1
            x(i) = d(i)
            x(i) = x(i)-dot_product(U(i,i+1:m1),x(i+1:m1))
            x(i) = x(i)/U(i,i)
          end do

          Ainv(1:m1,k) = x(1:m1)
          b(k) = zero
        end do

      end if

      end subroutine inverseMat_std

! **********************************************************************

      module subroutine linSolve_std_LU(A,b,x)
      ! this subroutine solves for linear system using LU decomposition

        use global_parameters, only: wp, zero, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(inout)  :: A(:,:), b(:)
        real(kind=wp), intent(out)    :: x(:)
        real(kind=wp)                 :: Ainv(size(A,1),size(A,2))
        integer                       :: m, n
        type(logger)                  :: msg

        m = size(A,1)
        n = size(A,2)

        if ( m .ne. n) then
          call msg%ferror(flag=error,src='linSolve_std_LU',
     &          msg ='Matrix is not square.', ivec=[m , n])
          return
        end if

        call inverse(A,Ainv)

        x = matmul(Ainv,b)

      end subroutine linSolve_std_LU

! **********************************************************************

      module subroutine luDecompose_std(A,L,U)
      ! this subroutine computes the inverse of an arbitrary
      ! square matrix of size nxn using LU decomposition

        use global_parameters, only: wp, zero, one, error, warn
        use error_logging
        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        real(kind=wp), intent(out)  :: L(:,:), U(:,:)
        real(kind=wp)               :: coeff, detA
        integer                     :: m1, n1, m2, n2, m3, n3
        integer                     :: i, j
        type(logger)                :: msg

        m1 = size(A,1)
        n1 = size(A,2)
        m2 = size(L,1)
        n2 = size(L,2)
        m3 = size(U,1)
        n3 = size(U,2)

        if (m1 .ne. n1) then
          call msg%ferror(flag=error,src='luDecompose_std',
     &          msg ='Matrix is not square.', ivec=[m1 , n1])
          return
        end if

        !! TODO: add dimensional checking for the requested output arguments

        L = zero
        U = A

        do j = 1, m1-1
          do i = j+1, m1
            coeff = U(i,j)/U(j,j)
            L(i,j) = coeff
            U(i,j+1:m1) = U(i,j+1:m1) - coeff*U(j,j+1:m1)
          end do
        end do

        do i=1,m1
          L(i,i) = one
        end do

      end subroutine luDecompose_std

      end submodule computational

! **********************************************************************
! **********************************************************************
!    lapack submodule includes wrapper subroutines for multiple LAPACK
!     subroutines to perform linear algebraic calculation on matrices
! **********************************************************************
! **********************************************************************
!       TODO: add the following subroutines to the lapack submodule
!             normMat, QRdecompose_lapack, eigen_lapack_general
! **********************************************************************

      submodule (linear_algebra) lapack

      contains

! **********************************************************************

      module function cond(A,norm) result(con)

        use global_parameters, only: wp

        implicit none 

        real(kind=wp), intent(in)       :: A(:,:)
        character(len=1), intent(inout) :: norm
        real(kind=wp)                   :: con

        call condition(A,con,norm)

      end function cond
      
! **********************************************************************

      module subroutine condition(A,cond,norm)
      ! this subroutine calculates condiiton of a matrix based on
      ! '1' or 'I' norm. there is no subroutine in LAPACK for 2-norm

        use global_parameters, only: wp, one, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(in)       :: A(:,:)
        real(kind=wp), intent(out)      :: cond
        character(len=1), intent(inout) :: norm
        integer                         :: m, n, lda, info
        integer                         :: ipiv(min(size(A,1),
     &                                              size(A,2)))
        integer                         :: iwork(size(A,2))
        real(kind=wp)                   :: work(4*size(A,1))
        real(kind=wp)                   :: anorm, rcond, dlange
        real(kind=wp)                   :: mat(size(A,1),size(A,2))
        real(kind=wp), parameter        :: tol = 1.0e-6_wp
        type(logger)                    :: msg

        m = size(A,1)
        n = size(A,2)
        lda = m
        mat = A

        ! norm of A
        anorm = dlange(norm, m, n, mat, lda, work)

        ! LU factorization of A
        call dgetrf(m, n, mat, lda, ipiv, info)

        if (info .eq. 0) then
          call dgecon(norm, n, mat, lda, anorm, rcond, work,
     &                iwork, info)

          if (info .lt. 0) then
            call msg%ferror(flag=error, src='condition',
     &            msg='Illegal argument in DGECON.', ia=info)
            return
          end if

          cond = one/rcond

          if (rcond .lt. tol) then
            call msg%ferror(flag=warn, src='condition',
     &            msg='Ill-conditioned matrix', ra=cond)
          end if

        else if(info .gt. 0) then
          call msg%ferror(flag=error, src='condition',
     &            msg='Singular upper triangular matrix.', ia=info)

        else
          call msg%ferror(flag=error, src='condition',
     &            msg='Illegal argument in DGETRF.', ia=info)
          return
        end if

      end subroutine condition

! **********************************************************************

      module function det_lapack(A,lib) result(detA)

        use global_parameters, only: wp

        implicit none 

        real(kind=wp), intent(in)   :: A(:,:)
        character(len=*), intent(in):: lib
        real(kind=wp)               :: detA

        call detMat_lapack(A,detA,lib)

      end function det_lapack

! **********************************************************************

      module subroutine detMat_lapack(A,detA,lib)
      ! calculates determinant of matrix using LAPACK's DGETRF

        use global_parameters, only: wp, one, zero, error, warn
        use error_logging
        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        character(len=*), intent(in):: lib
        real(kind=wp), intent(out)  :: detA
        integer                     :: ipiv(min(size(A,1),size(A,2)))
        integer                     :: m, n, lda, i, info
        real(kind=wp)               :: mat( size(A,1),size(A,2) )
        type(logger)                :: msg

        if ( trim(lib) .ne. 'LAPACK' ) then
          call msg%ferror(flag=warn, src='detMat_lapack',
     &           msg='Wrong keyword for library.', ch=lib)
          return
        end if

        m = size(A,1)
        n = size(A,2)
        lda = m
        mat = A

        call dgetrf(m, n, mat, lda, ipiv, info)

        if(info .gt. 0) then
          call msg%ferror(flag=error, src='detMat_lapack',
     &            msg='Singular upper triangular matrix.', ia=info)
          return

        else if (info .lt. 0) then
         call msg%ferror(flag=error, src='detMat_lapack',
     &            msg='Illegal argument in DGETRF.', ia=info)
          return
        end if

        ! returned 'mat' contains L (with 1 in the diagonal) and U
        detA = one
        do i = 1, n
          detA = detA * mat(i,i)
          if ( ipiv(i) .ne. i ) detA = - detA
        end do

        if (detA .le. zero) then
          call msg%ferror(flag=warn, src='detMat_lapack',
     &          msg='Negative determinant.', ra=detA)
        end if

      end subroutine detMat_lapack

! **********************************************************************

      module subroutine eigen_lapack_sym(A,eVal,eVec,lib)
      ! this subroutine calculates the eigen value and vectors
      ! of a symmetric matrix using DGEEV subroutine from LAPACK

        use global_parameters, only: wp, zero, error, warn
        use error_logging

        implicit none
        real(kind=wp), intent(in)     :: A(:,:)
        real(kind=wp), intent(out)    :: eVal(:), eVec(:,:)
        character(len=*), intent(in)  :: lib
        integer                       :: info, lda, lwork, i, j
        integer                       :: dummy(1)
        integer, allocatable          :: work(:)
        integer, parameter            :: nb = 64
        integer                       :: n
        type(logger)                  :: msg


        if ( trim(lib) .ne. 'LAPACK' ) then
          call msg%ferror(flag=warn, src='detMat_lapack',
     &           msg='Wrong keyword for library.', ch=lib)
          return
        end if

        n   = size(A,1)
        lda = n

        !! TODO:  add dimensional checking for eVal and eVec
        !! TODO:  add checking for symmetry of the matrix

        !form upper triangular matrix from the input
        eVec = triU(A)

        ! first call to dsyev to query workspace
        lwork = -1
        call dsyev('Vectors', 'Upper', n, eVec, lda, eVal,
     &              dummy, lwork, info)

        ! call to dsyev to solve for eigen problem
        lwork = max((nb+2)*n, int(dummy(1)))
        allocate( work(lwork) )
        call dsyev('Vectors', 'Upper', n, eVec, lda, eVal,
     &              work, lwork, info)

        if (info .lt. 0) then
          call msg%ferror(flag=error, src='eigen_lapack_sym',
     &          msg='Illegal argument DSYEV.', ia=info)
        else if(info .gt. 0) then
          call msg%ferror(flag=error,src='eigen_lapack_sym',
     &          msg='Did not cnvergence.', ia=info)
        end if

      end subroutine eigen_lapack_sym

! **********************************************************************

      module function inv_lapack(A,lib) result(Ainv)

        use global_parameters, only: wp

        implicit none

        real(kind=wp), intent(in)     :: A(:,:)
        character(len=*), intent(in)  :: lib
        real(kind=wp)                 :: Ainv(size(A,1),size(A,2))

        call inverseMat_lapack(A,Ainv,lib)

      end function inv_lapack

! **********************************************************************

      module subroutine inverseMat_lapack(A,Ainv,lib)
      ! this subroutine calculates inverse of a square matrix
      ! it uses DGETRF and DGETRI subroutines from LAPACK

        use global_parameters, only: wp, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(in)     :: A(:,:)
        character(len=*), intent(in)  :: lib
        real(kind=wp), intent(out)    :: Ainv(:,:)
        integer                       :: ipiv(min(size(A,1),size(A,2)))
        integer                       :: m, n, lda, lwork, info
        integer, parameter            :: nb = 64
        real(kind=wp), allocatable    :: work(:)
        type(logger)                  :: msg

        if (trim(lib) .ne. 'LAPACK') then
          call msg%ferror(flag=warn, src='inverseMat_lapack',
     &           msg='Wrong keyword for library.', ch=lib)
          return
        end if


        m = size(A,1)
        n = size(A,2)
        lda = m
        lwork = nb*n
        allocate( work(lwork) )
        Ainv  = A             ! preventing A from being overwritten

        ! LAPACK subroutine: returns pivot indices from LU decomposition
        call dgetrf(m, n ,Ainv, lda, ipiv, info)

        if(info .gt. 0) then
          call msg%ferror(flag=error, src='inverseMat_lapack: dgetrf',
     &            msg='Singular upper triangular matrix.', ia=info)
          return

        else if (info .lt. 0) then
          call msg%ferror(flag=error, src='inverseMat_lapack: dgetrf',
     &            msg='Illegal argument in DGETRI.', ia=info)
          return
        end if

        ! LAPACK subroutine: calculates inverse from LU
        call dgetri(n, Ainv, lda, ipiv, work, m, info)

        if(info .gt. 0) then
          call msg%ferror(flag=error, src='inverseMat_lapack: dgetri',
     &            msg='Singular upper triangular matrix.', ia=info)
          return

        else if (info .lt. 0) then
          call msg%ferror(flag=error, src='inverseMat_lapack: dgetri',
     &            msg='Illegal argument in DGETRI.', ia=info)
          return
        end if

      end subroutine inverseMat_lapack

! **********************************************************************

      module subroutine linSolve_lapack_LU(A,b,x,lib)
      ! this subroutine solves for linear system using LAPACK subroutine

        use global_parameters, only: wp, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(inout)    :: A(:,:), b(:)
        character(len=*), intent(in)    :: lib
        real(kind=wp), intent(out)      :: x(:)
        integer                         :: ipiv(size(A,2))
        integer                         :: m, n, lda, ldb
        integer                         :: info
        real(kind=wp)                   :: mat(size(A,1),size(A,2))
        integer, parameter              :: nrhs = 1
        type(logger)                    :: msg

        if (trim(lib) .ne. 'LAPACK') then
          call msg%ferror(flag=warn, src='linSolve_lapack_LU',
     &           msg='Wrong keyword for library.', ch=lib)
          return
        end if

        m = size(A,1)
        n = size(A,2)
        lda   = m
        ldb   = size(b)     ! no of rows of the rhs vector

        if ( m .ne. n) then
          call msg%ferror(flag=error, src='linSolve_lapack_LU',
     &          msg ='Matrix is not square.', ivec=[m, n])
          return
        end if

        if (ldb .ne. m) then
          call msg%ferror(flag=error, src='linSolve_lapack_LU',
     &      msg ='Incorrect dimension of rhs vector.', ivec=[m, ldb])
          return
        end if

        ! DGESV solves the linear system and returns in b
        mat = A     ! copying to avoid overwriting
        x = b       ! copying to avoid overwriting
        call dgesv(n, nrhs, mat, lda, ipiv, x, ldb, info)

        if (info .lt. 0) then
          call msg%ferror(flag=error, src='linSolve_lapack_LU: dgesv',
     &            msg='Illegal argument in DGESV.', ia=info)
          return
        else if(info .gt. 0) then
          call msg%ferror(flag=error, src='linSolve_lapack_LU: dgesv',
     &            msg='Singular upper triangular matrix.', ia=info)
          return
        end if

      end subroutine linSolve_lapack_LU

! **********************************************************************

      module subroutine linSolve_lapack_QR(A,b,x,lib,method)
      ! this subroutine solves for linear system using LAPACK subroutine
      ! optional method is QR approach

        use global_parameters, only: wp, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(inout)    :: A(:,:), b(:)
        character(len=*), intent(in)    :: lib, method
        real(kind=wp), intent(out)      :: x(:)
        integer                         :: m, n, lda, ldb, lwork, info
        integer, allocatable            :: work(:)
        real(kind=wp)                   :: mat(size(A,1),size(A,2))
        integer, parameter              :: nb = 64, nrhs = 1
        type(logger)                    :: msg

        if ( trim(lib) .ne. 'LAPACK' ) then
          call msg%ferror(flag=warn, src='linSolve_lapack_QR',
     &           msg='Wrong keyword for library.', ch=lib)
          return
        end if

        !!!! TODO: add dimensional check for the arguments

        if ( trim(method) .eq. 'QR') then
          ! use DGELS from LAPACK to solve A*x = b using QR factorization
          m = size(A,1)         ! no of rows of the matrix
          n = size(A,2)
          lda = m
          ldb = size(b)           ! no of rows of the rhs vector
          lwork = n + nb*m

          allocate( work(lwork) )

          mat = A     ! swaping to prevent over-writing
          x = b       ! same reason swap
          call dgels('No transpose', m, n, nrhs, mat, lda, x, ldb,
     &                work, lwork, info)

        else
          call msg%ferror(flag=warn, src='linSolve_lapack_QR',
     &           msg='Wrong keyword for method.', ch=method)
          return
        end if

      end subroutine linSolve_lapack_QR

! **********************************************************************

      module subroutine luDecompose_lapack(A,P,L,U,lib)
      ! calculates A = P * L* U using LAPACK routines DGETRF

        use global_parameters, only: wp, zero, one, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        character(len=*), intent(in):: lib
        real(kind=wp), intent(out)  :: P(size(A,1),size(A,1))
        real(kind=wp), intent(out)  :: L(size(A,1),size(A,1))
        real(kind=wp), intent(out)  :: U(size(A,1),size(A,2))
        integer                     :: ipiv(min(size(A,1),size(A,2)))
        integer                     :: m, n, lda, info, i
        real(kind=wp)               :: row(size(A,1))
        type(logger)                :: msg

        if ( trim(lib) .ne. 'LAPACK' ) then
          call msg%ferror(flag=warn, src='luDecompose_lapack',
     &            msg='Wrong keyword for library.', ch=lib)
          return
        end if


        m = size(A,1)
        n = size(A,2)
        lda = m
        L = A

        call dgetrf( m, n, L, lda, ipiv, info )

        if(info .gt. 0) then
          call msg%ferror(flag=error, src='luDecompose_lapack: dgetrf',
     &            msg='Singular upper triangular matrix.', ia=info)
          return

        else if (info .lt. 0) then
          call msg%ferror(flag=error, src='luDecompose_lapack: dgetrf',
     &            msg='Illegal argument in DGETRF.', ia=info)
          return
        end if

        U = zero
        P = zero

        do i = 1, m
          U(i,i:n) = L(i,i:n)
          L(i,i:n) = zero
          L(i,i)   = one
          P(i,i)   = one
        end do

        !... Assuming that P = P[ipiv(n),n] * ... * P[ipiv(1),1]
        !... where P[i,j] is a permutation matrix for i- and j-th rows.
        do i = 1, m
          row = P(i,:)
          P(i,:) = P(ipiv(i), :)
          P(ipiv(i), :) = row
        end do

      end subroutine luDecompose_lapack

! **********************************************************************

      module subroutine polarDecompose(A,V,U,R)
      ! this subroutine compute left and right polar decompositions
      ! of matrix A (used for deformation gradient)

        use global_parameters, only: wp, zero, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent (in)    :: A(:,:)
        real(kind=wp), intent (out)   :: R(:,:), U(:,:), V(:,:)
        real(kind=wp)                 :: Vinv(size(V,1),size(V,2))
        real(kind=wp)                 :: detV
        integer                       :: m, n
        type(logger)                  :: msg

        m = size(A,1)
        n = size(A,2)

        if (m .ne. n) then
          call msg%ferror(flag=error, src='polarDecompose',
     &          msg='Matrix is not square.', ivec=[m, n])
          return
        end if

        !! TODO: add dimensional checking for the output arguments
        Vinv = zero

        !  Decompose A into A=RU=VR  where U,V are symmetric and R is orthogonal
        R = matmul(A,transpose(A))        ! R is just temporary variable here
        call sqrtMat(R,V)                 ! V= sqrt(A*A^T)
        call inverse(V,Vinv)
        R = matmul(Vinv,A)                ! R = V^-1*A
        U = matmul(transpose(R),A)        ! U = R^T*A

      end subroutine polarDecompose

! **********************************************************************

      module function norm(A,normType) result(anorm)

        use global_parameters, only: wp

        implicit none

        real(kind=wp), intent(in)       :: A(:,:)
        character(len=1), intent(inout) :: normType
        real(kind=wp)                   :: anorm
        integer                         :: m, n, lda, info
        integer                         :: ipiv(min(size(A,1),
     &                                              size(A,2)))
        real(kind=wp)                   :: work(4*size(A,1))
        real(kind=wp)                   :: dlange
        real(kind=wp)                   :: mat(size(A,1),size(A,2))

        m = size(A,1)
        n = size(A,2)
        lda = m
        mat = A

        ! norm of A
        anorm = dlange(normType, m, n, mat, lda, work)

      end function norm

! **********************************************************************

      module function solve(A,b,lib,method) result(x)

        use global_parameters, only: wp, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(in)               :: A(:,:)
        real(kind=wp), intent(in)               :: b(:)
        character(len=*), intent(in), optional  :: lib
        character(len=*), intent(in), optional  :: method
        real(kind=wp)                           :: x(size(b))
        integer                                 :: m, n
        type(logger)                            :: msg
        
        real(kind=wp)       :: mat(size(A,1),size(A,2)), vec(size(b))

        m = size(A,1)
        n = size(A,2)

        if ( m .ne. n ) then
          call msg%ferror(flag=error, src='solve',
     &          msg='Matrix is not square.', ivec=[m, n])
          return
        end if

        mat = A
        vec = b

        if ( (.not. present(lib)) .and. (.not. present(method)) ) then
          call linSolve(mat,vec,x)
        else if ( present(lib) .and. (.not. present(method)) ) then
          call linSolve(mat,vec,x,lib)
        else if ( present(lib) .and. (method .eq. 'QR') ) then 
          call linSolve(mat,vec,x,lib,method)
        else
          call msg%ferror(error, src='solve', msg='Illegal argument.')
          return
        end if

      end function solve

! **********************************************************************

      module function sqrtm(A) result(sqrtA)

        use global_parameters, only: wp

        implicit none

        real(kind=wp), intent(in)   :: A(:,:)
        real(kind=wp)               :: sqrtA(size(A,1),size(A,2))
        real(kind=wp)               :: mat(size(A,1),size(A,2))

        mat = A
        call sqrtMat(mat,sqrtA)

      end function sqrtm

! **********************************************************************

      module subroutine sqrtMat(A,B)
      ! this subroutines computes square root of a symmetric matrix

        use global_parameters, only: wp, zero, error, warn
        use error_logging

        implicit none

        real(kind=wp), intent(inout)  :: A(:,:)
        real(kind=wp), intent(out)    :: B(:,:)
        real(kind=wp)                 :: D(size(A,1),size(A,2))
        real(kind=wp)                 :: eVal(size(A,1))
        real(kind=wp)                 :: eVec(size(A,1),size(A,2))
        character(len=*), parameter   :: lib = 'LAPACK'
        integer                       :: m1, n1, m2, n2, i
        type(logger)                  :: msg

        m1 = size(A,1)
        n1 = size(A,2)
        m2 = size(B,1)
        n2 = size(B,2)

        if ( (m1 .ne. n1) .and. (m2 .ne. n2) ) then
          call  msg%ferror(flag=error, src='sqrtMat',
     &          msg='Matrix is not square.', ivec=[m1, n1, m2, n2])
          return
        end if

        D = zero
        eVal = zero
        eVec = zero

        call eigen(A,eVal,eVec,lib)

        D = zero
        do i = 1,m1
          D(i,i) = dsqrt( eVal(i) )
        end do

        B = matmul( eVec, matmul(D,transpose(eVec)) )

      end subroutine sqrtMat

      end submodule lapack

! **********************************************************************
! **********************************************************************
! **********************************************************************