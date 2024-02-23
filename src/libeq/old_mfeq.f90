! Compile with $ f2py3 -llapack -c -m mfeq eq.f90

subroutine lsolve(A, b, x, n)
  ! Solve system of lineal equations using LAPACK

  implicit none

  integer, intent(in) :: n
  double precision, dimension(n,n), intent(in) :: A
  ! The nxn array of coefficients
  double precision, dimension(n), intent(in) :: b
  ! the n array of independent terms
  double precision, dimension(n), intent(out) :: x
  ! the n array where the result is returned

  integer, dimension(n) :: ipiv   ! pivot indices
  integer :: info, nrhs, lda, ldb

  ! External procedures defined in LAPACK
  external dgesv

  ! Store b in x to prevent it from being overwritten by LAPACK
  x = b

  ! Set up solve by LAPACK routine DGESV
  nrhs = 1 ! number of right hand sides in b
  lda = n  ! leading dimension of a
  ldb = n  ! leading dimension of b

  call dgesv(n, nrhs, A, lda, ipiv, x, ldb, info)

  if (info /= 0) then
     stop 'Linear system could not be solved!'
  end if
end subroutine lsolve


subroutine newtonr1(C, T, B, P, CTRL, N, E, S)
  ! Perform Newton-Raphson iterations per point when the beta matrix
  ! is 1D. The first initial value must be provided in the first row
  ! of C. The following points take as starting value the value of 
  ! the last iteration.

  implicit none
  integer, intent(in) :: N, E, S
  ! N -> number of experimental points
  ! E -> number of equilibria
  ! S -> number of species
  integer, intent(inout) :: CTRL
  ! control value. 0 on no error, 1 otherwise
  integer, dimension(E, S), intent(in) :: P
  ! The stoichiometry array
  double precision, dimension(N, E+S), intent(inout) :: C
  ! C -> (inout) initial estimates and refined values
  double precision, dimension(N, S), intent(in) :: T
  ! T -> The stoichiometry values
  double precision, dimension(E), intent(in) :: B
  integer i, j, ctrl2
  ctrl2 = ctrl

  do i=1,N
    call nr(C(i,:), T(i,:), B, P, CTRL, E, S)
    if (ctrl .eq. 0) then
      do j=1,S
        if (i .le. 3) then      ! Use previous point as new guess
          C(i+1, j) = C(i, j)   ! for the next except for the last one
        else
          call interp(C(i-3:i,:), S)
        end if
      end do
    end if
  end do
end subroutine newtonr1


subroutine newtonr2(C, T, B, P, CTRL, N, E, S)
  ! Perform Newton-Raphson iterations per point when the beta matrix
  ! is 2D. The first initial value must be provided in the first row
  ! of C. The following points take as starting value the value of 
  ! the last iteration.

  implicit none
  integer, intent(inout) :: CTRL
  ! CTRL -> (out) control value
  integer, dimension(E, S), intent(in) :: P
  double precision, dimension(N, E+S), intent(inout) :: C
  ! C -> (inout) initial estimates
  double precision, dimension(N, S), intent(in) :: T
  ! T -> The stoichiometry values
  double precision, dimension(N, E), intent(in) :: B
  integer, intent(in) :: N, E, S
  ! N -> (in) number of experimental points
  ! E -> (in) number of equilibria
  ! S -> (in) number of species
  integer i, j, ctrl2
  ctrl2 = CTRL

  write (*,*) N, E, S, shape(C), shape(T), shape(B)
  do i=1,N
    call nr(C(i,:), T(i,:), B(i,:), P, CTRL, E, S)
    if (ctrl .eq. 0) then
      do j=1,S
        if (i .le. 3) then      ! Use previous point as new guess
          C(i+1, j) = C(i, j)   ! for the next
        else
          call interp(C(i-3:i,:), S)
        end if
      end do
    end if
  end do
end subroutine newtonr2


subroutine cplussn(C, B, P, N, E, S)
  ! Compute the free concentrations for all the components for all
  !  titration points.
  integer N, E, S
  double precision C(N, S+E)
  integer P(E, S)
  double precision B(E)
  intent(in) N, E, S, B, P
  intent(inout) C
  integer i
  DO i=1,N
    call cpluss(C(i,:), B, P, E, S)
  END DO
END subroutine cplussn


subroutine cpluss(C, B, P, E, S)
  ! Compute the free concentrations for all the components.
  integer, intent(in) :: E
  integer, intent(in) :: S
  double precision, dimension(S+E), intent(inout) :: C
  integer, dimension(E, S), intent(in):: P
  double precision, dimension(E) :: B
  double precision AUX
  integer i, j
  do i=1,E
    AUX = 1D0
    do j=1,S
      AUX = AUX*C(j)**P(i,j)
    end do
    C(i+S) = B(i)*AUX
  end do
end subroutine cpluss


subroutine fobj(F, C, P, T, E, S)
  ! Compute the objective function
  ! f_i = c_i + \sum_{j=1}^E { p_{ji}c_{j+S} } - T_i
  INTEGER, intent(in) :: E, S
  double precision, dimension(S), intent(inout) :: F
  double precision, dimension(S), intent(in) :: T
  integer, dimension(E, S), intent(in):: P
  double precision, dimension(S+E), intent(in) :: C
  double precision AUX
  integer i, j
  DO i=1,S
    AUX = C(i)-T(i)
    DO j=1,E
      AUX = AUX+P(j,i)*C(j+S)
    END DO
    F(i)=AUX
  END DO
END subroutine fobj


subroutine jacobian(Jac, C, P, E, S)
  ! Compute the free concentrations for all the components.
  ! J_{ij} = \delta_{ij} + c_j^{-1} \sum_{k=1}^E {
  !        p_{ki}p_{kj}c_{k+S} }
  integer, intent(in) :: E, S
  integer, dimension(E, S), intent(in) :: P
  double precision, dimension(S+ E), intent(in) :: C
  double precision, dimension(S, S), intent(inout) :: Jac
  integer i, j, k
  double precision aux
  do i=1,S
    do j=1,S
      aux = 0.0
      do k=1,E
        aux=aux+P(k,i)*P(k,j)*C(k+s)
      end do
      Jac(i,j)=aux/C(j)
      if (j .eq. i) Jac(i,j) = 1.0 + Jac(i,j)
    end do
  end do
end subroutine jacobian


subroutine nr(X, T, B, P, ctrl, E, S)
  ! Perform the actual Newton-Raphson iterations. The initial value is supplied
  ! in the X variable and it is replaced with the refined value upon exit.
  implicit none
  integer, intent(in) :: E, S
  integer, intent(inout) :: ctrl
  ! ctrl -> control value. 0 = normal exit on TOLF, 1 = normal exit on TOLX,
  !         2 = normal exit upon initial guess being a root
  !         99 = abnormal exit on MAXITS
  integer, dimension(E, S), intent(in) :: P
  double precision, dimension(S+E), intent(inout) :: X
  double precision, dimension(S), intent(in) :: T
  double precision, dimension(E), intent(in) ::  B

  integer its, i, j
  double precision, dimension(S) :: fvec, jdiag
  double precision, dimension(S, S) ::  jac
  double precision :: f, fold
  double precision, dimension(S) :: xold
  double precision, dimension(S) :: g
  double precision, dimension(S) :: dx
  double precision sump, test, temp
  logical check
  integer MAXITS
  double precision TOLX, TOLF
  parameter (MAXITS=1000, TOLF=1.d-20, TOLX=1.d-15)
  
  f = 0.0
  do i=1,S
    xold(i) = X(i)
  end do
  call cpluss(X, B, P, E, S)
  call fobj(fvec, X, P, T, E, S)

  ! Test for initial guess being a root. Use more stringent test than simply TOLF.
  test=0.
  do i=1,S
    if (abs(fvec(i)) .gt. test) test=abs(fvec(i))
  end do
  if (test .lt. .01*TOLF) then
    check=.false.
    ctrl=2
    return
  end if

  ! initial damping is helpful
  !call damping(X, fvec, T, B, P, E, S)
  call pcf(X, T, B, P, E, S)

  ! Main loop
  do its=1,MAXITS
    call jacobian(jac, X, P, E, S)

    ! compute gradient for lnsearch
    do i=1,S
      sump=0.
      do j=1,S
        sump=sump+jac(j,i)*fvec(j)
      end do
      g(i)=sump
    end do

    ! Store x and f into old
    fold = 0.
    do i=1,S
      xold(i)=x(i)
      fold = fold + fvec(i)*fvec(i)
    end do
    fold=0.5*fold

    ! Scale
    do i=1,S
      jdiag(i) = jac(i,i)
    end do
    do i=1,S
      do j=1,S
        jac(i, j) = jac(i, j) / sqrt(jdiag(i)*jdiag(j))
      end do
      fvec(i) = fvec(i) / sqrt(jdiag(i))
    end do

    ! Solve linear system. Includes negation of fvec
    call lsolve(jac, fvec, dx, S)

    ! De-scale to obtain the actual 'dx'
    do i=1,S
      dx(i) = - dx(i) / sqrt(jdiag(i)) 
    end do

    ! lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.
    ! output                    v    v   v  v             v
    ! call lnsrch(xold, fold, g, dx, fvec, x, f, B, T, P, check, S, E)

    ! Test for convergence on function values.
    test=0.
    do i=1,S
      if (abs(fvec(i)) .gt. test) test = abs(fvec(i))
    end do
    if (test .lt. TOLF) then
      check=.false.
      ctrl=0
      ! write (*,*) "Exit on TOLF"
      return
    endif

    test=0.
    ! Test for convergence on delta_x.
    do i=1,S
      !temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
      temp = (abs(dx(i))/abs(x(i)))
      if (temp .gt. test) test = temp
    end do
    if (test .lt. TOLX) then
      ! write (*,*) "Exit on TOLX"
      ctrl = 1
      return
    end if
    call damping(X, fvec, T, B, P, E, S)
  end do
  ctrl=99
  ! stop ’MAXITS exceeded in newt’
  ! write (*,*) "Exit on MAXITS"
end subroutine nr


subroutine lnsrch(xold, fold, g, dx, fvec, x, f, B, T, P, check, S, E)
  integer, intent(in) :: S, E
  double precision, dimension(S), intent(in) :: xold, g
  double precision, dimension(S+E), intent(inout) :: x
  double precision, dimension(S), intent(inout) :: dx
  double precision, dimension(S), intent(inout) :: fvec
  double precision, intent(in) :: fold
  double precision, intent(inout) :: f      ! f = 0.5*F**2
  double precision, dimension(E), intent(in) :: B
  double precision, dimension(S), intent(in) :: T
  integer, dimension(E, S), intent(in):: P
  logical, intent(inout) :: check
  ! Given an n-dimensional point xold(1:n), the value of the function and gradient there,
  ! fold and g(1:n), and a direction dx(1:n), finds a new point x(1:n) along the direction
  ! p from xold where the function func has decreased sufficiently. The new function value
  ! is returned in f. stpmax is an input quantity that limits the length of the steps so that you
  ! do not try to evaluate the function in regions where it is  unde ned  or subject to overflow.
  ! p is usually the Newton direction. The output quantity check
  ! is false on a normal exit. It is true  when x is  too  close  to xold.
  ! In a minimization  algorithm, this usually signals
  ! convergence and can be  ignored.  However, in  a zero-finding algorithm the calling  program
  ! should check  whether  the  convergence  is  spurious.
  ! Parameters:
  !  ALF       ensures  sufficient  decrease  in  function  value;
  !  TOLX is  the  convergence criterion  on Delta x .
  double precision ALF, TOLX
  parameter (ALF=1.e-4,TOLX=1.e-15)
  integer i
  double precision aux_a,alam,alam2,alamin,aux_b,disc,f2,rhs1,rhs2,slope,temp,test,tmplam,stpmax

  check=.false.
  f2 = 0.0
  alam2 = 0.0

  ! Scale if necessary to forbid negative concentrations
  stpmax=1.
  do i=1,S
    if (xold(i) + dx(i) .lt. 0.0) then
      stpmax = min(-xold(i) / dx(i), stpmax)
    end if
  enddo
  do i=1,S
    dx(i)=dx(i)*stpmax/(1.0+TOLX)
  enddo

  slope=0.
  do i=1,S
    slope=slope+g(i)*dx(i)
  enddo

  if(slope.ge.0.) stop 'roundoff problem in lnsrch'

  test=0.
  ! Compute lambda min.
  do i=1,S
    temp=abs(dx(i))/max(abs(xold(i)),1.)
    if (temp.gt.test) test=temp
  enddo
  alamin=TOLX/test

  alam=1.  ! Always try full Newton step first.
  do
    ! Start of iteration loop.
    do i=1,S
      x(i)=xold(i)+alam*dx(i)
    enddo

    ! evaluate f
    call cpluss(x, B, P, E, S)
    call fobj(fvec, x, P, T, E, S)
    f = 0.0
    do i=1,S
      f = f + fvec(i)*fvec(i)
    end do
    f = 0.5*f
    ! write (*,*) alam, f

    if (alam .lt. alamin) then
      ! Convergence  on Delta x. For  zero finding, the calling program
      ! should verify the convergence.
      do i=1,S
        x(i)=xold(i)
      enddo

      check=.true.
      return
    else if(f .le. fold+ALF*alam*slope) then
      ! Sufficient function decrease.
      return
    else
      ! Backtrack.
      if(alam.eq.1.)then
        ! First  time.
        tmplam=-slope/(2.*(f-fold-slope))
        f2=0.0
      else
        ! Subsequent backtracks.
        rhs1=f-fold-alam*slope
        rhs2=f2-fold-alam2*slope
        aux_a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
        aux_b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
        if(aux_a.eq.0.)then
          tmplam=-slope/(2.*aux_b)
        else
          disc=aux_b*aux_b-3.*aux_a*slope
          if(disc.lt.0.)then
            tmplam=.5*alam
          else if(aux_b.le.0.)then
            tmplam=(-aux_b+sqrt(disc))/(3.*aux_a)
          else
            tmplam=-slope/(aux_b+sqrt(disc))
          endif
        endif
        if(tmplam.gt..5*alam)tmplam=.5*alam  ! lambda <= 0.5 lambda_1
      endif
    endif
    alam2=alam
    f2=f
    alam=max(tmplam,.1*alam) ! lambda >= 0.1 lambda_1
  ! Try  again.
  end do
end


subroutine damping(C, FVEC, T, B, P, E, S)
  integer, intent(in) :: S, E
  double precision, dimension(S), intent(inout) :: FVEC
  double precision, dimension(S+E), intent(inout) :: C
  double precision, dimension(S), intent(in) :: T
  double precision, dimension(E), intent(in) :: B
  integer, dimension(E, S), intent(in):: P

  double precision R
  integer i
  double precision g
  logical guard
  integer dampcount, MAXDAMPS
  parameter (MAXDAMPS=20)

  do dampcount=1,MAXDAMPS
    guard = .true.
    g = .5
    do i=1,S
      if (T(i) + FVEC(i) .gt. 0.0) then
        R = T(i) / (T(i) + FVEC(i))
        if (R .gt. 0.0) then 
          if (R < 0.1 .or. R > 10.0) then
            guard = .false.
            C(i) = C(i) * R**g
          end if
        end if
      end if
    end do
    if (guard) return
    call cpluss(C, B, P, E, S)
    call fobj(FVEC, C, P, T, E, S)
  end do
end subroutine damping


subroutine pcf(C, T, B, P, E, S)
  ! Implementation of the Positive Continuous Fraction (PCF) as described
  ! by Carrayrou (2017)
  integer, intent(in) :: S, E
  double precision, dimension(S+E), intent(inout) :: C
  double precision, dimension(S), intent(in) :: T
  double precision, dimension(E), intent(in) :: B
  integer, dimension(E, S), intent(in):: P

  double precision, dimension(S) :: expn, convaux
  integer, dimension(E+S, S) :: morel
  integer :: aux, i, j
  logical :: converged
  double precision :: theta, sump, sumr

  ! find exponents
  aux = 1
  do i=1,S
    do j=1,E
      if ((P(j,i) < aux) .and. (P(j,i) > 0)) then
        aux = P(j,i)
      end if
    end do
    expn(i) = 1/aux
  end do

  ! build Morel table. REMOVE IN THE FUTURE
  do i=1,S
    do j=1,S
      if (i .eq. j )then
        morel(i,j) = 1
      else
        morel(i,j) = 0
      end if
    end do
    do j=1,E
      morel(j+S, i) = P(j, i)
    end do
  end do

  converged = .false.
  do while (converged)
    call cpluss(C, B, P, E, S)
    ! calculate sump and sumr
    do i=1,S
      if (T(i) .ge. 0.0) then
        sump = T(i)
        sumr = 0.0
      else
        sumr = -T(i)
        sump = 0.0
      end if
      do j=1,S+E
        if (P(j,i) > 0) then
          sumr = sumr + morel(j,i) * C(j)
        else
          sump = sump - morel(j,i) * C(j)
        end if
      end do
      if (sump < sumr) then
        if (sumr .ne. 0.0) then
          theta = 0.9 - 0.8 * sump / sumr
        else
          theta = 0.9
        end if
      else
        if (sump .ne. 0.0) then
          theta = 0.9 - 0.8 * sumr / sump
        else
          theta = 0.9
        end if
      end if
      ! update concentrations
      C(i) = theta * C(i) * (sump/sumr)**expn(i) + (1-theta)*C(i)
      ! check convergence
      convaux(i) = abs(sump-sumr)/(sump + sump)
    end do
    converged = maxval(convaux, 1) .lt. 0.25
  end do
end subroutine pcf


subroutine interp(y, n)
  integer, intent(in) :: n
  double precision, dimension(4, n), intent(inout) :: y

  integer :: i
  do i=1,n
    y(4,n) = y(3,n) + (y(3,n)-y(2,n))**2/(y(2,n)-y(1,n))
  end do
end subroutine interp
