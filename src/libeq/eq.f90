! EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
!
! Compile with $ f2py3 -llapack -c -m mfeq eq.f90
!
! control codes
! 0 -> normal exit, no errors
! 1 -> 

subroutine newtonr1(C, T, B, P, CTRL, N, E, S, ES)
  ! Perform Newton-Raphson iterations per point when the beta matrix
  ! is 1D. The first initial value must be provided in the first row
  ! of C. The following points take as starting value the value of 
  ! the last iteration.

  implicit none
  integer, intent(in) :: N, E, S, ES
  ! N -> number of experimental points
  ! E -> number of equilibria
  ! S -> number of species
  integer, intent(inout) :: CTRL
  ! control value. 0 on no error, 1 otherwise
  integer, dimension(E, S), intent(in) :: P
  ! The stoichiometry array
  double precision, dimension(N, ES), intent(inout) :: C
  ! C -> (inout) initial estimates and refined values
  double precision, dimension(N, S), intent(in) :: T
  ! T -> The stoichiometry values
  double precision, dimension(E), intent(in) :: B
  integer i, ctrl2

  ctrl2 = ctrl

  do i=1,N
    if (i > 1) then
      C(i,:) = C(i-1,:)
    end if
    if (i == 1) then
      call pcf(C(i,:), T(i,:), B, P, 1E-9, E, S)
    else
      call pcf(C(i,:), T(i,:), B, P, 0.25, E, S)
    end if
    call nr(C(i,:), T(i,:), B, P, CTRL, E, S)
  end do
end subroutine newtonr1


subroutine newtonr2(C, T, B, P, CTRL, N, E, S, ES)
  ! Perform Newton-Raphson iterations per point when the beta matrix
  ! is 2D. The first initial value must be provided in the first row
  ! of C. The following points take as starting value the value of 
  ! the last iteration.

  implicit none
  integer :: N, E, S, ES
  !fpy intent(hide) :: N=shape(T,0), E=shape(P,0), S=shape(P,1)
  ! N -> (in) number of experimental points
  ! E -> (in) number of equilibria
  ! S -> (in) number of species
  integer, intent(inout) :: CTRL
  ! CTRL -> (out) control value
  integer, dimension(E, S), intent(in) :: P
  double precision, dimension(N, S), intent(in) :: T
  ! T -> The stoichiometry values
  double precision, dimension(N, E), intent(in) :: B
  ! C -> (inout) initial estimates
  double precision, dimension(N, ES), intent(inout) :: C
  integer i, ctrl2
  ctrl2 = CTRL

  !DEBUG! print *, "-----"
  !DEBUG! print *, N, E, S
  !DEBUG! print *, shape(C)
  do i=1,N
    if (i > 1) then
      C(i,:) = C(i-1,:)
    end if
    if (i == 1) then
      call pcf(C(i,:), T(i,:), B(i,:), P, 1E-9, E, S)
    else
      call pcf(C(i,:), T(i,:), B(i,:), P, 0.25, E, S)
    end if
    call nr(C(i,:), T(i,:), B(i,:), P, CTRL, E, S)
  end do
end subroutine newtonr2


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


subroutine cplussn(C, B, P, N, E, S)
  ! Compute the free concentrations for all the components.
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
  integer, intent(in) :: E, S
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
        aux = aux + P(k,i)*P(k,j)*C(k+s)
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
  ! ctrl -> control value. 
  !         0 = normal exit
  !         2 = normal exit upon initial guess being a root
  !         99 = abnormal exit on MAXITS
  integer, dimension(E, S), intent(in) :: P
  double precision, dimension(S+E), intent(inout) :: X
  double precision, dimension(S), intent(in) :: T
  double precision, dimension(E), intent(in) ::  B

  integer its, i, j
  double precision, dimension(S) :: fvec, jdiag
  double precision, dimension(S, S) ::  jac
  double precision :: f
  double precision, dimension(S) :: xold
  double precision, dimension(S) :: dx
  double precision test
  logical check
  integer MAXITS
  double precision relax
  parameter (MAXITS=100)
  
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
    dx(i) = 0.0
  end do
  if (test .lt. 1D-25) then
    check=.false.
    ctrl=2
    return
  end if

  ! Main loop
  do its=1,MAXITS
    call jacobian(jac, X, P, E, S)

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

    ! De-scale
    do i=1,S
      dx(i) = - dx(i) / sqrt(jdiag(i)) 
    end do

    do i=1,S
      relax = 1.0 / max(1.0, -2*dx(i)/x(i))
      x(i) = x(i) + relax*dx(i)
    end do
    call cpluss(X, B, P, E, S)
    call fobj(fvec, X, P, T, E, S)

    ! test for convergence on relative values
    test = 0.
    do i=1,S
       if (T(i) .ne. 0.0) then
         test = test + (fvec(i) / T(i))**2
       end if
    end do
    if ((test <= 1e-16) .and. (maxval(abs(fvec)) <= 1e-9)) then
       ctrl = 0
       !DEBUG! print *, "Normal exit", its, sum(fvec**2), x(1:2)
       return
    end if
  end do
  ctrl=99
  !! print *, "Exit on MAXITS"
end subroutine nr


subroutine pcf(C, T, B, P, THRESH, E, S)
  ! Implementation of the Positive Continuous Fraction (PCF) as described
  ! by Carrayrou (2017)
  integer, intent(in) :: S, E
  double precision, dimension(S+E), intent(inout) :: C
  double precision, dimension(S), intent(in) :: T
  double precision, dimension(E), intent(in) :: B
  real, intent(in) :: THRESH
  integer, dimension(E, S), intent(in):: P

  double precision, dimension(S) :: expn, convaux  !, DEBUG_F
  integer, dimension(E+S, S) :: morel
  integer :: aux, i, j, guard
  logical :: converged
  double precision :: theta, sump, sumr, backtrack

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
  guard = 0
  backtrack = 1.0
  do while (.not. converged)
    guard = guard + 1
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
        if (morel(j,i) > 0) then
          sumr = sumr + morel(j,i) * C(j)
        else
          sump = sump - morel(j,i) * C(j)
        end if
      end do
      if (sump < sumr) then
        theta = backtrack * (0.9 - 0.8 * sump / sumr)
      else
        theta = backtrack * (0.9 - 0.8 * sumr / sump)
      end if
      ! update concentrations
      C(i) = theta * C(i) * (sump/sumr)**expn(i) + (1-theta)*C(i)
      ! check convergence
      convaux(i) = abs(sump-sumr)/(sump + sump)
    end do

    converged = maxval(convaux) .lt. THRESH
    if (guard == 20 .or. guard == 40 .or. guard == 60) then
        backtrack = backtrack * 0.5
    end if
    if (guard >100) then
        converged = .true.
    end if
  end do
end subroutine pcf
