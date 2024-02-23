C FILE: nr.f
C Subroutines borrowed from Numerical Recipes
C
      SUBROUTINE ludcmp(a,n,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(n,n),TINY
Cf2py intent(inout) a
Cf2py intent(in) n, indx, NMAX, d
C     Largest expected n,  and  a small  number.
      PARAMETER (NMAX=500,TINY=1.0e-20)
C     Given a matrix a(1:n,1:n), with physical dimension np by np, this routine replaces it by
C     the LU decomposition of a rowwise permutation of itself. a and n
C     are input. a is  output, arranged  as  in  equation (2.3.14) above; indx(1:n)
C     is an output vector that records the row permutation effected by the partial pivoting;
C     d is output as +/-1 depending on whether the  number  of  row  interchanges  was  
C     even or odd, respectively. This routine is used in combination  with
C     lubksb to  solve linear  equations  or invert  a  matrix.
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)
C     vv stores the  implicit  scaling  of each row.  
      d=1.
C     No  row interchanges  yet.
      do 12 i=1,n
C       Loop  over  rows  to  get  the  implicit  scaling  information.
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      enddo
C       No nonzero largest element.
C       if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
C       Save  the  scaling.
12    enddo
      do 19 j=1,n
C       This is the loop over columns of Crout's method.
        do 14 i=1,j-1
C         This  is equation (2.3.12) except for i = j.
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        enddo
          a(i,j)=sum
14      enddo
      aamax=0.
C     Initialize  for the search for largest pivot element.
      do 16 i=j,n
C       This is i = j of equation (2.3.12) and i = j +1 :::N of  equation (2.3.13).
        sum=a(i,j)
        do 15 k=1,j-1
          sum=sum-a(i,k)*a(k,j)
15      enddo
        a(i,j)=sum
        dum=vv(i)*abs(sum)
C       Figure  of  merit  for  the pivot.
        if (dum.ge.aamax) then
C         Is it  better than  the  best so far?
          imax=i
          aamax=dum
        endif
16    enddo
      if (j.ne.imax)then
C       Do  we  need to  interchange rows?
        do 17 k=1,n
C         Yes,  do  so...
          dum=a(imax,k)
          a(imax,k)=a(j,k)
          a(j,k)=dum
17      enddo
        d=-d
C       ...and change the  parity  of d.
        vv(imax)=vv(j)
C       Also interchange the  scale factor.
      endif
      indx(j)=imax
      if(a(j,j).eq.0.)a(j,j)=TINY
C       If the pivot element is zero the matrix is singular (at least to the
C       precision of the algorithm). For some applications on singular matrices,
C       it is desirable to substitute TINY for  zero.
        if (j.ne.n) then
C         Now, finally, divide by the pivot element.
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        end do
        end if
19    end do
C     Go back for the next column in  the reduction.
      return
      END

C     Here is the routine for forward substitution and backsubstitution, implementing
C     equations  (2.3.6)  and  (2.3.7).
      SUBROUTINE lubksb(a,n,indx,b)
      INTEGER n,indx(n)
      REAL a(n,n),b(n)
Cf2py intent(in) a, n, np, indx
Cf2py intent(out) b
C     Solves the  set of n linear  equations AX = B. Here a is input,  not  as the  matrix
C     A but rather as its LU decomposition, determined by the routine ludcmp .
C     indx is input as the permutation vector returned by ludcmp.  b(1:n)
C     is input as the right-hand side vector B, and returns with the solution vector
C     X. a, n, np, and indx are not modi ed by this routine
C     and can be left in place for successive calls with  different right-hand sides
C     b.  This  routine takes into account the possibility that b
C     will begin with many zero elements, so it is efficient for use in matrix inversion.
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
C     When ii is set to a positive value, it will become the in- dex of
C     the rest nonvanishing element of b. We now do the forward substitution,
C     equation (2.3.6). The only new wrinkle is to unscramble the permutation as we go.
      do 22 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 21 j=ii,i-1
          sum=sum-a(i,j)*b(j)
21      enddo
        else if (sum.ne.0.) then
          ii=i
C         A nonzero element was encountered, so from now on we will
C         have  to do  the sums  in  the  loop above.
        endif
        b(i)=sum
22    enddo
      do 24 i=n,1,-1
C       Now we do the backsubstitution, equation  (2.3.7).
        sum=b(i)
        do 23 j=i+1,n
          sum=sum-a(i,j)*b(j)
23      enddo
        b(i)=sum/a(i,i)
C       Store a  component of the solution vector X.
24    enddo
      
C     To  summarize,  this is the preferred  way  to solve the linear  set  of  equations
C     A x = b :
C     call ludcmp(a,n,np,indx,d)
C     call lubksb(a,n,np,indx,b)
C     The  answer x will be returned in b.  Your original matrix A will have been destroyed.
C     If you subsequently want to solve a set  of equations with the same
C     A but a different  right-hand side b ,  you  repeat only
C     call lubksb(a,n,np,indx,b)
C     not,  of course,  with the original matrix A ,  but with a and indx
C     as  were  already returned  from ludcmp
      
      END

      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
      INTEGER n
      LOGICAL check
      REAL f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX
Cf2py intent(in) n, xolf, fold, g, p, x, f, stpmax, check, func
Cf2py intent(out) g, x
      PARAMETER (ALF=1.e-4,TOLX=1.e-7)
      EXTERNAL func
C  USES func
C     Given an n-dimensional point xold(1:n), the value of the function and gradient there,
C     fold and g(1:n), and a direction p(1:n), nds a new point x(1:n) along the direction
C     p from xold where the function func has decreased sufficiently."The new function value
C     is returned in f. stpmax is an input quantity that limits the length of the steps so that you
C     do not try to evaluate the function in regions where it is  unde ned  or subject to overflow.
C     p is usually the Newton direction. The output quantity check
C     is false on a normal exit. It is true  when x is  too  close  to xold.
C     In a minimization  algorithm, this usually signals
C     convergence and can be  ignored.  However, in  a zero-finding algorithm the calling  program
C     should check  whether  the  convergence  is  spurious.
C     Parameters:
C      ALF       ensures  sufficient  decrease  in  function  value;
C      TOLX is  the  convergence criterion  on Delta x .
      INTEGER i
      REAL a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,
     *     sum,temp,test,tmplam
      check=.false.
      sum=0.
      do 31 i=1,n
        sum=sum+p(i)*p(i)
31    enddo

      sum=sqrt(sum)
      if(sum.gt.stpmax)then
C       Scale if  attempted step is too big.
        do 32 i=1,n
          p(i)=p(i)*stpmax/sum
32      enddo
      endif
      slope=0.
      do 33 i=1,n
        slope=slope+g(i)*p(i)
33    enddo

C     if(slope.ge.0.) pause 'roundoff problem in lnsrch'
      test=0.
C     Compute lambda min.
      do 34 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.)
        if(temp.gt.test)test=temp
34    enddo

      alamin=TOLX/test
      alam=1.
C     Always try full Newton step  rst.
1     continue
C       Start of iteration loop.
        do 35 i=1,n
          x(i)=xold(i)+alam*p(i)
35      enddo

        f=func(x)
        if(alam.lt.alamin)then
C         Convergence  on Delta x. For  zero finding, the calling program
C         should verify the convergence.
          do 36 i=1,n
            x(i)=xold(i)
36        enddo

          check=.true.
          return
      else if(f.le.fold+ALF*alam*slope)then
C        Sufficient function decrease.
         return
      else
C       Backtrack.
        if(alam.eq.1.)then
C         First  time.
          tmplam=-slope/(2.*(f-fold-slope))
        else
C         Subsequent backtracks.
          rhs1=f-fold-alam*slope
          rhs2=f2-fold-alam2*slope
          a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
          b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/
     *      (alam-alam2)
          if(a.eq.0.)then
            tmplam=-slope/(2.*b)
          else
            disc=b*b-3.*a*slope
            if(disc.lt.0.)then
              tmplam=.5*alam
            else if(b.le.0.)then
              tmplam=(-b+sqrt(disc))/(3.*a)
            else
              tmplam=-slope/(b+sqrt(disc))
            endif
          endif
          if(tmplam.gt..5*alam)tmplam=.5*alam  ! lambda <= 0.5 lambda_1
        endif
      endif
      alam2=alam
      f2=f
      alam=max(tmplam,.1*alam) ! lambda >= 0.1 lambda_1
      goto 1  !   Try  again.
      END
