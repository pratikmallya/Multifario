
c   %W%
c   %D% %T%

c   Please refer to the LICENSE file in the top directory


c       author: Mike Henderson mhender@watson.ibm.com

      subroutine dbosl(abd,lda,n,ml,mu,b,ldb,nb,c,ldc,d,ldd,x,y,
     *                  r,s,xsave,ipvt,job)
      double precision ysave(20000),e
      integer lda,n,ml,mu,ipvt(n+nb),job
      integer ldb,nb,ldc,ldd
      double precision abd(lda,1)
      double precision b(ldb,1)
      double precision c(ldc,1)
      double precision d(ldd,1)
      double precision x(nb+ml+1,1)
      double precision xsave(nb+ml+1,1)
      double precision y(1)
      double precision r(1)
      double precision s(1)

c     m.henderson 12/86

c     dbosl solves the double precision band system
c     a * x=b  or  trans(a) * x=b
c     using the factors computed by dgbco or dgbfa.

c     on entry

c        abd     double precision(lda, n)
c                the output from dgbco or dgbfa.

c        lda     integer
c                the leading dimension of the array  abd.

c        n       integer
c                the order of the original matrix.

c        ml      integer
c                number of diagonals below the main diagonal.

c        mu      integer
c                number of diagonals above the main diagonal.

c        ipvt    integer(n)
c                the pivot vector from dgbco or dgbfa.

c        b       double precision(n)
c                the right hand side vector.

c        job     integer
c               =0         to solve  a*x=b ,
c               =nonzero   to solve  trans(a)*x=b , where
c                            trans(a)  is the transpose.

c     on return

c        b       the solution vector  x.

c     error condition

c        a division by zero will occur ifthe input factor contains a
c        zero on the diagonal. technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda. it will not occur ifthe subroutines are
c        called correctly and ifdgbco has set rcond.gt.0.0
c        or dgbfa has set info.eq.0.

      double precision ddot,t
      integer k,kb,l,la,lb,lm,m,nm1

      m=mu + ml + 1
      nm1=n - ml - 1
      if(job.eq.0)then

c       job=0 , solve  a * x=r
c       first solve l*y=r

         if(ml.ne.0.and.nm1.ge.1)then
            do k=1, nm1
               lm=min(ml,n-k-1)
               l=ipvt(k)
               t=r(l)
               if(l.ne.k)then
                  r(l)=r(k)
                  r(k)=t
               endif
               call daxpy(lm,t,abd(m+1,k),1,r(k+1),1)
               call daxpy(nb,t,    c(1,k),1,  s(1),1)
            enddo
         endif

         kb=1
         do k=n-ml,n
            y(kb)=r(k)
            ysave(kb)=y(kb)
            kb=kb+1
         enddo

         do k=1,nb
            y(kb)=s(k)
            ysave(kb)=y(kb)
            kb=kb+1
         enddo

         call dgeslmod(x,nb+ml+1,nb+ml+1,ipvt(n-ml),y,job)

c        write(6,*)'dgeslmod returned y, Xy should be b!'
c        write(6,*)' eq      Ax.    b  .|Ax-b| .   x   .'
c        do i=1,nb+ml+1
c         e=0.
c         do j=1,nb+ml+1
c          e=e+xsave(i,j)*y(j)
c         enddo
c         write(6,'(i3,4e10.3)')i,e,ysave(i),abs(e-ysave(i)),y(i)
c        enddo

         do k=1,ml+1
            r(n-ml-1+k)=y(k)
         enddo

         do k=1,nb
            s(k)=y(ml+k+1)
         enddo

c        now solve  u*x=y

         do k=1,nb
            t=-s(k)
            call daxpy(n-ml-1,t,b(1,k),1,r(1),1)
         enddo

         do k=n,n-ml,-1
            la=max(1,mu+ml+1-k)
            lm=min(mu+ml+1,n-k+mu) - la + 1
            lb=n - ml - lm
            t=-r(k)
            call daxpy(lm,t,abd(la,k),1,r(lb),1)
         enddo

         do k=n-ml-1,1,-1
            r(k)=r(k)/abd(mu+ml+1,k)
            la=max(1,mu+ml+2-k)
            lm=min(mu+ml+1,n-k+mu) - la
            lb=k-lm
            t=-r(k)
            call daxpy(lm,t,abd(la,k),1,r(lb),1)
         enddo

      else

         do k=1, nm1
            lm=min(k,m)-1
            la=m - lm
            lb=k - lm
            t=ddot(lm,abd(la,k),1,r(lb),1)
            r(k)=(r(k) - t)/abd(m,k)
         enddo

         do k=n,n-ml,-1
            lm=m-k+n-ml-1
            lb=k-m+1
            t= ddot(lm,abd(1,k),1,r(lb),1)
            r(k)=r(k)-t
         enddo

         do k=1,nb
            t=ddot(n-ml-1,b(1,k),1,r(1),1)
            s(k)=s(k)-t
         enddo

         kb=1
         do k=n-ml,n
            y(kb)=r(k)
            kb=kb+1
         enddo

         do k=1,nb
            y(kb)=s(k)
            kb=kb+1
         enddo

         call dgeslmod(x,nb+ml+1,nb+ml+1,ipvt(n-ml),y,job)

         do k=1,ml+1
            r(n-ml-1+k)=y(k)
         enddo

         do k=1,nb
            s(k)=y(ml+k+1)
    3    enddo

         if(ml.ne.0.and.nm1.ge.1)then
           do kb=1, nm1
             k=nm1 + 1 - kb
             lm=min0(ml,n-k)
             r(k)=r(k)+ddot(lm,abd(m+1,k),1,r(k+1),1)
     *                +ddot(nb,c(1,k),1,s(1),1)
               l=ipvt(k)
               if(l.ne.k)then
                  t=r(l)
                  r(l)=r(k)
                  r(k)=t
               endif
            enddo
         endif

      endif
      return
      end

      subroutine dgeslmod(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)

      double precision ddot,t
      integer k,kb,l,nm1

      nm1=n - 1
      if(job.eq.0)then
        call dlaswp(1,b,n,1,n,ipvt,1)
        call dtrsm('Left','Lower','No transpose','Unit',
     *                                  n,1,1.d0,a,lda,b,n)
        call dtrsm('Left','Upper','No transpose','Non-unit',
     *                                  n,1,1.d0,a,lda,b,n)
      else
        call dtrsm('Left','Upper','Transpose','Non-unit',
     *                                  n,1,1.d0,a,lda,b,n)
        call dtrsm('Left','Lower','Transpose','Unit',
     *                                  n,1,1.d0,a,lda,b,n)
        call dlaswp(1,b,n,1,n,ipvt,-1)
      endif
      return
      end
