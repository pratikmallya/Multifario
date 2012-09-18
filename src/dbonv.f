c
c   @(#)dbonv.f	1.1
c   00/08/25 10:00:21
c
c       (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
c       CORPORATION 9/1/99.  ALL RIGHTS RESERVED.
c
c       author: Mike Henderson mhender@watson.ibm.com

      subroutine dbonv(idim,abd,lda,n,ml,mu,b,ldb,nb,c,ldc,d,ldd,
     *                                           x,s,vr,phi)
      integer idim
      integer lda,n,ml,mu,ipvt(n+nb)
      integer ldb,ldc,ldd
      double precision abd(lda,1)
      double precision b(ldb,1)
      double precision c(ldc,1)
      double precision d(ldd,1)
      double precision x(nb+ml+1,1)
      double precision vr(nb+ml+1,1)
      double precision s(1)
      double precision phi(1)
c
c     modified by m. henderson 07/17/2003  Changed to use svd
c
c     dbosl modified by m. henderson 11/00
c
c     dgbsl modified by m. henderson 12/86
c
c     dbonl finds right null vectors of the double precision band system
c
      double precision ddot,t
      integer l,la,lb,lm,m
c
      j = 0
      m=0
      do i=1,nb+ml+1
        if(abs(s(i)).lt.1.e-7)then
         j=j+1
         if(j.eq.idim+1)m=i
        endif
      enddo

      do i=1,n
        phi(i)=0.
      enddo

      do i=1,nb+ml+1
         phi(n-ml-1+i)=vr(m,i)
      enddo

c        now solve  u*x=y

      do i=1,nb
         t=-phi(n+i)
         call daxpy(n-ml-1,t,b(1,i),1,phi(1),1)
      enddo

      do i=n,n-ml,-1
         la=max(1,mu+ml+1-i)
         lm=min(mu+ml+1,n-i+mu) - la + 1
         lb=n - ml - lm
         t=-phi(i)
         call daxpy(lm,t,abd(la,i),1,phi(lb),1)
      enddo

      do i=n-ml-1,1,-1
         phi(i)=phi(i)/abd(mu+ml+1,i)
         la=max(1,mu+ml+2-i)
         lm=min(mu+ml+1,n-i+mu) - la
         lb=i-lm
         t=-phi(i)
         call daxpy(lm,t,abd(la,i),1,phi(lb),1)
      enddo

      return
      end
