
c   %W%
c   %D% %T%

c  Please refer to the LICENSE file in the top directory*/

c       author: Mike Henderson mhender@watson.ibm.com 

      subroutine dbofa(abd,lda,n,ml,mu,b,ldb,nb,c,ldc,d,ldd,x,xsave,
     *                                                    ipvt,info)
      integer lda,n,ml,mu,ipvt(1),info                                          
      integer ldb,ldc,ldd                                                    
      double precision abd(lda,1)                                               
      double precision b(ldb,1)                                                 
      double precision c(ldc,1)                                                 
      double precision d(ldd,1)                                                 
      double precision x(nb+ml+1,nb+ml+1)                                       
      double precision xsave(nb+ml+1,nb+ml+1)                                       

c     m.henderson 12/86                   

c     dbofa factors a double precision band matrix with border                  
c     by elimination.                                                          

c     on entry                                                                  

c        abd     double precision(lda, n)                                       
c                contains the matrix in band storage. the columns              
c                of the matrix are stored in the columns of  abd  and           
c                the diagonals of the matrix are stored in rows                 
c                ml+1 through 2*ml+mu+1 of  abd.                              
c                see the comments below for details.                           

c        lda     integer                                                        
c                the leading dimension of the array  abd.                     
c                lda must be.ge.2*ml + mu + 1.                              

c        n       integer                                                        
c                the order of the original matrix.                             

c        ml      integer                                                        
c                number of diagonals below the main diagonal.                  
c                0.le.ml.lt.n.                                            

c        mu      integer                                                        
c                number of diagonals above the main diagonal.                  
c                0.le.mu.lt.n.                                            
c                more efficient if ml.le.mu.                               



c        b       double precision(ldb, nb)                                      
c                contains the right border                                      

c        ldb     integer                                                        
c                the leading dimension of the array  b.                        
c                ldb must be.ge.n                                             

c        nb      integer                                                        
c                the size of the border                                         


c        c       double precision(ldc, n)                                       
c                contains the lower left border                                 

c        ldc     integer                                                        
c                the leading dimension of the array  c.                        
c                ldc must be.ge.nb                                            


c        d       double precision(ldd, nb)                                      
c                contains the lower right border                                

c        ldd     integer                                                        
c                the leading dimension of the array  d.                        
c                ldd must be.ge.nb                                            


c        x       double precision(nb+mu+1,nb+mu+1)                              
c                work array                                                     

c     on return                                                                 

c        abd     an upper triangular matrix in band storage and                 
c         +      the multipliers which were used to obtain it.                 
c         x      the factorization can be written  a=l*u  where               
c                l  is a product of permutation and unit lower                  
c                triangular matrices and  u  is upper triangular.              

c        ipvt    integer(n+m)                                                   
c                an integer vector of pivot indices.                           

c        info    integer                                                        
c               =0  normal value.                                            
c               =k  if u(k,k).eq.0.0. this is not an error               
c                     condition for this subroutine, but it does                
c                     indicate that dgbsl will divide by zero if               
c                     called. use  rcond  in dgbco for a reliable              
c                     indication of singularity.                               


      double precision t                                                        
      integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1                     


      m=ml + mu + 1                                                           
      info=0                                                                  

      j0=mu + 2                                                               
      j1=min0(n,m) - 1                                                        
      if(j1.ge.j0)then
        do jz=j0, j1                                                         
          i0=m + 1 - jz                                                        
          do i=i0, ml                                                       
            abd(i,jz)=0.0d0                                                   
          enddo                                                               
        enddo
      endif
      jz=j1                                                                   
      ju=0                                                                    

      nm1=n - ml - 1                                                          
      if(nm1.ge.1)then
        do k=1, nm1                                                         
          kp1=k + 1                                                            

          jz=jz + 1                                                            
          if(jz.le.n.and.ml.ge.1)then
            do i=1, ml                                                     
              abd(i,jz)=0.0d0                                                
            enddo
          endif

          lm=min0(ml,n-k)                                                      
          l=idamax(lm+1,abd(m,k),1) + m - 1                                    
          ipvt(k)=l + k - m                                                    

          if(abd(l,k).ne.0.0d0)then

            if(l.ne.m)then
              t=abd(l,k)                                                     
              abd(l,k)=abd(m,k)                                              
              abd(m,k)=t                                                     
            endif

            t=-1.0d0/abd(m,k)                                                 
            call dscal(lm,t,abd(m+1,k),1)                                       
            call dscal(nb,t,c(1,k),1)                                           

            do j=1,nb                                                      
              t=b(ipvt(k),j)                                                 
              if(k.ne.ipvt(k))then
                b(ipvt(k),j)=b(k,j)                                         
                b(k,j)=t                                                    
              endif
              lm=min(ml,n-k-1)                                                 
              call daxpy(lm,t,abd(m+1,k),1,b(k+1,j),1)                         
              call daxpy(nb,t,    c(1,k),1,  d(1,j),1)                         
            enddo

            ju=min0(max0(ju,mu+ipvt(k)),n)                                    
            if(ju.ge.kp1)then
              mm=m                                                              
              do j=kp1, ju                                                   
                l=l - 1                                                        
                mm=mm - 1                                                      
                t=abd(l,j)                                                     
                if(l.ne.mm)then
                  abd(l,j)=abd(mm,j)                                          
                  abd(mm,j)=t                                                 
                endif
                call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)                      
                call daxpy(nb,t,    c(1,k),1,     c(1,j),1)                      
              enddo
            endif
         else
          info=k                                                            
         endif
        enddo
      endif

      do i=1,ml+1                                                           
        do j=1,ml+1                                                        
           x(i,j)=abd(ml+mu+i-j+1,n-ml-1+j)                                    
        enddo
        do j=1,nb                                                          
           x(i,ml+1+j)=b(n-ml-1+i,j)                                           
           x(ml+1+j,i)=c(j,n-ml-1+i)                                           
        enddo
      enddo

      do i=1,nb                                                             
        do j=1,nb                                                          
          x(ml+1+i,ml+1+j)=d(i,j)                                             
        enddo
      enddo

c     write(6,*)'dbofa, system is ',n+ldc,'x',n+ldc,
c    *              ' LR is ',nb+ml+1,'x',nb+ml+1
c     write(6,*)'unfactored lower right corner:'
c     do i=1,nb+ml+1
c      write(6,'(20f10.7)')(x(i,j),j=1,nb+ml+1)
c      write(6,'(1p20e7.0)')(x(i,j),j=1,nb+ml+1)
c     enddo

      if(loc(xsave).ne.0)then
        do i=1,nb+ml+1
         do j=1,nb+ml+1
           xsave(i,j)=x(i,j)
         enddo
        enddo
      endif

      call dgetrf(ml+1+nb,ml+1+nb,x,ml+1+nb,ipvt(n-ml),info)                             

      return                                                                    
      end                                                                       
