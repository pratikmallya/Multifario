
c   %W%
c   %D% %T%

c  Please refer to the LICENSE file in the top directory*/

c       author: Mike Henderson mhender@watson.ibm.com 

      subroutine dbosvd(abd,lda, n,  ml, mu,b,ldb, nb,c,ldc,d,ldd,x,s,
     *       vr,work,lwork,ipvt,info)

      integer lda,n,ml,mu,ipvt(1),info                                          
      integer ldb,nb,ldc,ldd                                                    
      double precision abd(lda,1)                                               
      double precision b(ldb,1)                                                 
      double precision c(ldc,1)                                                 
      double precision d(ldd,1)                                                 
      double precision x(nb+ml+1,nb+ml+1)                                       
      double precision s(nb+ml+1)
      double precision vr(nb+ml+1,nb+ml+1)                                       
      double precision work(lwork)
      double precision U(1)                                               

c     m.henderson 07/18/2003 modified dbofa to call svd instead of factor

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

      call dgesvd('N','A',ml+1+nb,ml+1+nb,x,ml+1+nb,s,U,ml+1+nb,
     *                                                vr,ml+1+nb,
     *        work,lwork,info)

      return                                                                    
      end                                                                       
