      subroutine fcurvp2(t,n,x,y,p,yp,sigma,result);

      double precision function curvp2
      double precision    t
      integer n
      double precision    x(n)
      double precision    y(n)
      double precision    p
      double precision    yp(n)
      double precision    sigma
      double precision    result

      result=curvp2(t,n,x,y,p,yp,sigma);

      return
      end

      subroutine fcurv2(t,n,x,y,yp,sigma,result);

      double precision function curv2
      double precision    t
      integer n
      double precision    x(n)
      double precision    y(n)
      double precision    yp(n)
      double precision    sigma
      double precision    result

      result=curv2(t,n,x,y,yp,sigma);

      return
      end

      subroutine fdcurvp2(t,n,x,y,p,yp,sigma,result);

      double precision function dcurvp2
      double precision    t
      integer n
      double precision    x(n)
      double precision    y(n)
      double precision    p
      double precision    yp(n)
      double precision    sigma
      double precision    result

      result=dcurvp2(t,n,x,y,p,yp,sigma);

      return
      end
