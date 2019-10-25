      double precision function Eloc(spin,Nv,Nh,wght,theta,hfield,Jex,a,b)

      implicit none
      integer k
      integer Nv, Nh
      integer spin(Nv)
      double precision hfield
      double precision Jex
      double precision theta(Nh)
      double precision a(Nv),b(Nh)
      double precision wght(Nh,Nv)
      double precision, external:: prd, logistic2

      Eloc = 0.D0
      do k = 1, Nv-1
         Eloc = Eloc - Jex*spin(k)*spin(k+1) 
      enddo

      Eloc = Eloc - Jex*spin(Nv)*spin(1)

      do k = 1, Nv
!        Eloc = Eloc - hfield*exp( -a(k)*spin(k)) &
!               *prd(theta,Nh,Nv,k,wght,spin(k)) 
         Eloc = Eloc - hfield*(1.D0/logistic2(a(k)*spin(k))-1.D0) &
                *prd(theta,Nh,Nv,k,wght,spin(k)) 
      enddo 
      return
      end function Eloc 

      double precision function prd(theta,Nh,Nv,ii,wght,sigma)

      implicit none
      integer ii
      integer k   ! index for counting
      integer sigma
      integer Nh, Nv
      double precision wght(Nh,Nv)
      double precision theta(Nh)
      double precision xx, yy
      double precision,external:: logistic2

!     lgs(x) = 1 / ( 1 + exp(-x)) ==> exp(-x) = 1/lgs(x) - 1
!     
!     cosh (x) / cosh (y) = (exp(x) + exp(-x))/(exp(y) + exp(-y))
!                         = exp(x-y)*( 1 + exp (-2x) )/(1 + exp(-y))
!                         = (1/log(y-x)-1)/log(2x)*log(2y)

!     prd = cosh(theta(1) - 2.D0*sigma*wght(1,ii))/cosh(theta(1))
      xx = theta(1) - 2.D0*sigma*wght(1,ii)
      yy = theta(1)
      prd = (1.D0/logistic2(yy-xx)-1.D0)/logistic2(2.D0*xx)*&
            logistic2(2.D0*yy)
      do k = 2, Nh
         xx = theta(k) - 2.D0*sigma*wght(k,ii)
         yy = theta(k)
!        prd = prd*cosh(theta(k) - 2.D0*sigma*wght(k,ii))/cosh(theta(k))
         prd = prd*(1.D0/logistic2(yy-xx)-1.D0)/logistic2(2.D0*xx)*&
                    logistic2(2.D0*yy)
      enddo

      prd = dsqrt(prd)

      return
      end function prd
