      double precision function Eloc()

      use parameters
      use spin_hidn

      implicit none
      integer k
      double precision, external:: prd, logistic2

      Eloc = 0.D0
      do k = 1, Nv-1
         Eloc = Eloc - 0.25D0*spin(k)*spin(k+1) 
      enddo

      Eloc = Eloc - 0.25D0*spin(Nv)*spin(1)

      do k = 1, Nv
         Eloc = Eloc - 0.5D0*hfield*(1.D0/logistic2(prmts(k)*spin(k))-1.D0) &
                *prd(k,spin(k)) 
      enddo 
      return
      end function Eloc 

      double precision function prd(ii,sigma)

      use parameters
      use spin_hidn, ONLY: theta
      implicit none
      integer ii
      integer k   ! index for counting
      integer sigma
      double precision xx, yy
      double precision,external:: logistic2

!     lgs(x) = 1 / ( 1 + exp(-x)) ==> exp(-x) = 1/lgs(x) - 1
!     
!     cosh (x) / cosh (y) = (exp(x) + exp(-x))/(exp(y) + exp(-y))
!                         = exp(x-y)*( 1 + exp (-2x) )/(1 + exp(-2y))
!                         = (1/log(y-x)-1)/log(2x)*log(2y)

      prd = cosh(theta(1) - 2.D0*sigma*prmts(kshift+ii))/cosh(theta(1))

!     xx = theta(1) - 2.D0*sigma*wght(1,ii)
!     xx = theta(1) - 2.D0*sigma*prmts(kshift+ii)
!     yy = theta(1)
!     prd = (1.D0/logistic2(yy-xx)-1.D0)/logistic2(2.D0*xx)*&
!           logistic2(2.D0*yy)
      do k = 2, Nh
         xx = theta(k) - 2.D0*sigma*prmts(kshift+ii+(k-1)*Nv)
         yy = theta(k)
         prd = prd*cosh(theta(k) - 2.D0*sigma*prmts(kshift+ii+(k-1)*Nv))/cosh(theta(k))
!        write(6,*) k,theta(k), prmts(kshift+ii+(k-1)*Nv),cosh(theta(k))
       
!        prd = prd*(1.D0/logistic2(yy-xx)-1.D0)/logistic2(2.D0*xx)*&
!                   logistic2(2.D0*yy)
!        write(6,*) k,logistic2(yy-xx),logistic2(2.D0*xx),logistic2(2.D0*yy)
      enddo

      prd = dsqrt(prd)

      return
      end function prd
