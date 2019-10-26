      double complex function Eloc()

      use parameters
      use spin_hidn

      implicit none
      integer k
      double complex,   external:: prd
      double complex, external:: zcosh,zexp
      double complex zqrt 
      parameter(zqrt = dcmplx(0.25D0,0.D0))
      double complex coshthm(Nh)

      do k = 1, Nh
         coshthm(k) = zcosh(theta(k)) 
      enddo

      Eloc = 0.D0
      do k = 1, Nv-1
         if(spin(k).and.spin(k+1)) then
            Eloc = Eloc + zqrt
         elseif(spin(k).and.(.not.spin(k+1))) then
            Eloc = Eloc - zqrt
         elseif((.not.spin(k)).and.spin(k+1)) then
            Eloc = Eloc - zqrt
         else
            Eloc = Eloc + zqrt
         endif
      enddo

      if(spin(Nv).and.spin(1)) then
         Eloc = Eloc + zqrt
      elseif(spin(Nv).and.(.not.spin(1))) then
         Eloc = Eloc - zqrt
      elseif((.not.spin(Nv)).and.spin(1)) then
         Eloc = Eloc - zqrt
      else
         Eloc = Eloc + zqrt
      endif

      do k = 1, Nv-1
         if((spin(k)).and.(.not.spin(k+1))) then
             Eloc = Eloc + 0.5D0*zexp(-2.0*prmts(k) &
                   + 2.0*prmts(k+1))*prd(k,k+1,coshthm)  
         elseif((.not.spin(k)).and.(spin(k+1))) then
             Eloc = Eloc + 0.5D0*zexp( 2.0*prmts(k) &
                   - 2.0*prmts(k+1))*prd(k+1,k,coshthm)
         endif
      enddo

      if((spin(Nv)).and.(.not.spin(1))) then
          Eloc = Eloc + 0.5D0*zexp(-2.0*prmts(Nv)  &
               + 2.0*prmts(1))*prd(Nv,1,coshthm) 
      elseif((.not.spin(Nv)).and.(spin(1))) then
          Eloc = Eloc + 0.5D0*zexp( 2.0*prmts(Nv) &
               - 2.0*prmts(1))*prd(1,Nv,coshthm) 
      endif
      
      return
      end function Eloc 

      double complex function prd(ii,kk,coshthm)

      use parameters
      use spin_hidn, ONLY: theta
      implicit none
      integer ii,kk   ! ii - spin up, kk - spin down
      integer k
      double complex coshthm(Nh)
      double complex, external:: zcosh

!     lgs(x) = 1 / ( 1 + exp(-x)) ==> exp(-x) = 1/lgs(x) - 1
!     
!     cosh (x) / cosh (y) = (exp(x) + exp(-x))/(exp(y) + exp(-y))
!                         = exp(x-y)*( 1 + exp (-2x) )/(1 + exp(-2y))
!                         = (1/log(y-x)-1)/log(2x)*log(2y)

      prd = 1.D0

      do k = 1, Nh
         prd = prd*zcosh( theta(k)-2.0*prmts(kshift+ii+(k-1)*Nv)    &
              + 2.0*prmts(kshift+kk+(k-1)*Nv) ) / coshthm(k)
      enddo

      return
      end function prd

      double complex function zcosh(z)

      implicit none
      double complex z
      double precision x, y
      double precision, external :: cosl, sinl, coshl, sinhl

      ! 2 cosh(x+iy) = exp(x)*(cos(y)+isin(y))+exp(-x)*(cos(y)-isin(y))
      !              = 2cosh(x)*cos(y) + 2isinh(x)*sin(y)

      x = real(z)
      y = aimag(z)

      zcosh =dcmplx(coshl(x)*cosl(y), sinhl(x)*sinl(y))

      return
      end function zcosh

      double complex function ztanh(z)
 
      implicit none
      double complex z
      double precision x, y
      double precision, external:: tanhl, tanl

      x = real(z)
      y = aimag(z)

      ztanh = dcmplx(tanhl(x),tanl(y))/dcmplx(1.D0,tanhl(x)*tanl(y))
      return
      end function ztanh
   

      double complex function zexp(z)

      implicit none
      double complex z 
      double precision, external:: sinl, cosl,expl
      
      zexp = expl(real(z))*dcmplx( cosl(aimag(z)), sinl(aimag(z)))
      return
      end function zexp

      double precision function tanhl(x)
      implicit none
      double precision x
      double precision, external:: sinhl, coshl
  
      tanhl = sinhl(x)/coshl(x)
      return
      end function tanhl

      double precision function tanl(x)
      implicit none
      double precision x
      double precision, external:: sinl, cosl

      tanl = sinl(x)/cosl(x)

      return
      end function tanl

