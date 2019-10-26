      double complex function Eloc()

      use parameters
      use spin_hidn

      implicit none
      integer k,i
      integer inn
      double complex,   external:: prd
      double complex, external:: zcosh,zexp
      double complex zqrt 
      parameter(zqrt = dcmplx(0.25D0,0.D0))
      double complex coshthm(Nh)

      do k = 1, Nh
         coshthm(k) = zcosh(theta(k)) 
      enddo

      Eloc = dcmplx(0.D0,0.D0)
      do k = 1, Nv
         do i = 1, nz
            inn = innd(k,i)
            if(spin(k).and.spin(inn)) then
               Eloc = Eloc + zqrt
            elseif(spin(k).and.(.not.spin(inn))) then
               Eloc = Eloc - zqrt
            elseif((.not.spin(k)).and.spin(inn)) then
               Eloc = Eloc - zqrt
            else
               Eloc = Eloc + zqrt
            endif
         enddo
      enddo

      Eloc = 0.5D0*Eloc

      do k = 1, Nv
         do i = 1, Nz
            inn = innd(k,i)
            if((spin(k)).and.(.not.spin(inn))) then
               Eloc = Eloc + 0.25D0*zexp(-2.0*prmts(k) &
                    + 2.0*prmts(inn))*prd(k,inn,coshthm)  
            elseif((.not.spin(k)).and.(spin(inn))) then
               Eloc = Eloc + 0.25D0*zexp( 2.0*prmts(k) &
                   - 2.0*prmts(inn))*prd(inn,k,coshthm)
            endif
         enddo
      enddo

      
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

      prd = dcmplx(1.D0,0.D0)

      do k = 1, Nh
         prd = prd*zcosh( theta(k)-2.D0*prmts(kshift+ii+(k-1)*Nv)    &
              + 2.D0*prmts(kshift+kk+(k-1)*Nv) ) / coshthm(k)
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

