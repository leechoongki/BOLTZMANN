      subroutine Eloc(egs,sxx)

      use parameters
      use spin_hidn

      implicit none
      integer k
      double complex egs, sxx, sxx_tmp
      double complex,   external:: prd
      double complex, external:: zcosh,zexp
      double complex zqrt 
      parameter(zqrt = dcmplx(0.25D0,0.D0))
      double complex coshthm(Nh)

      do k = 1, Nh
         coshthm(k) = zcosh(theta(k)) 
      enddo

      egs = dcmplx(0.D0,0.D0)
      sxx = dcmplx(0.D0,0.D0)

      do k = 1, Nv-1
         if(spin(k).and.spin(k+1)) then
            egs = egs - zqrt
         elseif(spin(k).and.(.not.spin(k+1))) then
            egs = egs + zqrt
         elseif((.not.spin(k)).and.spin(k+1)) then
            egs = egs + zqrt
         else
            egs = egs - zqrt
         endif
      enddo

      if(spin(Nv).and.spin(1)) then
         egs = egs - zqrt
      elseif(spin(Nv).and.(.not.spin(1))) then
         egs = egs + zqrt
      elseif((.not.spin(Nv)).and.spin(1)) then
         egs = egs + zqrt
      else
         egs = egs - zqrt
      endif

      do k = 1, Nv
         if((spin(k))) then
             sxx_tmp = 0.5D0*zexp(                           &
                     - 2.D0*dcmplx(prmts(2*k-1),prmts(2*k))) & 
                     * prd(k,coshthm)
             sxx = sxx + sxx_tmp
             egs = egs - hfield*sxx_tmp
         else                                     
             sxx_tmp = 0.5D0*zexp(                           &
                       2.D0*dcmplx(prmts(2*k-1),prmts(2*k))) & 
                     * prd(k,coshthm)
             sxx = sxx + sxx_tmp
             egs = egs - hfield*sxx_tmp
         endif
      enddo

      return
      end subroutine Eloc 

      double complex function prd(ii,coshthm)

      use parameters
      use spin_hidn
      implicit none
      integer ii      ! ii - spin up, kk - spin down
      integer k
      double complex coshthm(Nh)
      double complex, external:: zcosh

      prd = (1.D0,0.D0)
      do k = 1, Nh
         if(spin(ii)) then
            prd = prd*zcosh( theta(k)                          &
                - 2.0*dcmplx(prmts(kshift+2*ii-1+(k-1)*2*Nv)   &
                            ,prmts(kshift+2*ii  +(k-1)*2*Nv))) &
              /coshthm(k)
         else
            prd = prd*zcosh( theta(k)                          &
                + 2.0*dcmplx(prmts(kshift+2*ii-1+(k-1)*2*Nv)   &
                            ,prmts(kshift+2*ii  +(k-1)*2*Nv))) &
              /coshthm(k)
         endif
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

