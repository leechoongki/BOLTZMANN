!      double precision function Eloc()
       double complex function Eloc()

      use parameters
      use spin_hidn

      implicit none
      integer k
      double complex,   external:: prd1,prd2
      double precision, external:: logistic2
      double complex, external:: zcosh
      double complex coshthm(Nh)
      double complex a1, a2, b1, b2, w11, w12, w21, w22

      do k = 1, Nh
         coshthm(k) = zcosh(dcmplx(theta_real(k),theta_imag(k))) 
!        write(6,*) k, theta_real(k), theta_imag(k),coshthm(k)
      enddo
      Eloc = 0.D0

      if(spin(Nv).and.spin(1)) then
         Eloc = Eloc + 0.25D0
      elseif(spin(Nv).and.(.not.spin(1))) then
         Eloc = Eloc - 0.25D0 
      elseif((.not.spin(Nv)).and.spin(1)) then
         Eloc = Eloc - 0.25D0
      else
         Eloc = Eloc + 0.25D0
      endif

      if((spin(1)).and.(.not.spin(Nv))) then
          Eloc = Eloc + 0.5D0*exp(-2.0*dcmplx(prmts(1),prmts(Np+1))  &
               + 2.0*dcmplx(prmts(Nv),prmts(Np+Nv)))*prd1(1,Nv,coshthm) 
      elseif((.not.spin(1)).and.(spin(Nv))) then
          Eloc = Eloc + 0.5D0*exp( 2.0*dcmplx(prmts(1),prmts(Np+1)) &
               - 2.0*dcmplx(prmts(Nv),prmts(Np+Nv)))*prd1(Nv,1,coshthm) 
      endif
      
      return
      end function Eloc 

      double complex function prd1(ii,kk,coshthm)

      use parameters
      use spin_hidn, ONLY: theta_real, theta_imag
      implicit none
      integer ii,kk   
      integer k
      double complex coshthm(Nh)
      double complex, external:: zcosh

!     lgs(x) = 1 / ( 1 + exp(-x)) ==> exp(-x) = 1/lgs(x) - 1
!     
!     cosh (x) / cosh (y) = (exp(x) + exp(-x))/(exp(y) + exp(-y))
!                         = exp(x-y)*( 1 + exp (-2x) )/(1 + exp(-2y))
!                         = (1/log(y-x)-1)/log(2x)*log(2y)

      prd1 = 1.D0

      do k = 1, Nh
         prd1 = prd1*zcosh(  dcmplx(                           &
                theta_real(k)-2.0*prmts(   kshift+ii+(k-1)*Nv) &
              + 2.0*prmts(   kshift+kk+(k-1)*Nv),              &
                theta_imag(k)-2.0*prmts(Np+kshift+ii+(k-1)*Nv) &
              + 2.0*prmts(Np+kshift+kk+(k-1)*Nv))              &
              ) / coshthm(k)
      enddo

      return
      end function prd1

!     double complex function prd2(ii,kk,coshthm)

!     use parameters
!     use spin_hidn, ONLY: theta_real, theta_imag
!     implicit none
!     integer ii,kk
!     integer k
!     double complex coshthm(Nh)
!     double complex, external:: zcosh

!     lgs(x) = 1 / ( 1 + exp(-x)) ==> exp(-x) = 1/lgs(x) - 1
!     
!     cosh (x) / cosh (y) = (exp(x) + exp(-x))/(exp(y) + exp(-y))
!                         = exp(x-y)*( 1 + exp (-2x) )/(1 + exp(-2y))
!                         = (1/log(y-x)-1)/log(2x)*log(2y)

!     prd2 = 1.D0

!     do k = 1, Nh
!        prd2 = prd2*zcosh(   dcmplx(                             &
!               theta_real(k) + 2.0*prmts(kshift+ii+(k-1)*Nv)    &
!              -2.0*prmts(kshift+kk+(k-1)*Nv),                   &
!               theta_imag(k) + 2.0*prmts(Np+kshift+ii+(k-1)*Nv) &  
!              -2.0*prmts(Np+kshift+kk+(k-1)*Nv))                &
!               )/coshthm(k)
!     enddo

!     return
!     end function prd2

      double complex function zcosh(z)

      implicit none
      double complex z
      double precision x, y

      ! 2 cosh(x+iy) = exp(x)*(cos(y)+isin(y))+exp(-x)*(cos(y)-isin(y))
      !              = 2cosh(x)*cos(y) + 2isinh(x)*sin(y)

      x = real(z)
      y = aimag(z)

      zcosh =dcmplx(cosh(x)*cos(y), sinh(x)*sin(y))

      return
      end function zcosh

