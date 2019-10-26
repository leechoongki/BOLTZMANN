       double precision function Eloc()
!      double complex function Eloc()

      use parameters
      use spin_hidn

      implicit none
      integer k
      double complex,   external:: prd
      double precision, external:: logistic2
      double complex, external:: zcosh
      double complex coshthm(Nh)

      do k = 1, Nh
         coshthm(k) = zcosh(dcmplx(theta_real(k),theta_imag(k))) 
!        write(6,*) k, theta_real(k), theta_imag(k),coshthm(k)
      enddo

      Eloc = 0.D0
      do k = 1, Nv-1
         if(spin(k).and.spin(k+1)) then
            Eloc = Eloc + 0.25D0
         elseif(spin(k).and.(.not.spin(k+1))) then
            Eloc = Eloc - 0.25D0 
         elseif((.not.spin(k)).and.spin(k+1)) then
            Eloc = Eloc - 0.25D0
         else
            Eloc = Eloc + 0.25D0
         endif
      enddo

      if(spin(Nv).and.spin(1)) then
         Eloc = Eloc + 0.25D0
      elseif(spin(Nv).and.(.not.spin(1))) then
         Eloc = Eloc - 0.25D0 
      elseif((.not.spin(Nv)).and.spin(1)) then
         Eloc = Eloc - 0.25D0
      else
         Eloc = Eloc + 0.25D0
      endif

      do k = 1, Nv-1
         if((spin(k)).and.(.not.spin(k+1))) then
             Eloc = Eloc + 0.5D0*exp(-2.0*dcmplx(prmts(k),prmts(Np+k)) &
                   + 2.0*dcmplx(prmts(k+1),prmts(Np+k+1)))*prd(k,k+1,coshthm)  
         elseif((.not.spin(k)).and.(spin(k+1))) then
             Eloc = Eloc + 0.5D0*exp( 2.0*dcmplx(prmts(k),prmts(Np+k)) &
                   - 2.0*dcmplx(prmts(k+1),prmts(Np+k+1)))*prd(k+1,k,coshthm)
         endif
      enddo

      if((spin(Nv)).and.(.not.spin(1))) then
          Eloc = Eloc + 0.5D0*exp(-2.0*dcmplx(prmts(Nv),prmts(Np+Nv))  &
               + 2.0*dcmplx(prmts(1),prmts(Np+1)))*prd(Nv,1,coshthm) 
      elseif((.not.spin(Nv)).and.(spin(1))) then
          Eloc = Eloc + 0.5D0*exp( 2.0*dcmplx(prmts(Nv),prmts(Np+Nv)) &
               - 2.0*dcmplx(prmts(1),prmts(Np+1)))*prd(1,Nv,coshthm) 
      endif
      
      return
      end function Eloc 

      double complex function prd(ii,kk,coshthm)

      use parameters
      use spin_hidn, ONLY: theta_real, theta_imag
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
         prd = prd*zcosh(  dcmplx(                            &
                theta_real(k)-2.0*prmts(kshift+ii+(k-1)*Nv)    &
              + 2.0*prmts(kshift+kk+(k-1)*Nv),                 &
                theta_imag(k)-2.0*prmts(Np+kshift+ii+(k-1)*Nv) &
              + 2.0*prmts(Np+kshift+kk+(k-1)*Nv))              &
              ) / coshthm(k)
      enddo

      return
      end function prd

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

      double complex function ztanh(z)
 
      implicit none
      double complex z
      double precision x, y

      x = real(z)
      y = aimag(z)

      ztanh = dcmplx(tanh(x),tan(y))/dcmplx(1.D0,tanh(x)*tan(y))
      return
      end function ztanh
   
