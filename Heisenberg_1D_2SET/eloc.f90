      subroutine Eloc(egs,ss_tmp)

      use parameters
      use spin_hidn

      implicit none
      integer k
      double complex   egs
      double precision ss_tmp
      double precision cosh_norm(Nhn)
      double precision cosh_phas(Nhp)
      double precision,external::prd_norm,prd_phas

      do k = 1, Nhn
         cosh_norm(k) = cosh(theta_norm(k)) 
      enddo
      do k = 1, Nhp
         cosh_phas(k) = cosh(theta_phas(k)) 
      enddo

      egs = dcmplx(0.D0,0.D0)
      ss_tmp = 0.D0

      do k = 1, Nv-1
         if(spin(k).and.spin(k+1)) then
            ss_tmp = ss_tmp+1.0
         elseif(spin(k).and.(.not.spin(k+1))) then
            ss_tmp = ss_tmp-1.0
         elseif((.not.spin(k)).and.spin(k+1)) then
            ss_tmp = ss_tmp-1.0
         else
            ss_tmp = ss_tmp+1.0
         endif
      enddo

      if(spin(Nv).and.spin(1)) then
         ss_tmp = ss_tmp+1.0
      elseif(spin(Nv).and.(.not.spin(1))) then
         ss_tmp = ss_tmp-1.0
      elseif((.not.spin(Nv)).and.spin(1)) then
         ss_tmp = ss_tmp-1.0
      else
         ss_tmp = ss_tmp+1.0
      endif

      egs = egs + Jz*ss_tmp

      do k = 1, Nv-1
         if((spin(k)).and.(.not.spin(k+1))) then
             egs=egs+2.0*exp(dcmplx(-prmts(k)+prmts(k+1),&
                            -prmts(Nv+k)+prmts(Nv+k+1)   &
                            +prd_phas(k,k+1,cosh_phas))) &
                 *prd_norm(k,k+1,cosh_norm)
         elseif(.not.(spin(k)).and.spin(k+1)) then
             egs=egs+2.0*exp(dcmplx( prmts(k)-prmts(k+1),&
                             prmts(Nv+k)-prmts(Nv+k+1)   & 
                            +prd_phas(k+1,k,cosh_phas))) &
                 *prd_norm(k+1,k,cosh_norm)
         endif
      enddo

      if((spin(Nv)).and.(.not.spin(1))) then
          egs=egs+2.0*exp(dcmplx(-prmts(Nv)+prmts(1),&
                         -prmts(Nv+Nv)+prmts(Nv+1)   &
                         +prd_phas(Nv,1,cosh_phas))) &
              *prd_norm(Nv,1,cosh_norm)
      elseif(.not.(spin(Nv)).and.spin(1)) then
          egs=egs+2.0*exp(dcmplx( prmts(Nv)-prmts(1),          &
                          prmts(Nv+Nv)-prmts(Nv+1)             &
                        + prd_phas(1,Nv,cosh_phas)))           &
              *prd_norm(1,Nv,cosh_norm)
      endif

      return
      end subroutine Eloc 

      double precision function prd_norm(ii,kk,coshthm)

      use parameters
      use spin_hidn
      implicit none
      integer ii,kk   ! ii - spin up, kk - spin down
      integer k
      double precision coshthm(Nhn)

      prd_norm=1.D0
      do k = 1, Nhn
         prd_norm=prd_norm                                             &
                 *cosh(theta_norm(k)-2.0*(prmts(kshift_norm+ii+(k-1)*Nv)&
                 -prmts(kshift_norm+kk+(k-1)*Nv)))/coshthm(k)
      enddo
      prd_norm=dsqrt(prd_norm)

      return
      end function prd_norm


      double precision function prd_phas(ii,kk,coshthm)

      use parameters
      use spin_hidn

      implicit none
      integer ii, kk
      integer k
      double precision coshthm(Nhp)
 
      prd_phas=1.0
      do k = 1, Nhp
         prd_phas=prd_phas                                            &
              *cosh(theta_phas(k)-2.0*(prmts(kshift_phas+ii+(k-1)*Nv) & 
                   -prmts(kshift_phas+kk+(k-1)*Nv)))/coshthm(k)       
      enddo            
      prd_phas = 0.5*log(prd_phas)
      
      return
      end function prd_phas
