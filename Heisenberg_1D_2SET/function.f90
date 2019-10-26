      subroutine find_theta()

      use parameters
      use spin_hidn
!!    theta(j) = \sum_i wght(j,i)*spin(i) + b(j)
!
      implicit none
      integer k,j,ii,nn
!     
      nn = 2*Nv+Nhn

      do j = 1, Nhn
         theta_norm(j) = 0.0
         ii = (j-1)*Nv
         do k = 1, Nv
            theta_norm(j)=theta_norm(j) &
                         +merge(1,-1,spin(k))*prmts(kshift_norm+k+ii)
         enddo
         theta_norm(j)=theta_norm(j)+prmts(2*Nv+j)
      enddo

      do j = 1, Nhp
         theta_phas(j) = 0.0
         ii = (j-1)*Nv
         do k = 1, Nv
            theta_phas(j)=theta_phas(j) &
                         +merge(1,-1,spin(k))*prmts(kshift_phas+k+ii)
         enddo
         theta_phas(j)=theta_phas(j)+prmts(nn+j)
      enddo

      return
      end subroutine find_theta

      double precision function expl(x)
 
      use lgs
      implicit none
      integer n, k
      double precision x

      k = int((x+xmax)/dxx) + 1
      if(k.lt.0.or.k.gt.nexp) then
         expl = exp(x)
         return
      else
         expl = exp_save(k)
         return
      endif

      end function expl

      subroutine explee_save()

      use lgs
      implicit none
      integer k
     
      dxx = 2.D0*xmax/dfloat(nexp)

      do k = 1, nexp
         exp_save(k) = exp( -xmax + (k-1)*dxx )      
      enddo

      return
      end subroutine explee_save

      double precision function coshl(x)

      use lgs
      implicit none
      double precision x
      double precision, external:: expl

      if(abs(x).gt.xmax) then
         coshl = cosh(x)
         return
      else
         coshl = (expl(x) + expl(-x))/2.D0
         return
      endif

      return
      end

      double precision function sinhl(x)

      use lgs
      implicit none
      double precision x
      double precision, external:: expl

      if(abs(x).gt.xmax) then
         sinhl = sinh(x)
         return
      else
         sinhl = (expl(x) - expl(-x))/2.D0
         return
      endif

      return
      end

