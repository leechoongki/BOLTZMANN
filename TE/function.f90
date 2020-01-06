      subroutine find_theta()

      use parameters
      use spin_hidn

!!    theta(j) = \sum_i wght(j,i)*spin(i) + b(j)

      implicit none
      integer k,i

      do i = 1, Nh
         theta(i) = dcmplx(0.D0,0.D0)
         do k = 1, Nv
            if(spin(k)) then
               theta(i) = theta(i) +  &
                       dcmplx(prmts(kshift+2*k-1+(i-1)*2*Nv), &
                              prmts(kshift+2*k  +(i-1)*2*Nv)) 
            else
               theta(i) = theta(i) -  &
                       dcmplx(prmts(kshift+2*k-1+(i-1)*2*Nv), &
                              prmts(kshift+2*k  +(i-1)*2*Nv))
            endif
         enddo
         theta(i) = theta(i) + &
                    dcmplx(prmts(2*Nv+2*i-1),prmts(2*Nv+2*i))
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

      double precision function cosl(x)
 
      use lgs
      implicit none
      integer n, k
      double precision x

      k = int((x+smax)/dcs) + 1
      if(k.lt.0.or.k.gt.ncos) then 
         cosl = cos(x)
         return
      else
         cosl = cos_save(k)
         return
      endif

      return
      end function cosl

      double precision function sinl(x)
 
      use lgs
      implicit none
      integer n, k
      double precision x

      k = int((x+smax)/dcs) + 1
      if(k.lt.0.or.k.gt.ncos) then
         sinl = sin(x)
         return
      else
         sinl = sin_save(k)
         return
      endif

      return
      end function sinl

      subroutine explee_save()

      use lgs
      implicit none
      integer k
     
      dxx = 2.D0*xmax/dfloat(nexp)
      dcs = 2.D0*smax/dfloat(ncos)

      do k = 1, nexp
         exp_save(k) = exp( -xmax + (k-1)*dxx )      
      enddo

      do k = 1, ncos
         cos_save(k) = cos( -smax + (k-1)*dcs )
         sin_save(k) = sin( -smax + (k-1)*dcs )
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

