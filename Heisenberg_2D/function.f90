      subroutine find_theta()

      use parameters
      use spin_hidn
!!    theta(j) = \sum_i wght(j,i)*spin(i) + b(j)
!
      implicit none
      integer k,i
!     
      do i = 1, Nh 
         theta(i) = dcmplx(0.D0,0.D0)
         do k = 1, Nv
            if(spin(k)) then
               theta(i) = theta(i) + prmts(kshift+k+(i-1)*Nv)
            else
               theta(i) = theta(i) - prmts(kshift+k+(i-1)*Nv)
            endif
         enddo!;stop  "visible layer"
         theta(i) = theta(i) + prmts(Nv+i)
      enddo
!!    write(6,*) (theta(i),i=1,Nh);stop

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

       subroutine find_nn(Nv,Nx,Ny,innd)
       implicit none
       integer Nv,Nx,Ny,nz
       parameter(nz = 4)
       integer innd(Nv,nz)
       integer i,k

       do i = 1, Nv
          innd(i,1) = i + 1      
          if(mod(i,Nx).eq.0) innd(i,1) = i - Nx + 1
          innd(i,2) = i - 1
          if(mod(i-1,Nx).eq.0) innd(i,2) = i + Nx - 1 
          innd(i,3) = i - Nx
          if(i-Nx.le.0) innd(i,3) = i + Nx*Ny - Nx
          innd(i,4) = i + Nx
          if(i+Nx.gt.Nx*Ny) innd(i,4) = i - (Nx*Ny - Nx)
       enddo 

       return
       end subroutine find_nn
