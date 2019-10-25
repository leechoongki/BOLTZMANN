      subroutine mcstep(Nh,Nv,wght,a,b,spin,hidn,gmm,theta)

      implicit none

      integer i,k
      integer Nh, Nv
      double precision wght(Nh,Nv)
      double precision a(Nv)
      double precision b(Nh)
      integer spin(Nv), hidn(Nh), spin_tmp(Nv)

      real rs
      double precision theta(Nh), gmm(Nv)
      double precision, external:: logistic2

      call find_gamma(Nv,Nh,wght,a,hidn,gmm)

!     do k = 1, 10
!        call random_number(rs)
!        write(6,*) rs
!     enddo; stop

      do k = 1, Nv
         call random_number(rs)
         if(logistic2(2.D0*gmm(k)).gt.rs) then 
            spin(k) =  1 
         else
            spin(k) = -1
         endif
      enddo
!     write(6,*) (spin(k),k=1,Nv);stop

      call find_theta(Nv,Nh,wght,b,spin,theta)

      do k = 1, Nh
         call random_number(rs)
         if(logistic2(2.D0*theta(k)).gt.rs) then 
            hidn(k) = 1 
         else
            hidn(k) = -1
         endif
      enddo

      return
      end
