      subroutine mcstep()

      use parameters
      use spin_hidn

      implicit none

      integer i,k
      real rs
      double precision, external:: logistic2

      call find_gamma()

!     It looks that random_number() generates different random numbers
!     when it is called by subroutines ...

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

      call find_theta()

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
