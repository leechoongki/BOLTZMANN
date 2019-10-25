      module lgs

      double precision dxx
      integer          nlgs
      double precision xmax
      double precision, allocatable:: logisave(:)
 
      end module lgs 


      program test

      use lgs

      implicit none
      integer k
      double precision x, dx
      double precision t1, t2
      double precision, external:: logistics 

      nlgs = 1000000
      allocate(logisave(nlgs))

      call logistic_save()

      dx = 0.037
      do k = 1, 1000
         x = -10.D0 + (k-1)*dx
         write(6,*) x, logistics(x), 1.D0/(1.D0+exp(-x))
      enddo

      end program test



