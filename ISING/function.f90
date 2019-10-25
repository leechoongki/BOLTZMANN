      subroutine find_theta(Nv,Nh,wght,b,spin,theta)

!
!     theta(j) = \sum_i wght(j,i)*spin(i) + b(j)
!
 
      implicit none
      integer Nv, Nh
      integer spin(Nv)
      double precision wght(Nh,Nv)
      double precision b(Nh)
      double precision theta(Nh)
      integer k,i
      
      do i = 1, Nh 
         theta(i) = sum(wght(i,:)*spin(:))+b(i)
!        do k = 1, Nv
!           write(6,*) k, wght(i,k),spin(k)
!        enddo;stop
      enddo
!     write(6,*) (theta(i),i=1,Nh);stop

      return
      end subroutine find_theta

      subroutine find_gamma(Nv,Nh,wght,a,hidn,gmm)
!
!     gmm(j) = \sum_i wght(j,i)*hidn(i) + a(j)
!
 
      implicit none
      integer Nv, Nh
      integer hidn(Nh)
      double precision wght(Nh,Nv)
      double precision a(Nv)
      double precision gmm(Nv)
      integer k,i
      
      do i = 1, Nv 
         gmm(i) = sum(wght(:,i)*hidn(:))+a(i)
      enddo

      return
      end subroutine find_gamma

      double precision function logistic(x)
 
      implicit none
      double precision x

      if(x.lt.-20.D0) then 
        logistic = 0.D0
        return
      elseif(x.gt.20.D0) then
        logistic = 1.D0
        return 
      else
        logistic = 1.D0/ ( 1.D0 + exp(-x))
        return
      endif

      end function logistic


      subroutine logistic_save()

      use lgs
      implicit none
      integer k
     
      xmax = 20.D0
      dxx = 2.D0*xmax/dfloat(nlgs) 

      do k = 1, nlgs
         logisave(k) = 1.D0/(1.D0 + exp( xmax - (k-1)*dxx))      
      enddo

      return
      end subroutine logistic_save

      double precision function logistic2(x)

      use lgs

      implicit none
      integer k
      double precision x
   
      if(x.le.-20.D0) then
        logistic2 = 0.D0
      elseif(x.gt.20.D0) then
        logistic2 = 1.D0
      else
        k = int((x+xmax)/dxx)+1 
        logistic2 = logisave(k)
      endif

      end function logistic2

      double precision function tanhb(x)

      implicit none
      double precision x
      double precision, external:: logistic2
      
      ! tanh x = ( 1 - exp (-2x) ) / ( 1 + exp (-2x) )
      !        = ( 2 - ( 1 + exp (-2x) ) / ( 1 + exp (-2x) )
      !        = 2*lgs(2x) - 1

      tanhb = 2.D0*logistic2(x) - 1.D0

      return
      end function tanhb
