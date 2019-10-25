      subroutine find_twopwr(n,twopwr)

      implicit none
      integer n
      integer twopwr(n)
      integer k
      
      twopwr(1) = 2
      do k = 2, n
         twopwr(k) = 2*twopwr(k-1) 
      enddo

      return 
      end subroutine find_twopwr

      integer function indx(n,spin,twopwr)
      implicit none
      integer n, k
      integer spin(n)
      integer twopwr(n)

      indx = 0 
      if(spin(1).eq.1) indx = 1
      do k = 2, n 
         if(spin(k).eq.1) then        
            indx = indx + twopwr(k-1)
         endif
      enddo 

      return
      end function indx

