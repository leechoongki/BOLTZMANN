      program minresql_test

      use matrix
      use minresqlpModule

      implicit none
      
      integer n, nout,itn
      double precision, allocatable:: b(:), x(:)
      integer i
      external aprod
      real rs

      n = 100
      shift = 0.0D0
      nout = 10

      allocate(b(n),x(n))

      call setup_matrix(n)

!     b(1) = 1.0
!     b(2) = 2.0
      b(1) =  0.2
      b(2) =  0.1

!     do i = 1, n
!        call random_number(rs)
!        b(i) = rs*1.D0
!        write(6,*) b(i)
!     enddo
 

      open(unit=nout,file="out.dat",form="formatted")
 
      call minresqlp(n=n,Aprod=Aprod,b=b,shift=shift,x=x,nout=nout,&
                     itn=itn) 

      close(10)

      do i = 1, n
         write(6,*) x(i)
      enddo

      write(6,*) "check the solution"
      do i = 1, n
         write(6,*) sum(A(i,:)*x(:))-b(i) 
      enddo
  

      end program minresql_test

      subroutine Aprod(n,x,y)

      use matrix, ONLY: A, ai
  
      implicit none

      integer n
      double precision x(n), y(n)
      double precision gob
      integer i, k
     
!     gob = sum(bi(:)*x(:))

      do i = 1, n
         y(i) = sum(A(i,:)*x(:))
      enddo

      return

      end subroutine Aprod

      subroutine setup_matrix(n)

      use matrix
      implicit none
      integer i,k,n
      real rs

      allocate(A(n,n),ai(n))

      
      do i = 1, n
         call random_number(rs)
         ai(i) = rs
      enddo

      do k = 1, n
         do i = 1, n
            A(i,k) = ai(i)*ai(k) 
!           A(i,k) = float(i+k)
         enddo
      enddo
      do k = 1, n
         A(k,k) = A(k,k)  - 0.01
      enddo
      

!     write(6,*) "Det=",A(1,1)*A(2,2) - A(1,2)*A(2,1)

      return
      end
