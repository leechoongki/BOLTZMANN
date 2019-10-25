      program ISING_ED

!--------------------------------------------------------------
!  
!  Exact diagnolization code for transverse field Ising model 
!  on one dimensional chain with periodic bounday condition.
!  You may need LAPACK to build executable file.
!
!                                                    2018.5.24
!
!                  Author: Lee CHOONGKI (leechoongki@gmail.com)
!
!---------------------------------------------------------------

      implicit none

      integer i,k
      integer nsite, nba
      logical, allocatable:: ibasis(:,:)
      integer, allocatable:: twopwr(:)
      double precision, allocatable:: Ha(:,:)
      double precision  hfield
      integer, external:: idx

      write(6,*)
      write(6,'(a)') "This program calculates the energy spectrum"
      write(6,'(a)') "of transverse field Ising model - cklee    "
      write(6,*)

      write(6,'(a)') "Enter number of sites - "
      read(5,*) nsite
      write(6,'(a)') "Enter the value of hfield - "
      read(5,*) hfield
      write(6,'(a)') "Calculating ..."

      allocate(twopwr(0:nsite))
      call find_twopwr(nsite,twopwr)
      nba = twopwr(nsite)
      allocate(ibasis(nsite,nba))

      call basis(nba,nsite,ibasis)
!     do k = 1, nba
!        write(6,*) (ibasis(i,k),i=1,nsite),idx(nsite,ibasis(:,k),twopwr)
!     enddo;stop

      allocate(Ha(nba,nba)) 

      call set_hamiltonian(hfield,nba,nsite,Ha,ibasis,twopwr)

      call lapack_diag(Ha,nba)

!     do k = 1, nba
!        write(6,*) Ha(k,k),(ibasis(i,k),i=1,nsite)
!     enddo

      end program ISING_ED

      subroutine basis(nba,nsite,ibasis)

      implicit none
      integer ii,i,k,nba,nsite
      logical ibasis(nsite,nba)

      ii = 0
      do k = 1, nba
         ibasis(:,k) = .false.
         do i = 0, nsite-1
            if(Btest(ii,i)) ibasis(i+1,k) = .true.
         enddo
         ii = ii + 1
      enddo

      return
      end subroutine basis

      subroutine find_twopwr(n,twopwr)

      implicit none
      integer n
      integer twopwr(0:n)
      integer k
      
      twopwr(0) = 1
      twopwr(1) = 2
      do k = 2, n
         twopwr(k) = 2*twopwr(k-1) 
      enddo

      return 
      end subroutine find_twopwr

      subroutine set_hamiltonian(hfield,nba,nsite,Ha,ibasis,twopwr)

      implicit none
      double precision hfield
      integer nba,nsite
      double precision Ha(nba,nba)
      logical ibasis(nsite,nba)
      logical iba_tmp(nsite)
      integer i,k,twopwr(0:nsite)
      double precision, external:: spin_gob
      integer, external:: idx

      Ha(:,:) = 0.D0
      do k = 1, nba
         do i = 1, nsite-1
            Ha(k,k) = Ha(k,k) - spin_gob(ibasis(i,k),ibasis(i+1,k))
!           write(6,*) ibasis(i,k), ibasis(i+1,k), spin_gob(ibasis(i,k),ibasis(i+1,k))
         enddo 
         Ha(k,k) = Ha(k,k) - spin_gob(ibasis(nsite,k),ibasis(1,k))
      enddo

      do k = 1, nba
         iba_tmp(:) = ibasis(:,k)
         do i = 1, nsite
            if(ibasis(i,k)) then
               iba_tmp(i) = .false.
               Ha(k,idx(nsite,iba_tmp,twopwr)) = &
                  Ha(k,idx(nsite,iba_tmp,twopwr))-hfield/2.0
               iba_tmp(i) = .true.
            endif
         enddo 

         do i = 1, nsite
            if(.not.ibasis(i,k)) then 
               iba_tmp(i) = .true.
               Ha(k,idx(nsite,iba_tmp,twopwr)) = &
                  Ha(k,idx(nsite,iba_tmp,twopwr))-hfield/2.0
               iba_tmp(i) = .false.
            endif
         enddo 
      enddo

      do k = 1, nba
         do i = 1, nba
            if(Ha(i,k).ne.Ha(k,i)) then
               write(6,*) i,k,Ha(i,k), Ha(k,i)
            endif
         enddo
      enddo

      return
      end subroutine set_hamiltonian

      double precision function spin_gob(f1,f2)

      logical f1,f2

      if(f1) then
         if(f2) spin_gob = 0.25
         if(.not.f2) spin_gob = -0.25
      else
         if(f2) spin_gob = -0.25
         if(.not.f2) spin_gob = 0.25
      endif
      
      return
      end function spin_gob 

      integer function idx(nsite,iba,twopwr)
      implicit none
      integer nsite, twopwr(0:nsite)
      logical iba(nsite)
      integer k

      idx = 1 
      do k = 1, nsite
         if(iba(k)) idx = idx + twopwr(k-1)
      enddo

      return
      end

      subroutine lapack_diag(Ha,nc)

      implicit none
      
      integer nc
      double precision Ha(nc,nc), elv(nc), p0(nc,nc)
      integer i,j,k

! ########################################################
 
      double precision A(nc,nc),ptmp(nc)
      integer lda, info
      double precision work(3*nc+10)

      lda = nc
      A(:,:) = Ha(:,:)

      call dsyev('v','u',nc,A,lda,elv,work,3*nc+10,info) 
      p0(:,:) = A(:,:)

      write(6,'(a)') "10 eigenvalues"
      do k = 1, 10
         write(6,'(i2,2x,f10.5)') k, elv(k)
      enddo

      do i = 1,nc
         ptmp(i) = sum(Ha(i,:)*p0(:,1))
      enddo

      write(6,*)
      write(6,'(a,x,f10.5)') "Ground state energy:", elv(1)

      return
      end
 
