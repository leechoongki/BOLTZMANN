      program Heisenberg_1D_ED

!
!     CHOONGKI LEE
!
!     Exact Diagonalization : Heisenberg chain 
!

      implicit none

      integer i,k
      integer nsite, nba
      logical, allocatable:: ibasis(:,:)
      integer, allocatable:: twopwr(:)
      double precision, allocatable:: Ha(:,:)
      integer, external:: idx

      nsite = 5
      allocate(twopwr(0:nsite))
      call find_twopwr(nsite,twopwr)
      nba = twopwr(nsite)
      allocate(ibasis(nsite,nba))

      call basis(nba,nsite,ibasis)
!     do k = 1, nba
!        write(6,*) (ibasis(i,k),i=1,nsite),idx(nsite,ibasis(:,k),twopwr)
!     enddo;stop

      allocate(Ha(nba,nba)) 

      call set_hamiltonian(nba,nsite,Ha,ibasis,twopwr)

      call lapack_diag(Ha,nba,nsite)

!     do k = 1, nba
!        write(6,*) Ha(k,k),(ibasis(i,k),i=1,nsite)
!     enddo

      end program Heisenberg_1D_ED

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

      subroutine set_hamiltonian(nba,nsite,Ha,ibasis,twopwr)

! 
!     H = -J \sum_k ( Sz_k Sz_{k+1} + 0.5*(S+_k*S-_{k+1} + S-_k*S+_{k+1} ) )
!

      implicit none
      integer nba,nsite
      double precision Ha(nba,nba)
      logical ibasis(nsite,nba)
      logical iba_tmp(nsite)
      integer i,k,twopwr(0:nsite)
      integer, external::  idx
      double precision, external:: spin_gob

      Ha(:,:) = 0.D0
      do k = 1, nba
         do i = 1, nsite-1
            Ha(k,k) = Ha(k,k) + spin_gob(ibasis(i,k),ibasis(i+1,k))
         enddo 
         if(nsite.gt.2) then
            Ha(k,k) = Ha(k,k) + spin_gob(ibasis(nsite,k),ibasis(1,k))
         endif
      enddo

      do k = 1, nba
         iba_tmp(:) = ibasis(:,k)
         do i = 1, nsite-1
            if(ibasis(i,k).and.(.not.ibasis(i+1,k))) then
               iba_tmp(i) = .false.
               iba_tmp(i+1) = .true.
               Ha(k,idx(nsite,iba_tmp,twopwr)) = &
                  Ha(k,idx(nsite,iba_tmp,twopwr))+2.0*0.25
               iba_tmp(i) = .true.
               iba_tmp(i+1) = .false.
            endif
         enddo 
         if(nsite.gt.2) then
           if(ibasis(nsite,k).and.(.not.ibasis(1,k))) then
              iba_tmp(nsite) = .false.
              iba_tmp(1) = .true.
              Ha(k,idx(nsite,iba_tmp,twopwr)) = &
              Ha(k,idx(nsite,iba_tmp,twopwr))+2.0*0.25
              iba_tmp(nsite) = .true.
              iba_tmp(1) = .false.
           endif
         endif

         do i = 1, nsite-1
            if((.not.ibasis(i,k)).and.ibasis(i+1,k)) then 
               iba_tmp(i) = .true.
               iba_tmp(i+1) = .false.
               Ha(k,idx(nsite,iba_tmp,twopwr)) = &
                  Ha(k,idx(nsite,iba_tmp,twopwr))+2.0*0.25
               iba_tmp(i) = .false.
               iba_tmp(i+1) = .true.
            endif
         enddo 
         if(nsite.gt.2) then
            if((.not.ibasis(nsite,k)).and.ibasis(1,k)) then
               iba_tmp(nsite) = .true.
               iba_tmp(1) = .false.
               Ha(k,idx(nsite,iba_tmp,twopwr)) = &
                  Ha(k,idx(nsite,iba_tmp,twopwr))+2.0*0.25
               iba_tmp(nsite) = .false.
               iba_tmp(1) = .true.
            endif
         endif
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
         if(f2) spin_gob = 0.25D0
         if(.not.f2) spin_gob = -0.25D0
      else
         if(f2) spin_gob = -0.25D0
         if(.not.f2) spin_gob = 0.25D0
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

      subroutine lapack_diag(Ha,nc,nsite)

      implicit none
      
      integer nc,nsite
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

      write(6,*) elv(1)
      write(6,*) "CHECK EIGENVECTOR"

      do i = 1,nc
         ptmp(i) = sum(Ha(i,:)*p0(:,1))
      enddo
      write(6,*) sum(ptmp(:)*p0(:,1))

      do i = 1, nc 
         write(8,'(f10.5)') p0(i,1)
      enddo
      write(6,*) "1. ENERGY/SITES=",elv(1)/float(nsite)
      write(6,*) "2. ENERGY/SITES=",elv(2)/float(nsite)
      write(6,*) "3. ENERGY/SITES=",elv(3)/float(nsite)
      write(6,*) "4. ENERGY/SITES=",elv(4)/float(nsite)
      write(6,*) "5. ENERGY/SITES=",elv(5)/float(nsite)
      write(6,*) "6. ENERGY/SITES=",elv(6)/float(nsite)
      write(6,*) "7. ENERGY/SITES=",elv(7)/float(nsite)
      write(6,*) "8. ENERGY/SITES=",elv(8)/float(nsite)
      write(6,*) "9. ENERGY/SITES=",elv(9)/float(nsite)
      write(6,*) "9. ENERGY/SITES=",elv(10)/float(nsite)
      write(6,*) "9. ENERGY/SITES=",elv(11)/float(nsite)
      write(6,*) "9. ENERGY/SITES=",elv(12)/float(nsite)
      write(6,*) "9. ENERGY/SITES=",elv(13)/float(nsite)
      write(6,*) "9. ENERGY/SITES=",elv(14)/float(nsite)
      write(6,*) "9. ENERGY/SITES=",elv(15)/float(nsite)
      write(6,*) "9. ENERGY/SITES=",elv(16)/float(nsite)
      write(6,*) "9. ENERGY/SITES=",elv(17)/float(nsite)
      write(6,*) "9. ENERGY/SITES=",elv(18)/float(nsite)
      write(6,*) "9. ENERGY/SITES=",elv(19)/float(nsite)
      write(6,*) "9. ENERGY/SITES=",elv(20)/float(nsite)
      write(6,*) "9. ENERGY/SITES=",elv(21)/float(nsite)

      return
      end
 
