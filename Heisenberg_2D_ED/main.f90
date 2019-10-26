      program Heisenberg_ED

      implicit none

      integer i,k
      integer nsite, nba, nx, ny,nz
      logical, allocatable:: ibasis(:,:)
      integer, allocatable:: twopwr(:)
      double precision, allocatable:: Ha(:,:)
      double precision Jz, Jxy
      integer, external:: idx
      integer, allocatable:: innd(:,:)

      nx = 2
      ny = 2
      nsite = nx*ny
      nz = 4
      Jz = 1.D0
      Jxy = 1.0D0

      allocate(innd(nsite,nz))
      allocate(twopwr(0:nsite))
      call find_twopwr(nsite,twopwr)
      nba = twopwr(nsite)
      allocate(ibasis(nsite,nba))

      call find_nn(nsite,nx,ny,innd)
      call basis(nba,nsite,ibasis)

!     do k = 1, nba
!        write(6,*) (ibasis(i,k),i=1,nsite),idx(nsite,ibasis(:,k),twopwr)
!     enddo;stop

      allocate(Ha(nba,nba)) 

      call set_hamiltonian(Jz,Jxy,nba,nsite,Ha,ibasis,twopwr,innd,nz)

      call lapack_diag(Ha,nba)

!     do k = 1, nba
!        write(6,*) Ha(k,k),(ibasis(i,k),i=1,nsite)
!     enddo

      end program HEISENBERG_ED

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

      subroutine set_hamiltonian(Jz,Jxy,nba,nsite,Ha,ibasis,twopwr,innd,nz)

      implicit none
      double precision Jz,Jxy
      integer nba,nsite,nz
      double precision Ha(nba,nba)
      logical ibasis(nsite,nba)
      logical iba_tmp(nsite)
      integer i,k,j,m,twopwr(0:nsite),isum,innd(nsite,nz)
      integer, external::  idx
      double precision, external:: spin_gob

      Ha(:,:) = 0.D0

      do k = 1, nba
         do i = 1, nsite
            do j = 1, nz
               Ha(k,k) = Ha(k,k) + spin_gob(ibasis(i,k),ibasis(innd(i,j),k))
            enddo
         enddo 
         Ha(k,k) = Jz*Ha(k,k)*0.5D0
      enddo

      do k = 1, nba
         iba_tmp(:)=ibasis(:,k) 
         do i = 1, nsite
            if(ibasis(i,k)) then
               iba_tmp(i) = .false.
               do j = 1, nz
                  if(.not.ibasis(innd(i,j),k)) then
                     iba_tmp(innd(i,j)) = .true. 
!                    write(6,*) "+",i,innd(i,j),(ibasis(m,k),m=1,nsite),k,(iba_tmp(m),m=1,nsite),idx(nsite,iba_tmp,twopwr)
                     Ha(k,idx(nsite,iba_tmp,twopwr)) = &
                     Ha(k,idx(nsite,iba_tmp,twopwr)) + Jxy*0.25D0
                     iba_tmp(innd(i,j)) = .false. 
                  endif  
               enddo
               iba_tmp(i) = .true.
            else
               iba_tmp(i) = .true.
               do j = 1, nz
                  if(ibasis(innd(i,j),k)) then
                     iba_tmp(innd(i,j)) = .false.
!                    write(6,*) "-",i,innd(i,j),(ibasis(m,k),m=1,nsite),k,(iba_tmp(m),m=1,nsite),idx(nsite,iba_tmp,twopwr)
                     Ha(k,idx(nsite,iba_tmp,twopwr)) = &
                     Ha(k,idx(nsite,iba_tmp,twopwr)) + Jxy*0.25D0
                     iba_tmp(innd(i,j)) = .true.
                  endif
               enddo
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

! ########################################################

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

! ########################################################

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

! ########################################################

      subroutine lapack_diag(Ha,nc)

      implicit none
      
      integer nc
      double precision Ha(nc,nc), elv(nc), p0(nc,nc)
      integer i,j,k
 
      double precision A(nc,nc),ptmp(nc)
      integer lda, info
      double precision work(3*nc+10)

      lda = nc
      A(:,:) = Ha(:,:)

      call dsyev('v','u',nc,A,lda,elv,work,3*nc+10,info) 
      p0(:,:) = A(:,:)

      do k = 1, 10
         write(6,*) elv(k)
      enddo

      write(6,*) "CHECK EIGENVECTOR"

      do i = 1,nc
         ptmp(i) = sum(Ha(i,:)*p0(:,1))
      enddo
      write(6,*) sum(ptmp(:)*p0(:,1))

      do i = 1, nc 
         do k = 1, 10
            write(10+k,'(f10.5)') p0(i,k)
         enddo
      enddo

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
