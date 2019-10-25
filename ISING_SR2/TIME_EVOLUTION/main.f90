      program ising_vmc

!-------------------------------------------------------------------------
!     
! This is the fortran code for variational Monte Carlo (VMC) calculation 
! using network representation of quantum state (NQS) for transverse field
! Ising model
! based on the method described in "Science, Vol. 355, Issue 6325, pp. 602-606"
!
!                                                                2018.5.24
!
!                             Author: Lee CHOONGKI (leechoongki@gmail.com)
!
!------------------------------------------------------------------------

      use lgs
      use spin_hidn
      use parameters

      implicit none
 
      integer i,j,k,iL ! counting index
      real rs, memory_max

      double precision egs ! ground state energy
      double precision egs_tmp
      double precision qxx
      double precision eta,eta0 ! learning rate
      double precision gradnorm, gradnorm_tol
      double precision, external:: Eloc, tanhb
      double precision,allocatable:: xgr(:)

      integer, allocatable:: twopwr(:) ! powers of 2
      integer, allocatable:: icount(:)
      integer ii, iisum, ik
      integer itg,MAXITR
      integer, external:: indx
      double precision,allocatable :: distribution(:)
      double precision prd,distsum

      Nv = 10
      Nh = 10
      Nw = 100000
      MAXITR = 500
      gradnorm_tol = 1.D-3
       L = 100
      Jex = 1.D0
      hfield = 4.0D0
      nlgs = 1000000
      eta0 = 0.5
      kshift = Nv + Nh
      Np = kshift + Nv*Nh
      npi = 0

      allocate(logisave(nlgs))
      call logistic_save()

      allocate(twopwr(Nv))
      call find_twopwr(Nv,twopwr)

      allocate(  prmts(Np))
      allocate( vprmts(Np))
      allocate( dprmts(Np))
      allocate(    xgr(Np))
      allocate(spin(Nv), hidn(Nh))
      allocate(gmm(Nv),theta(Nh))
      allocate(icount(0:twopwr(Nv)-1))
      allocate(distribution(twopwr(Nv)))
      allocate(Dpsi(Np,Nw))

      call init_random_seed()
!     initialize the network

      do k = 1, Np
         call random_number(rs)
         prmts(k) = (0.5-rs)
      enddo

!     write(88) (spin(k),k=1,Nv),(hidn(k),k=1,Nh),&
!               (prmts(k),k=1,Np)
 
      itg = 0

 20   continue

      spin(:) = 1
      hidn(:) = 1

      do k = 1, Nv
         call random_number(rs)
         if(rs.gt.0.5D0) spin(k) = -1
      enddo

      do k = 1, Nh
         call random_number(rs)
         if(rs.gt.0.5D0) hidn(k) = -1
      enddo

!           1 - Nv : a(i)
!     Nv+1 - Nv+Nh : b(i)
!     Nv+Nh+1 ~ Nv+Nh+Nv*Nh : w(i,j) - visual layer first 
!     prmts - visual layer first


      do k = 1, Nw
         call mcstep()
      enddo

      egs = 0.D0

      vprmts(:) = 0.D0  ! <  D_k  > 
      dprmts(:) = 0.D0  ! < H D_k > and \partial_k < H >

      icount(:) = 0

      iL = 0
      do k = 1, Nw*L
         call mcstep()
         if(mod(k,L).eq.0) then
            iL = iL + 1
            egs_tmp = Eloc()
            egs = egs + egs_tmp
            do i = 1, Nv
              Dpsi(i,iL) = 0.5D0*spin(i)
               vprmts(i) = vprmts(i) + Dpsi(i,iL)
               dprmts(i) = dprmts(i) + egs_tmp*Dpsi(i,iL)
            enddo
            do i = Nv+1, Nv+Nh
              Dpsi(i,iL) = 0.5D0*tanhb(theta(i-Nv))
               vprmts(i) = vprmts(i) + Dpsi(i,iL)
               dprmts(i) = dprmts(i) + egs_tmp*Dpsi(i,iL)
            enddo
            do i = 1, Nh
               do j = 1, Nv
                  ik = kshift+j+(i-1)*Nv
                  Dpsi(ik,iL) = 0.5D0*tanhb(theta(i))*spin(j) 
                  vprmts(ik)=vprmts(ik) + Dpsi(ik,iL)
                  dprmts(ik)=dprmts(ik) + egs_tmp*Dpsi(ik,iL)
               enddo
            enddo
         endif
         ii = indx(Nv,spin,twopwr)
         icount(ii) = icount(ii) + 1
      enddo


      egs = egs/float(Nw)
!     write(6,*) npi,"Egs = ", egs, "(", regul_factor, ")"

      vprmts(:) = vprmts(:)/float(Nw)

!
!     \partial_k < H > = < H D_k > - < H > < D_k >
!
      dprmts(:) = dprmts(:)/float(Nw)
      dprmts(:) = dprmts(:) - egs*vprmts(:)

      gradnorm = dsqrt(sum(dprmts(:)*dprmts(:)))
      write(6,*) npi,"Egs = ", egs, gradnorm
      call flush(6)

      call Sto_Reconf(xgr)
      prmts(:) = prmts(:) - eta*xgr(:)
!     prmts(:) = prmts(:) - eta*dprmts(:)

!     iisum = sum(icount(:))
!     do k = 0, twopwr(Nv)-1
!        write(100+itg,*) k+1, icount(k)/float(iisum)
!     enddo

!     ii = 0   
!     do k = 1, twopwr(Nv)

!        spin(:) = -1
!        do i = 0, Nv-1
!           if(Btest(ii,i)) spin(i+1) = 1
!        enddo

!        call find_theta()
!        prd = 2.D0*cosh(theta(1))
!        do i = 2, Nh
!           prd = prd*2.D0*cosh(theta(i))
!        enddo
!        ii = ii + 1
!        distribution(k)=prd*exp(sum(prmts(1:Nv)*spin(1:Nv)))
!     enddo
!
!     distsum = sum(distribution(:))
!     do k = 1, twopwr(Nv)
!        write(200+itg,*) k, distribution(k)/distsum
!     enddo

      if((itg.gt.MAXITR).or.(gradnorm.lt.gradnorm_tol)) goto 30 
      itg = itg + 1
      eta = eta0/sqrt(float(itg))
      goto 20

 30   continue

      deallocate(spin,hidn,prmts)
      deallocate(gmm,theta)
      deallocate(twopwr)
      deallocate(icount)
      deallocate(Dpsi)

      end program ising_vmc

      
