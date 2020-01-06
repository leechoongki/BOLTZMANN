      program ISING_TE

      use lgs
      use spin_hidn
      use parameters

      implicit none
 
      integer i,j,k,iL ! counting index
      real rs1, rs2, memory_max, cputime1, cputime2

!     double precision egs,egs_tmp ! ground state energy
      double complex   egs,egs_tmp ! ground state energy
      double complex   sxx,sxx_tmp ! < \sigma_x >
!     double precision vsum
      double precision vsum, lambda
      double precision dt
      double complex cmpone
      parameter(cmpone=dcmplx(0.D0,1.D0))
      double complex, external:: ztanh

      double precision, allocatable:: xgr(:)

      integer, allocatable:: twopwr(:) ! powers of 2
      integer, allocatable:: icount(:)
      integer ii, iisum, ik
      integer MAXITR
      integer, external:: indx
      double precision,allocatable :: distribution(:)
      double precision prd,distsum
  
      double complex, external:: zcosh 
      namelist/input/hfield,dt

      pi = acos(-1.D0)

      Nv = 10
      Nh = 10
      Nw  = 100000
      MAXITR = 500
      L = 100
      hfield = 4.0D0
      nexp = 1000000
      ncos = 1000000
      dt = 0.001D0
      kshift = 2*(Nv + Nh)
      Np = kshift/2 + Nv*Nh
      npi = 0

      read(5,nml=input) 

      allocate(exp_save(nexp))
      allocate(cos_save(ncos))
      allocate(sin_save(ncos))

      call explee_save()

      allocate(twopwr(Nv))
      call find_twopwr(Nv,twopwr)

      allocate(  prmts(2*Np))  ! 
      allocate( vprmts(2*Np))
      allocate( dprmts(2*Np))
      allocate(    xgr(2*Np))
      allocate(spin(Nv), hidn(Nh))
      allocate(theta(Nh))
      allocate(icount(0:twopwr(Nv)-1))
      allocate(distribution(twopwr(Nv)))
      allocate(Dpsi(2*Np,Nw))

      Dpsi = dcmplx(0.D0,0.D0)

      call init_random_seed()
!     initialize the network

      open(unit=10,file="parameters",form="unformatted")   
      read(10) (prmts(2*k-1),k=1,Np)
      close(10)
      open(unit=11,file="parameters.read",form="formatted")
      do k = 1, Np
         write(11,*) prmts(2*k-1), prmts(2*k)
      enddo
      close(11)

      npi = 0

 20   continue

      spin(:) = .false.
      hidn(:) = .false.
      icount(:) = 0

      do k = 1, Nv
         call random_number(rs1)
         if(rs1.gt.0.5D0) spin(k) = .true.
      enddo

      do k = 1, Nh
         call random_number(rs1)
         if(rs1.gt.0.5D0) hidn(k) = .true.
      enddo

      call find_theta()
!           1 - Nv : a(i)
!     Nv+1 - Nv+Nh : b(i)
!     Nv+Nh+1 ~ Nv+Nh+Nv*Nh : w(i,j) - visual layer first 
!     prmts - visual layer first

      call cpu_time(cputime1)
      do k = 1, Nw
         call mcstep()
      enddo

      egs = dcmplx(0.D0,0.D0)
      sxx = dcmplx(0.D0,0.D0)

      vprmts(:) = dcmplx(0.D0,0.D0)  ! <  D_k  > 
      dprmts(:) = dcmplx(0.D0,0.D0)  ! < H D_k > and \partial_k < H >

      icount(:) = 0

      iL = 0
      do k = 1, Nw*L
         call mcstep()
         if(mod(k,L).eq.0) then
            iL = iL + 1
            call Eloc(egs_tmp,sxx_tmp)
            egs  = egs + egs_tmp
            sxx  = sxx + sxx_tmp

            do i = 1, Nv
              if(spin(i)) then
                 Dpsi(2*i-1,iL) = dcmplx( 1.D0,0.D0)
                 Dpsi(2*i  ,iL) = dcmplx( 0.D0,1.D0)
              else
                 Dpsi(2*i-1,iL) = dcmplx(-1.D0, 0.D0)
                 Dpsi(2*i  ,iL) = dcmplx( 0.D0,-1.D0)
              endif 
              vprmts(2*i-1)  = vprmts(2*i-1) + Dpsi(2*i-1,iL)
              vprmts(2*i  )  = vprmts(2*i  ) + Dpsi(2*i  ,iL)

              dprmts(2*i-1)  = dprmts(2*i-1) + egs_tmp*dconjg(Dpsi(2*i-1,iL))
              dprmts(2*i  )  = dprmts(2*i  ) + egs_tmp*dconjg(Dpsi(2*i,iL))
            enddo

            do i = Nv+1, Nv+Nh
               Dpsi(2*i-1,iL) = ztanh(theta(i-Nv))
               Dpsi(2*i  ,iL) = ztanh(theta(i-Nv))*cmpone
               vprmts(2*i-1)  = vprmts(2*i-1) + Dpsi(2*i-1,iL)
               vprmts(2*i)    = vprmts(2*i  ) + Dpsi(2*i,iL)
               dprmts(2*i-1)  = dprmts(2*i-1) + egs_tmp*dconjg(Dpsi(2*i-1,iL))
               dprmts(2*i)    = dprmts(2*i)   + egs_tmp*dconjg(Dpsi(2*i,iL))
            enddo

            do i = 1, Nh
               do j = 1, Nv
                  ik = kshift+2*j-1+(i-1)*Nv*2
                  if(spin(j)) then
                     Dpsi(ik  ,iL) = ztanh(theta(i))
                     Dpsi(ik+1,iL) = ztanh(theta(i))*cmpone
                  else
                     Dpsi(ik  ,iL) =-ztanh(theta(i))
                     Dpsi(ik+1,iL) =-ztanh(theta(i))*cmpone
                  endif
                  vprmts(ik  )    = vprmts(ik  ) + Dpsi(ik,  iL)
                  vprmts(ik+1)    = vprmts(ik+1) + Dpsi(ik+1,iL)
                  dprmts(ik  )    = dprmts(ik  ) + egs_tmp*dconjg(Dpsi(ik  ,iL))
                  dprmts(ik+1)    = dprmts(ik+1) + egs_tmp*dconjg(Dpsi(ik+1,iL))
               enddo
            enddo
         endif
         ii = indx(Nv,spin,twopwr)
         icount(ii) = icount(ii) + 1
      enddo

      sxx = sxx/float(Nw)
      egs = egs/float(Nw)

      iisum = sum(icount(:))
!     write(6,*) egs, sxx

      vprmts(:) = vprmts(:)/float(Nw)
      dprmts(:) = dprmts(:)/float(Nw)
      dprmts(:) = dprmts(:) - egs*dconjg(vprmts(:))

      call Sto_Reconf(xgr)
      
      do k = 1, Np
         prmts(2*k-1) = prmts(2*k-1)+dt*xgr(2*k  )
         prmts(2*k  ) = prmts(2*k  )-dt*xgr(2*k-1)
      enddo

      call cpu_time(cputime2)

      write(6,'(3f10.5)') (npi-1)*dt, sxx

      if((npi.gt.MAXITR)) goto 30 
      goto 20

 30   continue

      deallocate(spin,hidn,prmts)
      deallocate(theta)
      deallocate(xgr)
      deallocate(twopwr)
      deallocate(icount)
      deallocate(Dpsi)
 

      end program ISING_TE
