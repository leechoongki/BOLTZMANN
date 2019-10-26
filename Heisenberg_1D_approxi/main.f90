      program ising_vmc

      use lgs
      use spin_hidn
      use parameters

      implicit none
 
      integer i,j,k,iL ! counting index
      real rs1, rs2, memory_max, cputime1, cputime2

!     double precision egs,egs_tmp ! ground state energy
      double complex   egs,egs_tmp ! ground state energy
!     double precision vsum
      double complex   vsum, lambda
      double precision eta,eta0 ! learning rate
      double precision gradnorm, gradnorm_tol
      double precision evariance, egsq
      double complex cmpone
      parameter(cmpone=dcmplx(0.D0,1.D0))
      double complex, external:: Eloc
      double complex, external:: ztanh

      double complex,allocatable:: xgr(:)

      integer, allocatable:: twopwr(:) ! powers of 2
      integer, allocatable:: icount(:)
      integer ii, iisum, ik
      integer MAXITR
      integer, external:: indx
      double precision,allocatable :: distribution(:)
      double precision prd,distsum
  
      double complex, external:: zcosh 

      pi = acos(-1.D0)

      Nv = 5
      Nh = 5
      Nw  = 100000
      Nwr = 50000*2
      MAXITR = 500
      gradnorm_tol = 1.D-3
      L = 100
      Jex = 1.D0
      hfield = 0.0D0
      nexp = 1000000
      ncos = 1000000
      eta0 = 0.5
      kshift = Nv + Nh
      Np = kshift + Nv*Nh
      npi = 0


      write(9,*) 
      write(9,'(a,3x,i8)') "     # of sites at visible layer : ", Nv
      write(9,'(a,3x,i8)') "      # of sites at hidden layer : ", Nh
      write(9,'(a,3x,i8)') " parameter dim.(Nv + Nh + Nv*Nh) : ", Np
      write(9,*)
      write(9,*) 
      write(9,'(a,3x,i8)') " # of measurements = ", L
      write(9,'(a,3x,f10.5)') " eta = ", eta0
      write(9,'(a,3x,f10.5,x,f10.5)') " Jex, Hfield       = ", Jex, Hfield
      write(9,'(a,3x,i8)') " # of grids for exp ftn  = ", nexp
      write(9,'(a,3x,i8)') " # of grids for cos ftn  = ", ncos

      memory_max = 8.0*Nw*(Nv*Nh + Nv + Nh)
      write(9,'(a,x,f10.5,x,a)') "MAXIMUM memory size:",memory_max/1024.0/1024.0,"MB"
      flush(9)

      allocate(exp_save(nexp))
      allocate(cos_save(ncos))
      allocate(sin_save(ncos))

      call explee_save()
      write(6,*) "MAKING EXP, SIN, COS, DB DONE..."

      allocate(twopwr(Nv))
      call find_twopwr(Nv,twopwr)

      allocate(  prmts(Np))  ! 

      allocate( vprmts(Np))
      allocate( dprmts(Np))

      allocate(    xgr(Np))
      allocate(spin(Nv), hidn(Nh))
      allocate(theta(Nh))
      allocate(icount(0:twopwr(Nv)-1))
      allocate(distribution(twopwr(Nv)))
      allocate(Dpsi(Np,Nw))

      Dpsi = dcmplx(0.D0,0.D0)

      call init_random_seed()
!     initialize the network

      do k = 1, Np
         call random_number(rs1)
         call random_number(rs2)
         prmts(k) = dcmplx(0.5-rs1,0.5-rs2)
      enddo

      npi = 0

      eta = 0.3 

 20   continue

!     Nwr = Nwr + npi*100

      spin(:) = .false.
      hidn(:) = .false.

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
      do k = 1, Nwr
         call mcstep()
      enddo

      egs = dcmplx(0.D0,0.D0)
      egsq = 0.D0
      evariance = 0.D0

      vprmts(:) = dcmplx(0.D0,0.D0)  ! <  D_k  > 
      dprmts(:) = dcmplx(0.D0,0.D0)  ! < H D_k > and \partial_k < H >

      icount(:) = 0

      iL = 0
      do k = 1, Nwr*L
         call mcstep()
         if(mod(k,L).eq.0) then
            iL = iL + 1
            egs_tmp = Eloc()
            egs  = egs + egs_tmp
            egsq = egsq + dconjg(egs_tmp)*egs_tmp

            do i = 1, Nv
              if(spin(i)) then
                 Dpsi(i,iL) = dcmplx( 1.D0,0.D0)
              else
                 Dpsi(i,iL) = dcmplx(-1.D0, 0.D0)
              endif 
              vprmts(i)    = vprmts(i) + Dpsi(i,iL)
              dprmts(i)    = dprmts(i) + egs_tmp*dconjg(Dpsi(i,iL))
            enddo

            do i = Nv+1, Nv+Nh
              Dpsi(i,   iL) = ztanh(theta(i-Nv))
              vprmts(i)     = vprmts(i) + Dpsi(i,iL)
              dprmts(i)     = dprmts(i) + egs_tmp*dconjg(Dpsi(i,iL))
            enddo
            do i = 1, Nh
               do j = 1, Nv
                  ik = kshift+j+(i-1)*Nv
                  if(spin(j)) then
                     Dpsi(ik,iL) = ztanh(theta(i))
                  else
                     Dpsi(ik,iL) =-ztanh(theta(i))
                  endif
                  vprmts(ik)    = vprmts(ik) + Dpsi(ik,iL)
                  dprmts(ik)    = dprmts(ik) + egs_tmp*dconjg(Dpsi(ik,iL))
               enddo
            enddo
         endif
!        ii = indx(Nv,spin,twopwr)
!        icount(ii) = icount(ii) + 1
      enddo


      egs = egs/float(Nwr)
      egsq = egsq/float(Nwr)
      evariance = (egsq - dconjg(egs)*egs)

!     write(6,*) npi,"Egs = ", egs, "(", regul_factor, ")"
!     write(6,*) npi,"Egs = ", egs

      vprmts(:) = vprmts(:)/float(Nwr)

!
!     \partial_k < H > = < H D_k > - < H > < D_k >
!
      dprmts(:) = dprmts(:)/float(Nwr)
      dprmts(:) = dprmts(:) - egs*dconjg(vprmts(:))

      gradnorm = abs(dot_product(dprmts,dprmts))

!     prmts(:) = prmts(:) - eta*dprmts(:)

!     if(mod(npi,20).eq.0) then
      eta = eta0/sqrt(float(npi+1))
!     endif

!     eta = eta0

      call Sto_Reconf(xgr)
      prmts(:) = prmts(:) - eta*xgr(:)

  
      lambda = 1.D0 + egs + sum(xgr(:)*vprmts(:))
      call cpu_time(cputime2)

!     write(6,'(i5,x,a,x,f10.5,x,f10.5,x,f10.5,x,f10.5,x,a,f10.5)') npi,"Egs = ", egs, evariance, gradnorm,"CPU TIME=", cputime2 - cputime1
      write(6,21) npi,"Egs = ", egs, evariance, gradnorm,eta,lambda,&
                  "CPU TIME=", cputime2 - cputime1
 21   format(i5,x,a,x,2f10.5,x,f10.5,x,f10.5,x,f10.5,x,2f10.5,x,a,x,f10.5)

!     iisum = sum(icount(:))
!     do k = 0, twopwr(Nv)-1
!        write(100+npi,*) k+1, icount(k)/float(iisum)
!     enddo

!     ii = 0   

!140  continue
      
!     egs = 0.D0
!     do k = 1, twopwr(Nv)

!        spin(:) = .false.
!        do i = 0, Nv-1
!           if(Btest(ii,i)) spin(i+1) = .true.
!        enddo

!        call find_theta()

!        prd = 1.0
!        do i = 1, Nh
!           prd = prd*(cosh(2.0*real(theta(i)))+cos(2.0*aimag(theta(i))))*2.0
!        enddo

!        ii = ii + 1
!        vsum = 0.D0
!        do i = 1, Nv
!           if(spin(i)) then
!              vsum = vsum + real(prmts(i))
!           else
!              vsum = vsum - real(prmts(i))
!           endif
!        enddo
!        distribution(k)=prd*exp(2.0*vsum)
!!       write(6,*) prd,exp(2.0*vsum),distribution(k)
!        egs = egs + eloc()*distribution(k)
!     enddo
!
!     distsum = sum(distribution(:))
!     write(6,*) egs/distsum
!     do k = 1, twopwr(Nv)
!        write(200+npi,*) k, distribution(k)/distsum
!     enddo

!     if((npi.gt.MAXITR).or.(gradnorm.lt.gradnorm_tol)) goto 30 

      if((npi.gt.MAXITR)) goto 30 
      goto 20

 30   continue

      deallocate(spin,hidn,prmts)
!     deallocate(gmm_real,gmm_imag)
      deallocate(theta)
      deallocate(xgr)
      deallocate(twopwr)
      deallocate(icount)
      deallocate(Dpsi)


      end program ising_vmc

      
