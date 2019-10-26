      program AFHeisenberg_vmc

      use lgs
      use spin_hidn
      use parameters

      implicit none
 
      integer i,j,k,iL ! counting index
      real rs1, rs2, cputime1, cputime2

      double complex   egs,egs_tmp ! ground state energy
      double complex   one,ione
      parameter( one=dcmplx(1.0,0.0))
      parameter(ione=dcmplx(0.0,1.0))
      double precision vsum, lambda
      double precision eta,eta0 ! learning rate
      double precision gradnorm, gradnorm_tol
      double precision evariance, egsq
      double precision ss, ss_tmp


      integer, allocatable:: twopwr(:) ! powers of 2
      integer, allocatable:: icount(:)
      integer ii, iisum, ik
      integer MAXITR
      integer, external:: indx
      double precision,allocatable:: distribution(:)
      double precision,allocatable:: xgr(:)
      double precision prd,distsum
      double complex   tantmp
  
!     double precision, external:: coshl,expl 
      logical file_exist

      namelist/input/eta0,Nv,ino,ipa,Nw

      pi = acos(-1.D0)

      Nv = 10
      ino = 1
      ipa = 1
      Nw  = 100000
      MAXITR = 500
      gradnorm_tol = 1.D-3
      L = 100
      Jz = 1.D0
      nexp = 2000000
      eta0 = 0.5
      npi = 0

      read(5,nml=input) 

      Nhn = ino*Nv
      Nhp = ipa*Nv

      kshift_norm = 2*Nv + Nhn + Nhp
      kshift_phas = 2*Nv + Nhn + Nhp + Nhn*Nv 
      Np = 2*Nv+Nhn+Nhp+(Nhp+Nhn)*Nv

      write(6,*) " Nw:",Nw
      write(6,*) "Nhn:",Nhn
      write(6,*) "Nhp:",Nhp
      write(6,*) " Np:",Np

      allocate(exp_save(nexp))

      call explee_save()

      allocate(twopwr(Nv))
      call find_twopwr(Nv,twopwr)

      allocate(  prmts(Np))  ! 
      allocate( vprmts(Np))
      allocate( dprmts(Np))

      allocate( xgr(Np))
      allocate(spin(Nv))
      allocate(theta_norm(Nhn),theta_phas(Nhp))
      allocate(icount(0:twopwr(Nv)-1))
      allocate(distribution(twopwr(Nv)))
      allocate(Dpsi(Np,Nw))

      Dpsi = dcmplx(0.D0,0.D0)

      call init_random_seed()
!     initialize the network

      inquire(file="parameters",exist=file_exist)
      if(file_exist) then
         write(6,*) "#Parameters read from the file "
         open(unit=10,file="parameters",form="unformatted")
              read(10) (prmts(k),k=1,Np)
         close(10)
      else
         do k = 1, Np
            call random_number(rs1)
            prmts(k) = (0.5-rs1)*0.5
        enddo
      endif

      npi = 0

 20   continue

      spin(:) = .false.

      do k = 1, Nv
         call random_number(rs1)
         if(rs1.gt.0.5D0) spin(k) = .true.
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
      egsq = 0.D0
      evariance = 0.D0

      vprmts(:) = dcmplx(0.D0,0.D0)  ! <  D_k  > 
      dprmts(:) = dcmplx(0.D0,0.D0)  ! < H D_k > and \partial_k < H >

      icount(:) = 0

      iL = 0
      ss = 0.D0
      do k = 1, Nw*L
         call mcstep()
         if(mod(k,L).eq.0) then
            iL = iL + 1
            call Eloc(egs_tmp,ss_tmp)

             ss  =  ss + ss_tmp
            egs  = egs + egs_tmp
            egsq = egsq + dconjg(egs_tmp)*egs_tmp

            ! Derivative w.r.t bias for visual layer - norm
            do i = 1, Nv
               Dpsi(i,iL) = 0.5*merge( one, -one,spin(i))
               vprmts(i) = vprmts(i)+Dpsi(i,iL)
               dprmts(i) = dprmts(i)+egs_tmp*dconjg(Dpsi(i,iL))
            enddo

            ! Derivative w.r.t bias for visual layer - phase
            do i = Nv+1, 2*Nv 
                 Dpsi(i,iL) = 0.5*merge(ione,-ione,spin(i-Nv))
               vprmts(i) = vprmts(i)+Dpsi(i,iL)
               dprmts(i) = dprmts(i)+egs_tmp*dconjg(Dpsi(i,iL))
            enddo

            ! Derivative w.r.t bias for hidden layer - norm
            do i = 2*Nv+1, 2*Nv+Nhn
               Dpsi(i,iL) = tanh(theta_norm(i-2*Nv))*0.5
               vprmts(i)  = vprmts(i) + Dpsi(i,iL)
               dprmts(i)  = dprmts(i) + egs_tmp*dconjg(Dpsi(i,iL))
            enddo

            ! Derivative w.r.t bias for hidden layer - phase
            do i = 2*Nv+Nhn+1, 2*Nv+Nhn+Nhp
               Dpsi(i,iL) = tanh(theta_phas(i-2*Nv-Nhn))*ione*0.5
               vprmts(i)  = vprmts(i) + Dpsi(i,iL)
               dprmts(i)  = dprmts(i) + egs_tmp*dconjg(Dpsi(i,iL))
            enddo

            ! Derivative w.r.t link for norm
            do j = 1, Nhn
               tantmp = tanh(theta_norm(j))*0.5
               do i = 1, Nv
                  ik = kshift_norm+i+(j-1)*Nv
                  Dpsi(ik,iL) = merge(tantmp,-tantmp,spin(i))
                  vprmts(ik) = vprmts(ik)+Dpsi(ik,iL)
                  dprmts(ik) = dprmts(ik)+egs_tmp*dconjg(Dpsi(ik,iL))
               enddo
            enddo

            ! Derivative w.r.t link for phase
            do j = 1, Nhp
               tantmp = tanh(theta_phas(j))*ione*0.5
               do i = 1, Nv
                  ik = kshift_phas+i+(j-1)*Nv
                  Dpsi(ik,iL) = merge(tantmp,-tantmp,spin(i))
                  vprmts(ik) = vprmts(ik)+Dpsi(ik,iL)
                  dprmts(ik) = dprmts(ik)+egs_tmp*dconjg(Dpsi(ik,iL))
               enddo
            enddo
            if(ik.ne.Np) stop "ik =\= Np"
         endif
!        ii = indx(Nv,spin,twopwr)
!        icount(ii) = icount(ii) + 1
      enddo

       ss  =   ss/float(Nw)
      egs  =  egs/float(Nw)
      egsq = egsq/float(Nw)

      evariance = (egsq - dconjg(egs)*egs)
      vprmts(:) = vprmts(:)/float(Nw)
!
!     \partial_k < H > = < H D_k > - < H > < D_k >
!
      dprmts(:) = dprmts(:)/float(Nw)
      dprmts(:) = dprmts(:) - egs*dconjg(vprmts(:))

      gradnorm = abs(dot_product(dprmts,dprmts))

      eta = eta0
!     eta = eta0/sqrt(float(npi+1))

      call Sto_Reconf(xgr)

      write(80,*) npi
      do k = 1, Np
         write(80,'(i3,x,f10.5,x,f10.5,x,f10.5,x,f10.5,x,f10.5)') &
         k, prmts(k), vprmts(k), dprmts(k)
      enddo
      write(80,*) " ------------------------------------------------ "
      prmts(:) = prmts(:) - eta*xgr(:)
  
!     lambda = 1.D0 + egs + real(sum(xgr(:)*vprmts(:)))
      call cpu_time(cputime2)

      write(6,21) npi,egs, evariance, gradnorm, ss,&
                  cputime2 - cputime1
 21   format(i5,x,x,2f10.5,x,f10.5,x,f10.5,x,f10.5,x,2f10.5,x,f10.5)
      call flush(6)

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
!        do i = 1, Nhn
!           prd = prd*cosh(theta_norm(i))*2.0
!        enddo

!        ii = ii + 1

!        vsum = 0.D0
!        do i = 1, Nv
!           vsum = vsum+merge(1,-1,spin(i))*prmts(i)
!        enddo
!        distribution(k)=prd*exp(vsum)
!        egs = egs + eloc()*distribution(k)
!     enddo
!
!     distsum = sum(distribution(:))
!     write(6,*) egs/distsum
!     do k = 1, twopwr(Nv)
!        write(200+npi,*) k, distribution(k)/distsum
!     enddo;stop

      open(unit=10,file="parameters",form="unformatted")
          write(10) (prmts(k),k=1,Np)
      close(10)
      rewind(10)

      if((npi.gt.MAXITR).or.(gradnorm.lt.gradnorm_tol)) goto 30 
      goto 20

 30   continue

      deallocate(spin,prmts)
      deallocate(theta_norm,theta_phas)
      deallocate(xgr)
      deallocate(twopwr)
      deallocate(icount)
      deallocate(Dpsi)
 

      end program AFheisenberg_vmc

      
