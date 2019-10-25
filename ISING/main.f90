      program ising_vmc

      use lgs
      implicit none

      integer Nv !# of sites in visible layer
      integer Nh !# of sites in hidden  layer
      integer Nw !# of mcstep for warm up moving
      integer L  !  steps for measuring

      integer, allocatable:: spin(:)  ! spin values at visible sites
      integer, allocatable:: hidn(:)  ! spin values at hidden  sites
      integer, allocatable:: pspin(:)
      double precision, allocatable:: wght(:,:) ! connection factor bet. visible and hidden sites.
 
      integer i,j,k ! counting index
      real rs

      double precision,allocatable::   a(:), b(:), c(:) ! bias at visible and hidden sites
      double precision,allocatable:: gmm(:), theta(:)
      double precision,allocatable:: da(:), db(:) ! bias at visible and hidden sites
      double precision,allocatable:: dwght(:,:)
      double precision,allocatable:: vda(:),vdb(:) ! bias at visible and hidden sites
      double precision,allocatable:: vdwght(:,:)

      double precision egs ! ground state energy
      double precision egs_tmp
      double precision qxx
      double precision eta ! learning rate
      double precision, external:: Eloc, tanhb
      double precision Jex, hfield  

      integer, allocatable:: twopwr(:) ! powers of 2
      integer, allocatable:: icount(:)
      integer ii, iisum, ik
      integer itg
      integer, external:: indx
      double precision,allocatable :: distribution(:)
      double precision prd,distsum

      Nv = 10
      Nh = 10
      Nw = 100000
       L = 100
      Jex = 1.D0
      hfield = 1.D0
      nlgs = 1000000
      eta = 0.2


      allocate(logisave(nlgs))
      call logistic_save()

      allocate(twopwr(Nv))
      call find_twopwr(Nv,twopwr)
!     do k = 1, Nv
!        write(6,*) k, twopwr(k)
!     enddo

      allocate(spin(Nv), hidn(Nh))
      allocate(pspin(Nv))
      allocate(wght(Nh,Nv)) ! cautious. dimension is  Nh x Nv
      allocate(a(Nv),b(Nh),c(Nv*Nh))
      allocate(gmm(Nv),theta(Nh))
      allocate(da(Nv),db(Nh))
      allocate(dwght(Nh,Nv))
      allocate(vda(Nv),vdb(Nh))
      allocate(vdwght(Nh,Nv))
      allocate(icount(0:twopwr(Nv)-1))
      allocate(distribution(twopwr(Nv)))


!     call init_random_seed()
!     call generate_seed(iseed(1))
!     call random_seed(put=iseed)

!     initialize the network


      goto 10

      ik = 0
      spin(:) = 1

      do i = 1, Nv
         call random_number(rs)
         if(rs.gt.0.5D0) spin(i) = -1
      enddo

      hidn(:) = 1
      do i = 1, Nh
         call random_number(rs)
         if(rs.gt.0.5D0) hidn(i) = -1
      enddo

      do i = 1, Nv
         call random_number(rs)
         a(i) = (0.5D0-rs)*2.D0
         ik = ik + 1
         write(6,*) ik, a(i)
      enddo

      do i = 1, Nv
         call random_number(rs)
         b(i) = (0.5D0-rs)*2.D0
         ik = ik + 1
         write(6,*) ik, b(i)
      enddo

      do i = 1, Nh
         do k = 1, Nv
            call random_number(rs)
            wght(i,k) = 2.D0*(rs-0.5)
            ik = ik + 1
            write(6,*) ik, wght(i,k)
         enddo
      enddo; stop

 10   continue
 
      read(88) (spin(k),k=1,Nv),(hidn(k),k=1,Nh),&
               (a(k),k=1,Nv),(b(k),k=1,Nh),&
               (c(k),k=1,Nv*Nh)

      do k = 1, Nv          
         do i = 1, Nh
            wght(i,k) = c( k + (i-1)*Nv )
         enddo
      enddo

      do k = 1, 10
         call mcstep(Nh,Nv,wght,a,b,spin,hidn,gmm,theta)
      enddo
!     write(6,*) (spin(k),k=1,Nv); stop
      write(6,*) "Warming up OK"

      itg = 1

 20   continue

      egs = 0.D0

      vda(:) = 0.D0
      vdb(:) = 0.D0
      vdwght(:,:) = 0.D0

      da(:) = 0.D0
      db(:) = 0.D0
      dwght(:,:) = 0.D0

      icount(:) = 0
      do k = 1, Nw*L
         call mcstep(Nh,Nv,wght,a,b,spin,hidn,gmm,theta)
         if(mod(k,L).eq.0) then
            egs_tmp = Eloc(spin,Nv,Nh,wght,theta,hfield,Jex,a,b)    
            egs = egs + egs_tmp
            do i = 1, Nv
                 qxx = 0.5D0*spin(i)
               vda(i) = vda(i) + qxx
                da(i) =  da(i) + egs_tmp*qxx
            enddo
            do i = 1, Nh
                  qxx = 0.5D0*tanhb(theta(i))
               vdb(i) = vdb(i) + qxx
                db(i) =  db(i) + egs_tmp*qxx
            enddo
            do i = 1, Nv
               do j = 1, Nh
                  qxx = 0.5D0*tanhb(theta(j))*spin(i) 
                  vdwght(j,i) = vdwght(j,i) + qxx
                   dwght(j,i) =  dwght(j,i) + egs_tmp*qxx 
               enddo
            enddo
         endif
         ii = indx(Nv,spin,twopwr)
         icount(ii) = icount(ii) + 1
      enddo

      iisum = sum(icount(:))
      do k = 0, twopwr(Nv)-1
         write(100+itg,*) k, icount(k)/float(iisum)
      enddo

      ii = 0   
      do k = 1, twopwr(Nv)
         pspin(:) = -1
         do i = 0, Nv-1
            if(Btest(ii,i)) pspin(i+1) = 1
         enddo
         call find_theta(Nv,Nh,wght,b,pspin,theta)
         prd = 2.D0*cosh(theta(1))
         do i = 2, Nh
            prd = prd*2.D0*cosh(theta(i))
         enddo
         distribution(k)=prd*exp(sum(a(:)*pspin(:)))
         ii = ii + 1
      enddo
 
      distsum = sum(distribution(:))
      do k = 1, twopwr(Nv)
         write(200+itg,*) k, distribution(k)/distsum
      enddo

      vda(:) = vda(:)/float(Nw)
      vdb(:) = vdb(:)/float(Nw)
      vdwght(:,:) = vdwght(:,:)/float(Nw)

      da(:) = da(:)/float(Nw)
      db(:) = db(:)/float(Nw)
      dwght(:,:) = dwght(:,:)/float(Nw)
 
      egs = egs/float(Nw)

      da(:) = da(:) - egs*vda(:)
      db(:) = db(:) - egs*vdb(:)
      dwght(:,:) = dwght(:,:) - egs*vdwght(:,:)

      write(6,*) "Egs = ", egs
      stop
 
      a(:) = a(:) - eta*da(:)
      b(:) = b(:) - eta*db(:)
      wght(:,:) = wght(:,:) - eta*dwght(:,:)

      if(itg.gt.20) goto 30 
      itg = itg + 1
      goto 20

 30   continue
!     write(6,*) 
!     do k = 1, Nv
!        write(6,'(i5,x,f10.5)') k, da(k)
!     enddo

!     write(6,*) 
!     do k = 1, Nh
!        write(6,'(i5,x,f10.5)') k, db(k)
!     enddo

!     write(6,*) 
!     do i = 1, Nv
!        do k = 1, Nh
!           write(6,'(i5,x,i5,x,f10.5)') k,i, dwght(k,i)
!        enddo
!     enddo

! 

!     deallocate(spin,hidn,wght)
!     deallocate(a,b)
!     deallocate(gmm,theta)
!     deallocate(twopwr)
!     deallocate(icount)

      end program ising_vmc

      
