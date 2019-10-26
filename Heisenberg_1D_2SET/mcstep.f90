      subroutine mcstep()

      use parameters
      use spin_hidn

      implicit none

      integer j,k
      real rs
      integer nselect,isig
      double precision tmp(Nhn)
      double precision aratio
      double precision prd, w
      double precision, external:: coshl, expl
 
      call random_number(rs)
      nselect = int(rs*Nv)+1
      isig = merge(-2,2,spin(nselect))

      do j = 1, Nhn 
         tmp(j)=theta_norm(j)+isig*prmts(kshift_norm+nselect+(j-1)*Nv)
      enddo
      aratio=exp(isig*prmts(nselect))

      prd = 1.D0
      do j = 1, Nhn
         prd = prd*cosh(tmp(j))/cosh(theta_norm(j))
      enddo

      w = aratio * prd 

      if(w.ge.1) then 
        theta_norm(:) = tmp(:)
        spin(nselect) = .not.spin(nselect)
        do j = 1, Nhp
           theta_phas(j)=theta_phas(j) &
                        +isig*prmts(kshift_phas+nselect+(j-1)*Nv)
        enddo
      else
        call random_number(rs)
        if(w.gt.rs) then
           theta_norm(:) = tmp(:)
           spin(nselect) = .not.spin(nselect)
           do j = 1,Nhp
              theta_phas(j)=theta_phas(j) &
                           +isig*prmts(kshift_phas+nselect+(j-1)*Nv)
           enddo
           return
        else 
           return
        endif
      endif

      end subroutine mcstep
