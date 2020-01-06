      subroutine mcstep()

      use parameters
      use spin_hidn

      implicit none

      integer j,k
      real rs
      integer nselect
      double   complex tmp(Nh)
      double precision aratio
      double precision prd, w
      double precision, external:: cosl, coshl, expl
 
      call random_number(rs)
      nselect = int(rs*Nv)+1

      if(spin(nselect)) then
         do j = 1, Nh 
            tmp(j) = theta(j)                                         &
                   - 2.D0*dcmplx(prmts(kshift+2*nselect-1+(j-1)*2*Nv),&
                                 prmts(kshift+2*nselect  +(j-1)*2*Nv))
         enddo
         aratio = expl(-4.D0*prmts(2*nselect-1))
      else
         do j = 1, Nh 
            tmp(j) = theta(j)                                          &
                   + 2.D0*dcmplx(prmts(kshift+2*nselect-1+(j-1)*2*Nv),&
                                 prmts(kshift+2*nselect  +(j-1)*2*Nv)) 
         enddo
         aratio = expl( 4.D0*prmts(2*nselect-1))
      endif

      prd = 1.D0
      do j = 1, Nh
         prd = prd*(coshl(2.D0*real(tmp(j)))+cosl(2.D0*aimag(tmp(j))))/&
                   (coshl(2.D0*real(theta(j)))+cosl(2.D0*aimag(theta(j)))) 
      enddo

      w = aratio * prd 

      if(w.ge.1) then 
        theta(:) = tmp(:)
        spin(nselect) = .not.spin(nselect)
        return
      else
        call random_number(rs)
        if(w.gt.rs) then
           theta(:) = tmp(:)
           spin(nselect) = .not.spin(nselect)
           return
        else 
           return
        endif
      endif

      return
      end
