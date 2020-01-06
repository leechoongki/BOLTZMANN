        subroutine Sto_Reconf(x)

        use parameters
        use spin_hidn
        use minresqlpModule

        implicit none  
       
        integer i,j,k,ii
        integer itn
        integer nout
        double precision x(2*Np)
        double precision shift
        external:: Aprod

        regul_factor = 0.D0

        shift = 0.D0
        nout = 999
       
        open(unit=999,file='minresql.dat',form="formatted")
        call minresqlp(n=2*Np,Aprod=aprod,b=real(dprmts),shift=shift,x=x,&
                      nout=nout,itn=itn)
        close(999)
                     
        npi = npi + 1
        return 
        end subroutine Sto_Reconf

        subroutine aprod(n,x,y)

        use parameters
        use spin_hidn
        implicit none

        integer n
        double precision x(n), y(n)
        integer i,j,k 
        double complex sumi, vxsum
!
!       S_kl = < O_k^* O_l > - < O_k >^* < O_l >
!
!       sum_l S_kl x_l = < conjg(O_k)*(sum_l O_l*x_l) > - conjg(< O_k >)*( sum_l < O_l >*x_l ) 

        vxsum = sum(vprmts(:)*x(:))

        y(:) = 0.D0
        do k = 1, Nw
           sumi = sum(Dpsi(:,k)*x(:))
           do i = 1, n
              y(i) = y(i) + real(dconjg(Dpsi(i,k))*sumi) &  
                   + regul_factor*real(dconjg(Dpsi(i,k)) &
                   * Dpsi(i,k))*x(i) 
           enddo
        enddo 

        y(:) = y(:)/dfloat(Nw)-real(dconjg(vprmts(:))*vxsum)       &
             - regul_factor*real(dconjg(vprmts(:))*vprmts(:)*x(:))

        return
        end
