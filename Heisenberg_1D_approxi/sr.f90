        subroutine Sto_Reconf(x)

        use parameters
        use spin_hidn
        use zminresqlpModule

        implicit none  
       
        integer i,j,k,ii
        integer itn
        integer nout
        double complex x(Np)
        double precision shift
        external:: Aprod
        double precision lambda0, bfactor, lambda_min
        parameter(lambda0=100.D0, bfactor=0.9D0, lambda_min=1.D-4)

        regul_factor = max(lambda0*bfactor**(npi),lambda_min) 

        shift = 0.D0
        nout = 999

        open(unit=999,file='minresql.dat',form="formatted")
        call zminresqlp(n=Np,Aprod=aprod,b=dprmts,shift=shift,x=x,&
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
        integer nseed
        double complex x(n), y(n)
        integer i,j,k 
        double complex sumi, vxsum

!
!       S_kl = < O_k^* O_l > - < O_k >^* < O_l > + regul_factor*S_kk*delta_kl 
!

!       sum_l S_kl x_l = < conjg(O_k)*(sum_l O_l*x_l) > - conjg(< O_k >)*( sum_l < O_l >*x_l ) 
!                      + regul_factor*( < |O_k|^2 > - |< O_k >|^2 )*x_k

        vxsum = sum(vprmts(:)*x(:))

        y(:) = 0.D0
        do k = 1, Nwr
           sumi = sum(Dpsi(:,k)*x(:))
           do i = 1, n
              y(i) = y(i) + dconjg(Dpsi(i,k))*sumi    &
                   + regul_factor*(dconjg(Dpsi(i,k))  &
                   * Dpsi(i,k))*x(i)
           enddo
        enddo 

!       y(:) = x(:)
        y(:) = y(:)/float(Nw) 
        y(:) = y(:)-dconjg(vprmts(:))*vxsum   &
             - regul_factor*dconjg(vprmts(:))*vprmts(:)*x(:)
!       do k = 1, n
!          write(6,*) k, y(k) 
!       enddo; stop

        return
        end
