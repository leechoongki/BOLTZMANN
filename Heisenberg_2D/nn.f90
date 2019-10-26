       program nn_test

       implicit none
       integer nx,ny,nv,nz,k,i
       integer,allocatable:: innd(:,:)

        nx = 3
        ny = 3
        nz = 4
        nv = nx*ny
        allocate(innd(nv,nz))
        call find_nn(nv,nx,ny,innd)

        do k = 1, nv
           write(6,*)k, (innd(k,i),i=1,nz)
        enddo
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
