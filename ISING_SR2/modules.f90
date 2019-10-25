      module lgs

      double precision dxx
      integer          nlgs
      double precision xmax
      double precision, allocatable:: logisave(:)
 
      end module lgs 

      module parameters

      integer Nv !# of sites in visible layer
      integer Nh !# of sites in hidden  layer
      integer Nw !# of mcstep for warm up moving
      integer L  !  steps for measuring
      integer Np
      integer kshift
      integer npi ! # of gradient-Descent step
      double precision regul_factor
      double precision,allocatable::   prmts(:) ! network parameters order of Nv + Nh + Nv*Nh
      double precision,allocatable::  dprmts(:) ! network parameters order of Nv + Nh + Nv*Nh
      double precision,allocatable::  vprmts(:)

      end module parameters

      module spin_hidn
      
      integer, allocatable:: spin(:), hidn(:)
      double precision, allocatable:: gmm(:), theta(:), Dpsi(:,:)
      double precision Jex, hfield

      end module spin_hidn
   
