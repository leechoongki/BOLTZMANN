      module lgs

      integer nexp
      double precision dxx
      double precision pi

      double precision xmax
      parameter(xmax = 20.D0)
      double precision, allocatable:: exp_save(:)
 
      end module lgs 

      module parameters

      integer Nv     !# of sites in visible layer
      integer Nhn    !# of sites in a hidden  layer for norm
      integer Nhp    !# of sites in a hidden  layer for phase
      integer ino    !# Nhn = ino * Nv
      integer ipa    !# Nhp = ipa * Nv
      integer Nw     !# of mcstep for warm up moving
      integer L      !  jumping steps for measuring
      integer Np     !# of wave function parameters = Nv + Nhn + Nhp + (Nhn+Nhp)*Nv
      integer kshift_norm ! = Nv + Nhn + Nhp
      integer kshift_phas ! = Nv + Nhn + Nhp + Nhn*Nv
      integer npi    !# of SR step
      double precision regul_factor
      double precision,allocatable::   prmts(:) ! network parameters order of Nv + Nhn + Nhp + Nhn*Nv + Nhp*Nv
      double   complex,allocatable::  dprmts(:) ! \partial_k < Eloc >
      double   complex,allocatable::  vprmts(:) ! < Eloc*D_k > - < Eloc > < D_k >

      end module parameters

      module spin_hidn
      
      logical, allocatable:: spin(:)
      double precision,allocatable:: theta_norm(:),theta_phas(:)
      double   complex,allocatable:: Dpsi(:,:)
      double precision Jz

      end module spin_hidn
   
