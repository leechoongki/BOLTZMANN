      module lgs

      double precision dxx
      double precision dcs
      double precision dsn
      double precision pi

      integer          nexp,ncos
      double precision xmax,smax
      parameter(xmax = 20.D0, smax = 20.D0)
      double precision, allocatable:: exp_save(:)
      double precision, allocatable:: cos_save(:)
      double precision, allocatable:: sin_save(:)
 
      end module lgs 

      module parameters

      integer Nv  !# of sites in visible layer
      integer Nh  !# of sites in hidden  layer
      integer Nw  !# of mcstep for warm up moving
      integer L  !  steps for measuring
      integer Np ! = Nv + Nh + Nv*Nh
      integer kshift ! = Nv + Nh
      integer npi ! # of gradient-Descent step
      double precision regul_factor
      double precision,allocatable::   prmts(:) ! network parameters order of Nv + Nh + Nv*Nh
      double complex,allocatable::  dprmts(:) ! \partial_k < Eloc >
      double complex,allocatable::  vprmts(:) ! < Eloc*D_k > - < Eloc > < D_k >

      end module parameters

      module spin_hidn
      
      logical(1), allocatable:: spin(:), hidn(:)
      double   complex, allocatable:: theta(:)
!     double precision, allocatable:: gmm_imag(:),gmm_real(:)
      double   complex, allocatable:: Dpsi(:,:)
      double precision hfield

      end module spin_hidn
   
