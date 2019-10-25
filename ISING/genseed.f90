!      subroutine generate_seed(iseed)
!      integer time, stime, t(9), iseed   
!      stime = time(%ref(0))
!      call gmtime(stime,t)
!      iseed = t(6) + 70*(t(5)+12*(t(4)+31*(t(3)+23*(t(2)+59*t(1)))))
!      if(mod(iseed,2).eq.0) iseed = iseed - 1
!      return
!      end

        SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
        END SUBROUTINE
