       program test

       implicit none
   
       integer(4) i
       real rs, rs2
       double precision sumdab,pi

       pi = acos(-1.0)

       sumdab = 0.0
       call cpu_time(rs)
       do i = 1, 10000
          sumdab = sumdab + 1.0/(float(i*i))
          write(6,*) sumdab, pi*pi/6.0
       enddo
       call cpu_time(rs2)

       write(6,*) rs2-rs, sumdab

       end
