set yr [0:1.0]
set xtics 1.0 
set ytics 0.2
set size 0.7,1.0
p 'h_2.0-1.0/out.dat' u 1:($2/5) title "4 -> 2" w p, 'h_0.25-0.5/out.dat' u 1:($2/5) title "1/2 -> 1" w p
pause -1

