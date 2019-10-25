set term postscript enhanced color
set output "a.eps"
set xr [1:1024]
set yr [:0.3]
set size 0.7, 0.7
set pointsize 0.5
set format y "%.2f"
plot "fort.7" u 1:2 ti "Gibbs sampling" w l,\
     "fort.8" u 1:2 ti "|{/Symbol Y}|^2" w p pt 7
pause -1
