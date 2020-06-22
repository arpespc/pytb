#!/usr/bin/gnuplot
#unset key
set term png
set output "band.png"
set yrange [-5:9]
plot './band.dat' u 1:2 w l lc 3, \
     './band.dat' u 1:3 w l lc 3, \
     './band.dat' u 1:4 w l lc 3
