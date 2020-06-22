#!/usr/bin/gnuplot
#unset key
set term png
set output "band.png"
plot './band.dat' u 1:2 w l lc 3, \
     './band.dat' u 1:3 w l lc 3
