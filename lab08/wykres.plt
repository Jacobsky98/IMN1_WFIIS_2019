#!/usr/bin/gnuplot
set term png

set output "wykresy/calka_xsr.png"
set xlabel "tn"
set ylabel "x_s_r, c(tn)"
set title "Calka c(tn) D=0.0"
plot "wyniki/calka_0.0.txt" u 1:2 w l t "c(tn) D=0.0", \
        "wyniki/sredniaPolozenia_0.0.txt" u 1:2 w l t "x_s_r(tn) D=0.0", \
        "wyniki/calka_0.1.txt" u 1:2 w l t "c(tn) D=0.1", \
        "wyniki/sredniaPolozenia_0.1.txt" u 1:2 w l t "x_s_r(tn) D=0.1"



set xl "x"
set yl "y"
set view map


set out 'wykresy/mapa_vx_0.0.png'
set title 'Mapa vx D=0.0'
splot 'wyniki/mapa_vx1_0.0.txt' u 1:2:3 w pm3d t ""

set out 'wykresy/mapa_vy_0.0.png'
set title 'Mapa vy D=0.0'
splot 'wyniki/mapa_vy1_0.0.txt' u 1:2:3 w pm3d t ""

set out 'wykresy/mapa_vx_0.1.png'
set title 'Mapa vx D=0.1'
splot 'wyniki/mapa_vx1_0.1.txt' u 1:2:3 w pm3d t ""

set out 'wykresy/mapa_vy_0.1.png'
set title 'Mapa vy D=0.1'
splot 'wyniki/mapa_vy1_0.1.txt' u 1:2:3 w pm3d t ""







