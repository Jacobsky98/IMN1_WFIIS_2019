#!/usr/bin/gnuplot
set term png

set output "calka_s_wszystkie.png"
set xlabel "i"
set ylabel "S(i)"
set title "Calka S(i)"
plot "wyniki/calka_16.txt" u 1:2 w l t "k=16", \
        "wyniki/calka_8.txt" u 1:2 w l t "k=8", \
        "wyniki/calka_4.txt" u 1:2 w l t "k=4", \
        "wyniki/calka_2.txt" u 1:2 w l t "k=2", \
        "wyniki/calka_1.txt" u 1:2 w l t "k=1"


set output "calka_16.png"
set xlabel "it"
set ylabel "S(it)"
set title "Calka S(i) k=16"
plot "wyniki/calka_16.txt" u 1:2 w l t ""



set output "calka_8.png"
set xlabel "it"
set ylabel "S(it)"
set title "Calka S(i) k=8"
plot "wyniki/calka_8.txt" u 1:2 w l t ""



set output "calka_4.png"
set xlabel "it"
set ylabel "S(it)"
set title "Calka S(i) k=4"
plot "wyniki/calka_4.txt" u 1:2 w l t ""



set output "calka_2.png"
set xlabel "it"
set ylabel "S(it)"
set title "Calka S(i) k=2"
plot "wyniki/calka_2.txt" u 1:2 w l t ""


set output "calka_1.png"
set xlabel "it"
set ylabel "S(it)"
set title "Calka S(i) k=1"
plot "wyniki/calka_1.txt" u 1:2 w l t ""


unset grid
set view map
set xrange [0:25.6]
set yrange [0:25.6]
set xlabel "x"
set ylabel "y"
set size square

set out 'mapa_16.png'
set title 'Mapa k=16'
splot 'wyniki/mapa_16.txt' u 1:2:3 w pm3d t ""

set out 'mapa_8.png'
set title 'Mapa k=8'
splot 'wyniki/mapa_8.txt' u 1:2:3 w pm3d t ""

set out 'mapa_4.png'
set title 'Mapa k=4'
splot 'wyniki/mapa_4.txt' u 1:2:3 w pm3d t ""

set out 'mapa_2.png'
set title 'Mapa k=2'
splot 'wyniki/mapa_2.txt' u 1:2:3 w pm3d t ""

set out 'mapa_1.png'
set title 'Mapa k=1'
splot 'wyniki/mapa_1.txt' u 1:2:3 w pm3d t ""
