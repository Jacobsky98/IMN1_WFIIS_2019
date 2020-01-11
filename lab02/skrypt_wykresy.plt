#!/usr/bin/gnuplot

set term png size 800,600
set encoding utf8

set output "wykres_zad_1.png"
set xlabel "t"
set ylabel "u(t), z(t)"
set title "Metoda niejawna trapezów - Picard"
plot "wynik_zad_1.txt" u 1:2 w l lw 2 t "u(t)", "" u 1:3 w l lw 2 t "z(t)"

set output "wykres_zad_2.png"
set xlabel "t"
set ylabel "u(t), z(t)"
set title "Metoda niejawna trapezów - Newton"
plot "wynik_zad_2.txt" u 1:2 w l lw 2 t "u(t)", "" u 1:3 w l lw 2 t "z(t)"



set output "wykres_zad_3.png"
set xlabel "t"
set ylabel "u(t), z(t)"
set title "Metoda niejawna RK2"
plot "wynik_zad_3.txt" u 1:2 w l lw 2 t "u(t)", "" u 1:3 w l lw 2 t "z(t)"
