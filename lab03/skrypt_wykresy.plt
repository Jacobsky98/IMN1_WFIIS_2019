#!/usr/bin/gnuplot

set term png size 800,600
set encoding utf8

    

    ### Metoda trapezow ###

set output "wykresy/wykres_metoda_trapezow_v(t).png"
set xlabel "t"
set ylabel "v"
set title "Metoda trapezow v(t)"
plot "metoda_trapezow_TOL-5.txt" u 1:4 w l lw 2 t "TOL=10^{-5}", "metoda_trapezow_TOL-2.txt" u 1:4 w l lw 2 t "TOL=10^{-2}"


set output "wykresy/wykres_metoda_trapezow_x(t).png"
set xlabel "t"
set ylabel "x"
set title "Metoda trapezow x(t)"
plot "metoda_trapezow_TOL-5.txt" u 1:3 w l lw 2 t "TOL=10^{-5}", "metoda_trapezow_TOL-2.txt" u 1:3 w l lw 2 t "TOL=10^{-2}"


set output "wykresy/wykres_metoda_trapezow_dt(t).png"
set xlabel "t"
set ylabel "dt"
set title "Metoda trapezow dt(t)"
plot "metoda_trapezow_TOL-5.txt" u 1:2 w l lw 2 t "TOL=10^{-5}", "metoda_trapezow_TOL-2.txt" u 1:2 w l lw 2 t "TOL=10^{-2}"


set output "wykresy/wykres_metoda_trapezow_v(x).png"
set xlabel "x"
set ylabel "v"
set title "Metoda trapezow v(x)"
plot "metoda_trapezow_TOL-5.txt" u 3:4 w l lw 2 t "TOL=10^{-5}", "metoda_trapezow_TOL-2.txt" u 3:4 w l lw 2 t "TOL=10^{-2}"


### Metoda RK2 ###

set output "wykresy/wykres_metoda_RK2_v(t).png"
set xlabel "t"
set ylabel "v"
set title "Metoda RK2 v(t)"
plot "metoda_RK2_TOL_-5.txt" u 1:4 w l lw 2 t "TOL=10^{-5}", "metoda_RK2_TOL_-2.txt" u 1:4 w l lw 2 t "TOL=10^{-2}"


set output "wykresy/wykres_metoda_RK2_x(t).png"
set xlabel "t"
set ylabel "x"
set title "Metoda RK2 x(t)"
plot "metoda_RK2_TOL_-5.txt" u 1:3 w l lw 2 t "TOL=10^{-5}", "metoda_RK2_TOL_-2.txt" u 1:3 w l lw 2 t "TOL=10^{-2}"


set output "wykresy/wykres_metoda_RK2_dt(t).png"
set xlabel "t"
set ylabel "dt"
set title "Metoda RK2 dt(t)"
plot "metoda_RK2_TOL_-5.txt" u 1:2 w l lw 2 t "TOL=10^{-5}", "metoda_RK2_TOL_-2.txt" u 1:2 w l lw 2 t "TOL=10^{-2}"


set output "wykresy/wykres_metoda_RK2_v(x).png"
set xlabel "x"
set ylabel "v"
set title "Metoda RK2 v(x)"
plot "metoda_RK2_TOL_-5.txt" u 3:4 w l lw 2 t "TOL=10^{-5}", "metoda_RK2_TOL_-2.txt" u 3:4 w l lw 2 t "TOL=10^{-2}"





