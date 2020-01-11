#!/usr/bin/gnuplot
set term png

set output "wykresy/mapa_zad1a.png"
set xlabel "x"
set ylabel "y"
set title "zad1a nx=ny=50"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:5][0:5] "mapa_zad1a.txt" i 0 u 2:1:3

reset

set output "wykresy/mapa_zad1b.png"
set xlabel "x"
set ylabel "y"
set title "zad1b nx=ny=100"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:10][0:10] "mapa_zad1b.txt" i 0 u 2:1:3

reset

set output "wykresy/mapa_zad1c.png"
set xlabel "x"
set ylabel "y"
set title "zad1c nx=ny=200"
set pm3d map
set palette defined (-10 "blue", 0 "white", 10 "red")
set size ratio -1

splot [0:20][0:20] "mapa_zad1c.txt" i 0 u 2:1:3

reset

set output "wykresy/mapa_zad2a.png"
set xlabel "x"
set ylabel "y"
set title "zad2a epsilon1=epsilon2= 1"
set pm3d map
set palette defined (-1 "blue", 0 "white", 1 "red")
set size ratio -1

splot [0:10][0:10][-0.8:0.8] "mapa_zad2a.txt" i 0 u 2:1:3


reset

set output "wykresy/mapa_zad2b.png"
set xlabel "x"
set ylabel "y"
set title "zad2b epsilon1=1 epsilon2= 2"
set pm3d map
set palette defined (-1 "blue", 0 "white", 1 "red")
set size ratio -1
set cbrange [-0.8:0.8]

splot [0:10][0:10][-0.8:0.8] "mapa_zad2b.txt" i 0 u 2:1:3


reset

set output "wykresy/mapa_zad2c.png"
set xlabel "x"
set ylabel "y"
set title "zad2c epsilon1=1 epsilon2= 10"
set pm3d map
set palette defined (-1 "blue", 0 "white", 1 "red")
set size ratio -1
set cbrange [-0.8:0.8]

splot [0:10][0:10][-0.8:0.8] "mapa_zad2c.txt" i 0 u 2:1:3
