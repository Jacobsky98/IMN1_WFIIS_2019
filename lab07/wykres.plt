#!/usr/bin/gnuplot
set term png




set xl "x"
set yl "y"
set view map



set out 'mapaPozioma-1000.png'
set title 'Mapa pozioma Q=-1000'
splot 'wyniki/mapaPozioma-1000.txt' u 1:2:3 w pm3d t ""

set out 'mapaPionowa-1000.png'
set title 'Mapa pionowa Q=-1000'
splot 'wyniki/mapaPionowa-1000.txt' u 1:2:3 w pm3d t ""

set out 'mapaPozioma-4000.png'
set title 'Mapa pozioma Q=-4000'
splot 'wyniki/mapaPozioma-4000.txt' u 1:2:3 w pm3d t ""

set out 'mapaPionowa-4000.png'
set title 'Mapa pionowa Q=-4000'
splot 'wyniki/mapaPionowa-4000.txt' u 1:2:3 w pm3d t ""

set out 'mapaPozioma+4000.png'
set title 'Mapa pozioma Q=4000'
splot 'wyniki/mapaPozioma+4000.txt' u 1:2:3 w pm3d t ""

set out 'mapaPionowa+4000.png'
set title 'Mapa pionowa Q=4000'
splot 'wyniki/mapaPionowa+4000.txt' u 1:2:3 w pm3d t ""


unset key
set contour 
set cntrparam levels 50

set out 'psi_1000.png'
splot 'wyniki/wykresKonturowyPsi-1000.txt' u 1:2:3 w pm3d

set out 'psi_4000.png'
splot 'wyniki/wykresKonturowyPsi+4000.txt' u 1:2:3 w pm3d

set out 'psi_-4000.png'
splot 'wyniki/wykresKonturowyPsi-4000.txt' u 1:2:3 w pm3d

set out 'ksi_1000.png'
splot 'wyniki/wykresKonturowyKsi-1000.txt' u 1:2:3 w pm3d

set out 'ksi_4000.png'
splot 'wyniki/wykresKonturowyKsi+4000.txt' u 1:2:3 w pm3d

set out 'ksi_-4000.png'
splot 'wyniki/wykresKonturowyKsi-4000.txt' u 1:2:3 w pm3d

