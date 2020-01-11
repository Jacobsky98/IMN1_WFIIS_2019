set term png

set out 'relaksacja_globalna_calka_s.png'
set logscale x
plot 'wyniki/relaksacja_globalna_0.6.txt' u 1:2 w l t'wG = 0.6',\
'wyniki/relaksacja_globalna_1.0.txt' u 1:2 w l t'wG = 1.0'

set out 'relaksacja_lokalna_calka_s.png'
set logscale x
plot 'wyniki/relaksacja_lokalna_1.0.txt' u 1:2 w l t'wL = 1.0',\
'wyniki/relaksacja_lokalna_1.4.txt' u 1:2 w l t'wL = 1.4',\
'wyniki/relaksacja_lokalna_1.8.txt' u 1:2 w l t'wL = 1.8',\
'wyniki/relaksacja_lokalna_1.9.txt' u 1:2 w l t'wL = 1.9'

set out 'mapa_globalna_0.6.png'
set title 'mapa zrelaksowanego potencjalu wg = 0.6'
plot 'wyniki/mapa_globalna_0.6.txt' matrix nonuniform with image

set out 'mapa_bledu_0.6.png'
set title 'mapa bledu potencjalu wg = 0.6'
plot 'wyniki/mapa_globalna_bledu_0.6.txt' matrix nonuniform with image


set out 'mapa_globalna_1.0.png'
set title 'mapa zrelaksowanego potencjalu wg = 1.0'
plot 'wyniki/mapa_globalna_1.0.txt' matrix nonuniform with image

set out 'mapa_bledu_1.0.png'
set title 'mapa bledu potencjalu wg = 1.0'
plot 'wyniki/mapa_globalna_bledu_1.0.txt' matrix nonuniform with image