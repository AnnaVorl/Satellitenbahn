# Flugbahn

Nachdem das Programm gelaufen ist und die entsprechenden Werte berechnet hat und in die File 'Ergebnisse.dat' geschrieben hat, 
kann die Flugbahn mit Gnuplot und folgenden Befehlen geplottet werden: 

set title 'Flugbahn der ISS'


set grid


set xlabel 'x-Position in m'


set ylabel 'y-Position in m'


set xrange [-8e6:8e6]


set yrange [-8e6:8e6]


plot 'Ergebnisse.dat' using 2:3 w lp lc 6 lw 2 t 'Bahnkurve'

