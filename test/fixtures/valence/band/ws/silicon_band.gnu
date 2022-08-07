set style data dots
set nokey
set xrange [0: 3.81842]
set yrange [ -6.72649 :  7.40924]
set arrow from  1.00811,  -6.72649 to  1.00811,   7.40924 nohead
set arrow from  2.17218,  -6.72649 to  2.17218,   7.40924 nohead
set arrow from  2.58374,  -6.72649 to  2.58374,   7.40924 nohead
set xtics ("L"  0.00000,"G"  1.00811,"X"  2.17218,"K"  2.58374,"G"  3.81842)
 plot "silicon_band.dat"
