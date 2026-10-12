set terminal pngcairo size 500,400 enhanced font 'Verdana,12'

# color definitions
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1 # --- blue
set style line 2 lc rgb '#ad6000' lt 1 lw 2 pt 5 ps 1 # --- red
set style line 6 lc rgb '#60ad00' lt 1 lw 2 pt 7 ps 1 # --- green
set style line 3 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 4 lc rgb '#ad6000' lt 1 lw 2 pt 7 ps 2 # --- red
set style line 5 lc rgb '#000000' lt 1 lw 2 pt 7 ps 2 # --- black

set lmargin 8
set bmargin 4
set xlabel "primary outbreak (%)"
set ylabel "mean secondary outbreak (%)"
set style fill  transparent solid 0.2 noborder
set style circle radius 1
unset key

#stats 'outfiles/secondary_complete' using 1 nooutput
#set xr [0:STATS_max*1.1] 
#stats 'outfiles/secondary_complete' using 2 nooutput
#set yr [0:STATS_max*1.1]
set xr [0:80]
set yr [0:20]
set title "complete symmetric network"
set output 'plots/secondary_complete.png'
plot 'outfiles/secondary_complete' using ($1):($2) with circles lc rgb "blue" t "each grid", 

set xr [0:80]
set yr [0:20]
set title "gamma network"
set output 'plots/secondary_gamma.png'
plot 'outfiles/secondary_gamma' using ($1):($2) with circles lc rgb "blue" t "each grid", 


set xr [0:80]
set yr [0:20]
set title "small-world ring"
set output 'plots/secondary_ring.png'
plot 'outfiles/secondary_ring' using ($1):($2) with circles lc rgb "blue" t "each grid",

set xr [0:80]
set yr [0:20]
set title "small-world grid"
set output 'plots/secondary_grid.png'
plot 'outfiles/secondary_grid' using ($1):($2) with circles lc rgb "blue" t "each grid", 

set xr [0:80]
set yr [0:20]
set title "random network"
set output 'plots/secondary_random.png'
plot 'outfiles/secondary_random' using ($1):($2) with circles lc rgb "blue" t "each grid", 

set xr [0:80]
set yr [0:20]
set title "ring-preferential network"
set output 'plots/secondary_ringpref.png'
plot 'outfiles/secondary_ringpref' using ($1):($2) with circles lc rgb "blue" t "each grid", 
