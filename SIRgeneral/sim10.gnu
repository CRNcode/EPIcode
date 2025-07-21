set terminal pngcairo size 800,600 enhanced font 'Verdana,18'
#set terminal epslatex color 10 size 3.5in,2.4in
#set terminal eps font 'Times,18'
set output "plots/sim_10.png"

# color definitions
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1 # --- blue
set style line 2 lc rgb '#ad6000' lt 1 lw 2 pt 7 ps 1 # --- red
set style line 6 lc rgb '#60ad00' lt 1 lw 2 pt 5 ps 1.5 # --- green
set style line 3 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5 # --- blue
set style line 4 lc rgb '#cc0000' lt 1 lw 4 pt 7 ps 2 # --- red
set style line 5 lc rgb '#000000' lt 1 lw 2 pt 2 ps 1.5 # --- black


set xrange [0:5]
set yrange [0:100]
#set ytics ("max" 1.0)
set lmargin 8
set bmargin 4
set xlabel "ϴ N_c"
set ylabel "mean outbreak size (normalised)"
set key bottom right

set title "complete, symmetric networks"
plot 'outfiles/sim_10_1000' using ($1):(100*$2) w l ls 4 t "theory (N_c = 1000)", \
	 'outfiles/sim_10_100' using ($1):(100*$3) w p ls 3 t "simulation (N_c = 100)", \
	 'outfiles/sim_10_200' using ($1):(100*$3) w p ls 6 t "simulation (N_c = 200)", \
         'outfiles/sim_10_1000' using ($1):(100*$3) w p ls 5 t "simulation (N_c = 1000)",
#'outfiles/sim_10_200' using ($1):(100*$2) w l ls 4 t "theoretical approximation", \
#'outfiles/sim_10_100' using ($1):(100*$2) w l ls 4 t "theoretical approximation", \

set xrange [0:0.12]
set yrange [0:100]
set title "Expected outbreak size, conditional on an outbreak occurring (50 compartments)"
set output 'plots/sim_50_50.png'
plot 'outfiles/sim_50' using ($1/49):(100*$2) w l ls 4 t "theoretical approximation", \
	 'outfiles/sim_50_50' using ($1/49):(100*$3) w p ls 3 t "Gillespie simulation (N_c = 50)", 


