set terminal pngcairo size 1000,500 enhanced font 'Verdana,14'


# color definitions
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1 # --- blue
set style line 2 lc rgb '#ad6000' lt 1 lw 2 pt 5 ps 1 # --- red
set style line 6 lc rgb '#60ad00' lt 1 lw 2 pt 7 ps 1 # --- green
set style line 3 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 4 lc rgb '#ad6000' lt 1 lw 2 pt 7 ps 2 # --- red
set style line 5 lc rgb '#000000' lt 1 lw 2 pt 7 ps 2 # --- black



set xrange [0:60]
#set yrange [0:20]
#set ytics ("max" 1.0)
set lmargin 10
set bmargin 4
set xlabel "prior immunity (%)"
set ylabel "mean outbreak size (%)"

#set label "infected" at -2,0.4 rotate left  font "Times,16"
#set label "time" at 22,-0.1 font "Times,16"


set title "Expected outbreak sizes as a function of prior immunity"
set output 'plots/random_10c.png'
plot 'outfiles/random_1_10000a' using ($1*100):($3*100) w l lc "black" lw 2 t "1 compartment theory", 'outfiles/random_1_10000a' using ($1*100):($5*100) w p pt 1 lc "blue" ps 1.5 t "1 compartment simulation", 'outfiles/random_2_5000a' using ($1*100):($3*100) w l lc "green" lw 2 t "2 compartment theory", 'outfiles/random_2_5000a' using ($1*100):($5*100) w p pt 5 lc "green" t "2 compartment simulation", 'outfiles/random_10_1000a' using ($1*100):($3*100) w l lc "red" lw 2 t "10 compartment theory", 'outfiles/random_10_1000a' using ($1*100):($5*100) w p pt 7 lc "red" t "10 compartment simulation"
