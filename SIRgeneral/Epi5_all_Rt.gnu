#Individual trajectories, deterministic mean, stochastic mean
set terminal pngcairo size 1000,500 enhanced font 'Verdana,14'
set output "plots/Epi5_all_Rt.png"
#data created with ./bin/SIRgeneral 1 2.0 2.0 1000 1 5

stats 'outfiles/Epi5_all_full_outbreakA' using 2 nooutput
set title "Tracking measures of population protection during an outbreak"
set ytics 0.5 nomirror
set ylabel 'R_t / expected outbreak size'
set y2tics 100 nomirror
set y2label 'infections'
set xlabel "time"
#set ylabel "R_t"
#set key top left
plot [0:50] 1 w l notitle, "outfiles/Epi5_all_full" using 1:($8) w l lc "purple" lw 1 t "compartments 1-5" axes x1y2, "outfiles/Epi5_all_full" using 1:($11) w l lc "purple" lw 2 notitle axes x1y2, "outfiles/Epi5_all_full" using 1:($14) w l lc "purple" lw 2 notitle axes x1y2, "outfiles/Epi5_all_full" using 1:($17) w l lc "purple" lw 2 notitle axes x1y2, "outfiles/Epi5_all_full" using 1:($20) w l lc "purple" lw 2 notitle axes x1y2, "outfiles/Epi5_all_full" using 1:($5) w l lc "blue" lw 2 t "total infections" axes x1y2, "outfiles/Epi5_all_full" using 1:($2) w l lc "black" lw 2 t "classical R_t", "outfiles/Epi5_all_full" using 1:($3) w l lc "red" lw 2 t "next generation R_t", "outfiles/Epi5_all_full_outbreakA" using 1:($2 * 2./STATS_max) w lp lc "black" pt 7 ps 0.5 lw 1 t "mean outbreak size"
