#Individual trajectories, deterministic mean, stochastic mean
#set terminal pngcairo size 1000,500 enhanced font 'Verdana,14'
set terminal epslatex standalone size 17cm,8.5cm 
set output "Epi5_all_Rt.tex"
#data created with ./bin/SIRgeneral 1 2.0 2.0 1000 1 5

stats '../outfiles/Epi5_all_full_outbreakA' using 2 nooutput
#set title "Tracking measures of population protection during an outbreak"
set ytics 0.5 nomirror
set ylabel '$\widehat{\mathcal{R}}_t$, $\mathcal{R}_t$ and mean outbreak size'
set y2tics 100 nomirror
set y2label 'infections'
set xlabel "time"
#set ylabel "R_t"
#set key top left
plot [0:50] 1 w l notitle, "../outfiles/Epi5_all_full" using 1:($8) w l lc "purple" lw 2 t "subpopulations 1-5" axes x1y2, "../outfiles/Epi5_all_full" using 1:($11) w l lc "purple" lw 2 notitle axes x1y2, "../outfiles/Epi5_all_full" using 1:($14) w l lc "purple" lw 2 notitle axes x1y2, "../outfiles/Epi5_all_full" using 1:($17) w l lc "purple" lw 2 notitle axes x1y2, "../outfiles/Epi5_all_full" using 1:($20) w l lc "purple" lw 2 notitle axes x1y2, "../outfiles/Epi5_all_full" using 1:($5) w l lc "blue" lw 4 t "total infections" axes x1y2, "../outfiles/Epi5_all_full" using 1:($2) w l lc "black" lw 4 t "$\\widehat{\\mathcal{R}}_t$", "../outfiles/Epi5_all_full" using 1:($3) w l lc "red" lw 4 t "$\\mathcal{R}_t$", "../outfiles/Epi5_all_full_outbreakA" using 1:($2 * 2./STATS_max) w lp lc "black" pt 7 ps 1.0 lw 2 t "mean outbreak size"

set output
