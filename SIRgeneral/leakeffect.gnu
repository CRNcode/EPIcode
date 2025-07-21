set terminal pngcairo size 1200,600 enhanced font 'Verdana,14'
#set terminal epslatex color 10 size 3.5in,2.4in
#set terminal eps font 'Times,18'


# color definitions
set style line 1 lc rgb '#0060ad' lt 1 lw 3 pt 5 ps 1.5 # --- blue
set style line 2 lc rgb '#ad6000' lt 1 lw 3 pt 4 ps 1.5 # --- red
set style line 3 lc rgb '#60ad00' lt 1 lw 3 pt 7 ps 1.5 # --- green
set style line 4 lc rgb '#000000' lt 1 lw 3 pt 2 ps 1.5 # --- black


set xrange [0:10]
set yrange [0:100]
#set ytics ("max" 1.0)
set lmargin 8
set bmargin 4
set xlabel "leak ϴ"
set ylabel "mean outbreak size (% of homogeneous)"
set key left top

#2D grid 5 X 3 Nc=50, rewiring = 10%
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_10.0_50_4_5_3_0.1 10000 2.0 0.0 10.0 50 4 5 3 0.1
set title "Expected outbreak size as a function of the leak"
set output 'plots/leak_2.0_0.0_10.0_50_4_5_3_0.1.png'
plot 'outfiles/leak_2.0_0.0_10.0_50_4_5_3_0.1' using ($1):(100*$2) w p ls 1 t "simulation",

#regular ring lattice n=15, Nc=50, rewiring = 10%
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_10.0_50_3_15_4_0.1 10000 2.0 0.0 10.0 50 3 15 4 0.1
set title "Expected outbreak size as a function of the leak"
set output 'plots/leak_2.0_0.0_10.0_50_3_15_4_0.1.png'
plot 'outfiles/leak_2.0_0.0_10.0_50_3_15_4_0.1' using ($1):(100*$2) w p ls 1 t "simulation",

#small "scale free" network; n=15, Nc=50, initial clique= 2
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_10.0_50_5_15_2 10000 2.0 0.0 10.0 50 5 15 2
set title "Expected outbreak size as a function of the leak"
set output 'plots/leak_2.0_0.0_10.0_50_5_15_2.png'
plot 'outfiles/leak_2.0_0.0_10.0_50_5_15_2' using ($1):(100*$2) w p ls 1 t "simulation",

set terminal pngcairo size 500,400 enhanced font 'Verdana,14'
set xrange [0:1]
set xlabel "ϴ/ϴ_{max}"
set ylabel "mean outbreak size (normalised)"
#complete, symmetric with 100 compartments (N_c=10)
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_9.9_10_1_100 10000 2.0 0.0 9.9 10 1 100
set title "complete, symmetric network"
set output 'plots/leak_2.0_0.0_9.9_10_1_100.png'
plot 'outfiles/leak_2.0_0.0_9.9_10_1_100' using ($1):(100*$2) w p ls 1 notitle

#ring (no rewiring, N_c=20)
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_16.0_20_3_100_4_0.0 1000 02.0 0.0 16.0 20 3 100 4 0.0
set title "rings"
set output 'plots/leak_2.0_0.0_16.0_20_3_100_4_0.0.png'
plot 'outfiles/leak_2.0_0.0_16.0_20_3_100_4_0.0' using ($1):(100*$2) w p ls 1 notitle

#ring (no rewiring, N_c=10,20,40)
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_8.0_10_3_100_4_0.0 10000 2.0 0.0 8.0 10 3 100 4 0.0
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_16.0_20_3_100_4_0.0 10000 2.0 0.0 16.0 20 3 100 4 0.0
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_32.0_40_3_100_4_0.0 10000 2.0 0.0 32.0 40 3 100 4 0.0
set title "rings"
set output 'plots/leak_2.0_0.0_X_X_3_100_4_0.0.png'
plot 'outfiles/leak_2.0_0.0_8.0_10_3_100_4_0.0' using ($1):(100*$2) w lp ls 1 t "N_c=10", 'outfiles/leak_2.0_0.0_16.0_20_3_100_4_0.0' using ($1):(100*$2) w lp ls 2 t "N_c=20", 'outfiles/leak_2.0_0.0_32.0_40_3_100_4_0.0' using ($1):(100*$2) w lp ls 3 t "N_c=40",

#2D grid (no rewiring)
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_4.0_5_4_10_10_0.0 10000 2.0 0.0 4.0 5 4 10 10 0.0
set title "grids"
set output 'plots/leak_2.0_0.0_4.0_5_4_10_10_0.0.png'
plot 'outfiles/leak_2.0_0.0_4.0_5_4_10_10_0.0' using ($1):(100*$2) w p ls 1 notitle,

#2D grid (no rewiring, N_c=10,20,40)
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_8.0_10_4_10_10_0.0 10000 2.0 0.0 8.0 10 4 10 10 0.0
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_16.0_20_4_10_10_0.0 10000 2.0 0.0 16.0 20 4 10 10 0.0
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_32.0_40_4_10_10_0.0 10000 2.0 0.0 32.0 40 4 10 10 0.0
set key bottom right
set title "grids"
set output 'plots/leak_2.0_0.0_X_X_4_10_10_0.0.png'
plot 'outfiles/leak_2.0_0.0_8.0_10_4_10_10_0.0' using ($1):(100*$2) w lp ls 1 t "N_c=10", 'outfiles/leak_2.0_0.0_16.0_20_4_10_10_0.0' using ($1):(100*$2) w lp ls 2 t "N_c=20", 'outfiles/leak_2.0_0.0_32.0_40_4_10_10_0.0' using ($1):(100*$2) w lp ls 3 t "N_c=40",

#ring preferential with larger compartments
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_66.67_100_7_100_10_2_1 10000 2.0 0.0 50 200 7 100 10 2 1
#set title "ring preferential"
set output 'plots/leak_2.0_0.0_66.67_100_7_100_10_2_1.png'
plot 'outfiles/leak_2.0_0.0_66.67_100_7_100_10_2_1' using ($1):(100*$2) w p ls 1 notitle,

#ring preferential (N_c=10,20,40)
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_6.67_10_7_100_10_2_1 10000 2.0 0.0 6.67 10 7 100 10 2 1
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_13.33_20_7_100_10_2_1 10000 2.0 0.0 13.33 20 7 100 10 2 1
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_26.67_40_7_100_10_2_1 10000 2.0 0.0 26.67 40 7 100 10 2 1
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_66.67_100_7_100_10_2_1 10000 2.0 0.0 50 200 7 100 10 2 1
set key bottom right
set title "ring-preferential networks"
set output 'plots/leak_2.0_0.0_X_X_7_100_10_2_1.png'
plot 'outfiles/leak_2.0_0.0_6.67_10_7_100_10_2_1' using ($1):(100*$2) w lp ls 1 t "N_c=10", 'outfiles/leak_2.0_0.0_13.33_20_7_100_10_2_1' using ($1):(100*$2) w lp ls 2 t "N_c=20", 'outfiles/leak_2.0_0.0_26.67_40_7_100_10_2_1' using ($1):(100*$2) w lp ls 3 t "N_c=40",
#'outfiles/leak_2.0_0.0_66.67_100_7_100_10_2_1' using ($1):(100*$2) w lp ls 3 t "N_c=100",

#random
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_4.0_5_6_100_4 10000 2.0 0.0 4.0 5 6 100 4
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_8.0_10_6_100_4 10000 2.0 0.0 8.0 10 6 100 4
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_16.0_20_6_100_4 10000 2.0 0.0 16.0 20 6 100 4
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_32.0_40_6_100_4 10000 2.0 0.0 32.0 40 6 100 4
set key bottom right
set title "random networks"
set output 'plots/leak_2.0_0.0_X_X_6_100_4.png'
plot 'outfiles/leak_2.0_0.0_4.0_5_6_100_4' using ($1):(100*$2) w lp ls 1 t "N_c=5", 'outfiles/leak_2.0_0.0_8.0_10_6_100_4' using ($1):(100*$2) w lp ls 2 t "N_c=10", 'outfiles/leak_2.0_0.0_16.0_20_6_100_4' using ($1):(100*$2) w lp ls 3 t "N_c=20",
#'SIRgeneral/outfiles/leak_2.0_0.0_32.0_40_6_100_4' using ($1):(100*$2) w lp ls 3 t "N_c=40",

#random network with mean degree 4 (N_c=10)
#./bin/SIRleakeffect SIRgeneral/outfiles/leak_2.0_0.0_8.0_10_6_100_4 10000 2.0 0.0 8.0 10 6 100 4
set title "random network"
set output 'plots/leak_2.0_0.0_8.0_10_6_100_4.png'
plot 'outfiles/leak_2.0_0.0_8.0_10_6_100_4' using ($1):(100*$2) w p ls 1 notitle,

#random network with mean degree 4 (N_c=10, constout=1)
#./bin/SIRleakeffect -c SIRgeneral/outfiles/leak_2.0_0.0_8.0_10_6_100_4 10000 2.0 0.0 8.0 10 6 100 4
#set title "random network"
#set output 'plots/leak_2.0_0.0_8.0_10_6_100_4c1.png'
#plot 'outfiles/leak_2.0_0.0_8.0_10_6_100_4c1' using ($1):(100*$2) w p ls 1 notitle,

set terminal pngcairo size 800,600 enhanced font 'Verdana,18'
set key top left
#put them together
set title "various architectures, and small compartments"
set output 'plots/leak_multi.png'
plot 'outfiles/leak_2.0_0.0_9.9_10_1_100' using ($1):(100*$2) w lp ls 4 t "1", 'outfiles/leak_2.0_0.0_8.0_10_6_100_4' using ($1):(100*$2) w lp ls 1 t "2", 'outfiles/leak_2.0_0.0_8.0_10_4_10_10_0.0' using ($1):(100*$2) w lp ls 2 t "3", 'outfiles/leak_2.0_0.0_8.0_10_3_100_4_0.0' using ($1):(100*$2) w lp ls 3 t "4", 
#1: complete, symmetric network
#2: random network
#3: 2D lattice
#4: regular ring lattice

