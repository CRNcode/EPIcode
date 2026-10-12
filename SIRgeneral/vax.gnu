set terminal pngcairo size 600,400 enhanced font 'Verdana,14'


# color definitions
set style line 1 lc rgb '#0060ad' lt 1 lw 3 pt 5 ps 1 # --- blue
set style line 2 lc rgb '#ad6000' lt 1 lw 3 pt 4 ps 1 # --- red
set style line 3 lc rgb '#60ad00' lt 1 lw 3 pt 7 ps 1 # --- green
set style line 4 lc rgb '#000000' lt 1 lw 3 pt 2 ps 1 # --- black
set style line 4 lc rgb '#ad0000' lt 1 lw 3 pt 2 ps 1 # --- red



set xrange [0:60]
#set yrange [0:20]
#set ytics ("max" 1.0)
set lmargin 10
set bmargin 4
set xlabel "prior random immunity (%)"
set ylabel "mean outbreak size (%)"

#set label "infected" at -2,0.4 rotate left  font "Times,16"
#set label "time" at 22,-0.1 font "Times,16"


#./bin/SIRvax SIRgeneral/outfiles/vax_1_10000a 1 10000 2.0 0.0 10000 1 1
#./bin/SIRvax SIRgeneral/outfiles/vax_2_5000a 1 10000 2.0 4.0 5000 1 2
#./bin/SIRvax SIRgeneral/outfiles/vax_10_1000a 1 10000 2.0 6.0 1000 1 10
#set title "Outbreak sizes as a function of random immunity"
set title "complete, symmetric networks"
set output 'plots/vax_10c.png'
plot 'outfiles/vax_1_10000a' using ($1*100):($3*100) w l lc "black" lw 2 t "theory (n=1)", 'outfiles/vax_1_10000a' using ($1*100):($5*100) w p pt 1 lc "blue" ps 1.5 t "simulation (n=1)", 'outfiles/vax_2_5000a' using ($1*100):($3*100) w l lc "green" lw 2 t "theory (n=2)", 'outfiles/vax_2_5000a' using ($1*100):($5*100) w p pt 5 lc "green" t "simulation (n=2)", 'outfiles/vax_10_1000a' using ($1*100):($3*100) w l lc "red" lw 2 t "theory (n=10)", 'outfiles/vax_10_1000a' using ($1*100):($5*100) w p pt 7 lc "red" t "simulation (n=10)"

#./bin/SIRvax SIRgeneral/outfiles/vax_2.0_15.0_100_4_10_10_0.05 1 10000 2.0 15.0 100 4 10 10 0.05
#./bin/SIRvax SIRgeneral/outfiles/vax_2.0_25.0_100_3_100_4_0.1 1 10000 2.0 25.0 100 3 100 4 0.1
#./bin/SIRvax SIRgeneral/outfiles/vax_2.0_20.0_100_6_100_6 1 10000 2.0 20.0 100 6 100 6
#set terminal pngcairo size 1000,500 enhanced font 'Verdana,14'
#set title "Outbreak sizes as a function of random immunity"
set title "small-world networks"
set output 'plots/vax_multi.png'
plot 'outfiles/vax_2.0_15.0_100_4_10_10_0.0' using ($1*100):($4*100) w lp ls 1 t "small-world grid", 'outfiles/vax_2.0_25.0_100_3_100_4_0.1' using ($1*100):($4*100) w lp ls 4 t "small-world ring", 'outfiles/vax_2.0_15.0_100_4_10_10_0.05' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case",
#'outfiles/vax_2.0_20.0_100_6_100_6' using ($1*100):($4*100) w lp ls 3 t "random",

#Can see the effect clearly in a weakly coupled random network
#./bin/SIRvax SIRgeneral/outfiles/vax_2.0_20.0_1000_6_25_5 10 1000 2.0 20.0 1000 6 25 5
set title "random network"
set output 'plots/vax_2.0_20.0_1000_6_25_5.png'
plot 'outfiles/vax_2.0_20.0_1000_6_25_5' using ($1*100):($4*100) w lp ls 1 t "simulations", 'outfiles/vax_2.0_20.0_1000_6_25_5' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case",

#Can see the effect clearly in a weakly coupled ring-preferential network
#./bin/SIRvax SIRgeneral/outfiles/vax_2.0_20.0_2000_7_50_6_2_1 10 1000 2.0 20.0 2000 7 50 6 2 1
set title "ring-preferential network"
set output 'plots/vax_2.0_20.0_2000_7_50_6_2_1.png'
plot 'outfiles/vax_2.0_20.0_2000_7_50_6_2_1' using ($1*100):($4*100) w lp ls 1 t "simulations", 'outfiles/vax_2.0_20.0_2000_7_50_6_2_1' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case",

#Here (strongly coupled random network) we see little difference with the homogeneous case
#./bin/SIRvax SIRgeneral/outfiles/vax_2.0_667_1000_6_20_6 10 1000 2.0 667 1000 6 20 6
#cp outfiles/network_properties outfiles/network_properties_2.0_667_1000_6_20_6
set title "random network"
set output 'plots/vax_2.0_667_1000_6_20_6.png'
plot 'outfiles/vax_2.0_667_1000_6_20_6' using ($1*100):($4*100) w lp ls 1 t "random network", 'outfiles/vax_2.0_667_1000_6_20_6' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case",

#Here (strongly coupled ring preferential) random immunity is less effective than in a homogeneous network
#./bin/SIRvax SIRgeneral/outfiles/vax_2.0_667_1000_7_20_4_2_1 10 1000 2.0 667 1000 7 20 4 2 1
#cp outfiles/network_properties outfiles/network_properties_2.0_667_1000_7_20_4_2_1
set title "ring preferential network"
set output 'plots/vax_2.0_667_1000_7_20_4_2_1.png'
plot 'outfiles/vax_2.0_667_1000_7_20_4_2_1' using ($1*100):($4*100) w lp ls 1 t "simulations", 'outfiles/vax_2.0_667_1000_7_20_4_2_1' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case",



################


#set output 'plots/vax.png'
#plot 'outfiles/vax' using ($1*100):($4*100) w p pt 7 lc "blue" ps 1.5 t "simulation", 'outfiles/vax' using ($1*100):($2*100) w l lc "black" lw 2 t "1 compartment theory",

#./bin/SIRvax SIRgeneral/outfiles/vax_2.0_15.0_100_4_13_9_0.1 10 1000 2.0 15.0 100 4 13 9 0.1
set terminal pngcairo size 500,400 enhanced font 'Verdana,14'
set title "2D lattice based small world network"
set output 'plots/vax_2.0_15.0_100_4_13_9_0.1.png'
plot 'outfiles/vax_2.0_15.0_100_4_13_9_0.1' using ($1*100):($4*100) w p pt 1 lc "black" ps 1.5 t "simulation", 'outfiles/vax_2.0_15.0_100_4_13_9_0.1' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case", 

#./bin/SIRvax SIRgeneral/outfiles/vax_2.0_15.0_100_4_10_10_0.0 1 10000 2.0 15.0 100 4 10 10 0.0
set title "2D lattice"
set output 'plots/vax_2.0_20.0_100_4_10_10_0.0.png'
plot 'outfiles/vax_2.0_20.0_100_4_10_10_0.0' using ($1*100):($4*100) w p pt 1 lc "black" ps 1.5 t "simulation", 'outfiles/vax_2.0_20.0_100_4_10_10_0.0' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case",

#./bin/SIRvax SIRgeneral/outfiles/vax_2.0_25.0_100_3_100_4_0.0 1 10000 2.0 25.0 100 3 100 4 0.0
set title "regular ring network"
set output 'plots/vax_2.0_25.0_100_3_100_4_0.0.png'
plot 'outfiles/vax_2.0_25.0_100_3_100_4_0.0' using ($1*100):($4*100) w p pt 1 lc "black" ps 1.5 t "simulation", 'outfiles/vax_2.0_25.0_100_3_100_4_0.0' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case",



#./bin/SIRvax SIRgeneral/outfiles/vax_2.0_15.0_100_6_100_5 10 1000 2.0 15.0 100 6 100 5
set terminal pngcairo size 500,400 enhanced font 'Verdana,14'
set title "random network (expected mean degree=5)"
set output 'plots/vax_2.0_15.0_100_6_100_5.png'
plot 'outfiles/vax_2.0_15.0_100_6_100_5' using ($1*100):($4*100) w p pt 1 lc "black" ps 1.5 t "simulation", 'outfiles/vax_2.0_15.0_100_6_100_5' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case", 

#Here we don't see the benefit: why? 
#./bin/SIRvax SIRgeneral/outfiles/vax_2.0_50.0_100_7_100_20_4_2 100 100 2.0 50.0 100 7 100 20 4 2
set title "ring-preferential network"
set output 'plots/vax_ringpref.png'
plot 'outfiles/vax_2.0_50.0_100_7_100_20_4_2' using ($1*100):($4*100) w lp ls 1 t "ring-preferential network", 'outfiles/vax_2.0_50.0_100_7_100_20_4_2' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case",


#A complete network with smallish compartments: Theta = Theta_max
#./bin/SIRvax SIRgeneral/outfiles/outfiles/vax_2.0_19.8_20_1_100 10 100 2.0 19.8 20 1 100
set title "complete network"
set output 'plots/vax_2.0_19.8_20_1_100.png'
plot 'outfiles/vax_2.0_19.8_20_1_100' using ($1*100):($4*100) w lp ls 1 t "complete network", 'outfiles/vax_2.0_19.8_20_1_100' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case",


#A complete network with smallish compartments: Theta = Theta_max
#./bin/SIRvax SIRgeneral/outfiles/outfiles/vax_2.0_9.9_10_1_100 10 100 2.0 9.9 10 1 100
set title "complete network"
set output 'plots/vax_2.0_9.9_10_1_100.png'
plot 'outfiles/vax_2.0_9.9_10_1_100' using ($1*100):($4*100) w lp ls 1 t "ring-preferential network", 'outfiles/vax_2.0_9.9_10_1_100' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case",


#A complete network with smallish compartments: Theta = Theta_max
#./bin/SIRvax SIRgeneral/outfiles/outfiles/vax_2.0_9.8_10_1_50 10 100 2.0 9.8 10 1 50
set title "complete network"
set output 'plots/vax_2.0_9.8_10_1_50.png'
plot 'outfiles/vax_2.0_9.8_10_1_50' using ($1*100):($4*100) w lp ls 1 t "complete network", 'outfiles/vax_2.0_9.8_10_1_50' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case",

#A complete network with a single compartment
#./bin/SIRvax SIRgeneral/outfiles/outfiles/vax_2.0_0_500_1_1 10 100 2.0 0 500 1 1
set title "complete network"
set output 'plots/vax_2.0_0_500_1_1.png'
plot 'outfiles/vax_2.0_0_500_1_1' using ($1*100):($5*100) w lp ls 1 t "1 compartment", 'outfiles/vax_2.0_9.8_10_1_50' using ($1*100):($4*100) w lp ls 2 t "50 compartments", 'outfiles/vax_2.0_0_500_1_1' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case",




#./bin/SIRvax SIRgeneral/outfiles/vax_2.0_667_1000_6_50_8 10 1000 2.0 667 1000 6 50 8
#cp outfiles/network_properties outfiles/network_properties_2.0_667_1000_6_50_8
set title "random network"
set output 'plots/vax_2.0_667_1000_6_50_8.png'
plot 'outfiles/vax_2.0_667_1000_6_50_8' using ($1*100):($4*100) w lp ls 1 t "random network", 'outfiles/vax_2.0_667_1000_6_50_8' using ($1*100):($2*100) w l lc "black" lw 2 t "homogeneous case",

