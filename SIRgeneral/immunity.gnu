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
set lmargin 8
set bmargin 4
set xlabel "prior immunity (%)"
set ylabel "mean outbreak size (%)"

#set label "infected" at -2,0.4 rotate left  font "Times,16"
#set label "time" at 22,-0.1 font "Times,16"


#./bin/SIRimmunity 1000 1000 2.0 2.0 1000 1 10
set title "complete, symmetric network (n = 10, N_c = 1000)"
set output 'plots/complete_10_1000.png'
plot 'outfiles/vax_protect_complete_10_1000' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     	'outfiles/vax_protect_complete_10_1000' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_complete_10_1000' using ($1):($3) w p ls 2 t "random immunity", \
	 'outfiles/vax_2.0_2.0_1000_1_10' using ($1*100):($2*100) w l ls 2 t "theory",



set terminal pngcairo size 600,400 enhanced font 'Verdana,14'
#Complete, symmetric
#./bin/SIRimmunity 1000 10000 2.0 2.5 50 1 30
set title "complete, symmetric network"
set output 'plots/complete_2.0_2.5_50_1_30.png'
plot 'outfiles/vax_protect_complete_2.0_2.5_50_1_30' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     	'outfiles/vax_protect_complete_2.0_2.5_50_1_30' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_complete_2.0_2.5_50_1_30' using ($1):($3) w p ls 2 t "random immunity",

#Gamma
#./bin/SIRimmunity -c2 1000 10000 2.0 2.5 50 2 30 0.1 0
set title "Gamma network"
set output 'plots/gamma_2.0_2.5_50_2_30_0.1_0.png'
plot 'outfiles/vax_protect_gamma_2.0_2.5_50_2_30_0.1' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_gamma_2.0_2.5_50_2_30_0.1' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_gamma_2.0_2.5_50_2_30_0.1' using ($1):($3) w p ls 2 t "random immunity",

#Gamma symmetric
#./bin/SIRimmunity -c2 5000 10 1000 2.0 2.5 50 2 30 0.1 1
set title "symmetric Gamma network"
set output 'plots/gamma_2.0_2.5_50_2_30_0.1_1_c2.png'
plot 'outfiles/vax_protect_2.0_2.5_50_2_30_0.1_1_c2' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_2.5_50_2_30_0.1_1_c2' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_2.5_50_2_30_0.1_1_c2' using ($1):($3) w p ls 2 t "random immunity",


#Gamma (constant infectivity)
#./bin/SIRimmunity -c1 5000 10 1000 2.0 2.5 50 2 30 0.1 0
set title "Gamma network (constant infectivity)"
set output 'plots/gamma_2.0_2.5_50_2_30_0.1_0_c1.png'
plot 'outfiles/vax_protect_2.0_2.5_50_2_30_0.1_0_c1' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_2.5_50_2_30_0.1_0_c1' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_2.5_50_2_30_0.1_0_c1' using ($1):($3) w p ls 2 t "random immunity",

#Gamma (constant susceptibility)
#./bin/SIRimmunity -c3 5000 10 1000 2.0 2.5 50 2 30 0.1 0
set title "Gamma network (constant susceptibility)"
set output 'plots/gamma_2.0_2.5_50_2_30_0.1_0_c3.png'
plot 'outfiles/vax_protect_2.0_2.5_50_2_30_0.1_0_c3' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_2.5_50_2_30_0.1_0_c3' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_2.5_50_2_30_0.1_0_c3' using ($1):($3) w p ls 2 t "random immunity",


#ring
#./bin/SIRimmunity 1000 10000 2.0 4.0 50 3 100 4 0.1
set title "small-world ring"
set output 'plots/sw2.0_4.0_50_3_100_4_0.1.png'
plot 'outfiles/vax_protect_sw2.0_4.0_50_3_100_4_0.1' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_sw2.0_4.0_50_3_100_4_0.1' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_sw2.0_4.0_50_3_100_4_0.1' using ($1):($3) w p ls 2 t "random immunity",

#ring (constant infectivity)
#./bin/SIRimmunity -c1 5000 10 1000 2.0 4.0 50 3 100 4 0.1
set title "small-world ring (constant infectivity)"
set output 'plots/sw_2.0_4.0_50_3_100_4_0.1_c1.png'
plot 'outfiles/vax_protect_2.0_4.0_50_3_100_4_0.1_c1' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_4.0_50_3_100_4_0.1_c1' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_4.0_50_3_100_4_0.1_c1' using ($1):($3) w p ls 2 t "random immunity",

#ring (constant susceptibility)
#./bin/SIRimmunity -c3 5000 10 1000 2.0 4.0 50 3 100 4 0.1
set title "small-world ring (constant susceptibility)"
set output 'plots/sw_2.0_4.0_50_3_100_4_0.1_c3.png'
plot 'outfiles/vax_protect_2.0_4.0_50_3_100_4_0.1_c3' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_4.0_50_3_100_4_0.1_c3' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_4.0_50_3_100_4_0.1_c3' using ($1):($3) w p ls 2 t "random immunity",

#2D
#./bin/SIRimmunity 1000 10000 2.0 3.5 50 4 19 13 0.05
#set title "grid network"
set title "small-world grid"
set output 'plots/swgrid_2.0_3.5_50_4_19_13_0.05.png'
plot 'outfiles/vax_protect_2.0_3.5_50_4_19_13_0.05' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_3.5_50_4_19_13_0.05' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_3.5_50_4_19_13_0.05' using ($1):($3) w p ls 2 t "random immunity",

#2D (constant infectivity)
#./bin/SIRimmunity -c1 5000 10 1000 2.0 3.5 50 4 19 13 0.05
set title "small-world grid (constant infectivity)"
set output 'plots/swgrid_2.0_3.5_50_4_19_13_0.05_c1.png'
plot 'outfiles/vax_protect_2.0_3.5_50_4_19_13_0.05_c1' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_3.5_50_4_19_13_0.05_c1' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_3.5_50_4_19_13_0.05_c1' using ($1):($3) w p ls 2 t "random immunity",

#2D (constant susceptibility)
#./bin/SIRimmunity -c3 5000 10 1000 2.0 3.5 50 4 19 13 0.05
set title "small-world grid (constant susceptibility)"
set output 'plots/swgrid_2.0_3.5_50_4_19_13_0.05_c3.png'
plot 'outfiles/vax_protect_2.0_3.5_50_4_19_13_0.05_c3' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_3.5_50_4_19_13_0.05_c3' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_3.5_50_4_19_13_0.05_c3' using ($1):($3) w p ls 2 t "random immunity",

#random
#./bin/SIRimmunity 2000 10000 2.0 3.0 50 6 30 
set title "random network"
set output 'plots/random_2.0_3.0_50_6_30_8.png'
plot 'outfiles/vax_protect_2.0_3.0_50_6_30_8' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_3.0_50_6_30_8' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_3.0_50_6_30_8' using ($1):($3) w p ls 2 t "random immunity",

#random (constant infectivity)
#./bin/SIRimmunity -c1 5000 10 1000 2.0 3.0 50 6 30 8
set title "random network (constant infectivity)"
set output 'plots/random_2.0_3.0_50_6_30_8_c1.png'
plot 'outfiles/vax_protect_2.0_3.0_50_6_30_8_c1' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_3.0_50_6_30_8_c1' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_3.0_50_6_30_8_c1' using ($1):($3) w p ls 2 t "random immunity",

#random (constant susceptibility)
#./bin/SIRimmunity -c3 5000 10 1000 2.0 3.0 50 6 30 8
set title "random network (constant susceptibility)"
set output 'plots/random_2.0_3.0_50_6_30_8_c3.png'
plot 'outfiles/vax_protect_2.0_3.0_50_6_30_8_c3' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_3.0_50_6_30_8_c3' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_3.0_50_6_30_8_c3' using ($1):($3) w p ls 2 t "random immunity",

#ring preferential
#./bin/SIRimmunity 1000 10000 2.0 2.5 5 7 200 10 4 1
set title "ring-preferential network"
set output 'plots/rsf_2.0_2.5_5_7_200_10_4_1.png'
plot 'outfiles/vax_protect_2.0_2.5_5_7_200_10_4_1' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_2.5_5_7_200_10_4_1' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_2.5_5_7_200_10_4_1' using ($1):($3) w p ls 2 t "random immunity",
#	 outfiles/vax_protect_2.0_2.5_5_7_200_10_4_1

#ring preferential (small compartments)
#./bin/SIRimmunity -c2 5000 50 1000 2.0 2.0 5 7 200 20 4 1
set title "ring-preferential network (N_c=5)"
set output 'plots/rsf_2.0_2.0_5_7_200_20_4_1_c2.png'
plot 'outfiles/vax_protect_2.0_2.0_5_7_200_20_4_1_c2' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_2.0_5_7_200_20_4_1_c2' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_2.0_5_7_200_20_4_1_c2' using ($1):($3) w p ls 2 t "random immunity",

#ring preferential (medium compartments)
#./bin/SIRimmunity -c2 5000 50 1000 2.0 2.5 10 7 200 20 4 1
set title "ring-preferential network (N_c=10)"
set output 'plots/rsf_2.0_2.5_10_7_200_20_4_1_c2.png'
plot 'outfiles/vax_protect_2.0_2.5_10_7_200_20_4_1_c2' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_2.5_10_7_200_20_4_1_c2' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_2.5_10_7_200_20_4_1_c2' using ($1):($3) w p ls 2 t "random immunity",

#ring preferential (larger compartments)
#./bin/SIRimmunity -c2 5000 50 1000 2.0 3.0 20 7 200 20 4 1
set title "ring-preferential network (N_c=20)"
set output 'plots/rsf_2.0_3.0_20_7_200_20_4_1_c2.png'
plot 'outfiles/vax_protect_2.0_3.0_20_7_200_20_4_1_c2' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_3.0_20_7_200_20_4_1_c2' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_3.0_20_7_200_20_4_1_c2' using ($1):($3) w p ls 2 t "random immunity",

#ring preferential (even larger compartments)
#./bin/SIRimmunity -c2 5000 50 1000 2.0 3.0 40 7 200 20 4 1
set title "ring-preferential network (N_c=40)"
set output 'plots/rsf_2.0_3.0_40_7_200_20_4_1_c2.png'
plot 'outfiles/vax_protect_2.0_3.0_40_7_200_20_4_1_c2' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_3.0_40_7_200_20_4_1_c2' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_3.0_40_7_200_20_4_1_c2' using ($1):($3) w p ls 2 t "random immunity",


#Preferential attachment to ring, constant infectivity
#./bin/SIRimmunity -c1 5000 50 1000 2.0 4.0 5 7 200 10 6 3
set title "ring preferential network (constant infectivity)"
set output 'plots/rsf_2.0_4.0_5_7_200_10_6_3_c1.png'
plot 'outfiles/vax_protect_2.0_4.0_5_7_200_10_6_3_c1' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_4.0_5_7_200_10_6_3_c1' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_4.0_5_7_200_10_6_3_c1' using ($1):($3) w p ls 2 t "random immunity",

#Preferential attachment to ring, constant susceptibility
#./bin/SIRimmunity -c3 5000 50 1000 2.0 4.0 5 7 200 10 6 3
set title "ring preferential network (constant susceptibility)"
set output 'plots/rsf_2.0_4.0_5_7_200_10_6_3_c3.png'
plot 'outfiles/vax_protect_2.0_4.0_5_7_200_10_6_3_c3' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_4.0_5_7_200_10_6_3_c3' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_4.0_5_7_200_10_6_3_c3' using ($1):($3) w p ls 2 t "random immunity",

#Preferential attachment to ring, constant susceptibility
#./bin/SIRimmunity -c3 5000 50 1000 2.0 3.0 5 7 200 10 6 3
set title "ring preferential network (constant susceptibility)"
set output 'plots/rsf_2.0_3.0_5_7_200_10_6_3_c3.png'
plot 'outfiles/vax_protect_2.0_3.0_5_7_200_10_6_3_c3' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_2.0_3.0_5_7_200_10_6_3_c3' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_3.0_5_7_200_10_6_3_c3' using ($1):($3) w p ls 2 t "random immunity",



set terminal pngcairo size 1000,500 enhanced font 'Verdana,14'
#./bin/SIRimmunity 10000 10000 2.0 1.5 5 5 500 4
set title "preferential attachment (n = 500, N_c = 5, leak=1.5, initial clique=4)"
set output 'plots/sf_2.0_1.5_5_5_500_4.png'
plot 'outfiles/vax_protect_sf_2.0_1.5_5_5_500_4' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
     'outfiles/vax_protect_sf_2.0_1.5_5_5_500_4' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_sf_2.0_1.5_5_5_500_4' using ($1):($3) w p ls 2 t "random immunity",

#Gamma with large variance
set title "Gamma network"
set output 'plots/gamma_2.0_6.0_10_2_100_0.01_0_c2.png'
plot 	 'outfiles/vax_protect_2.0_6.0_10_2_100_0.01_0_c2' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
	 'outfiles/vax_protect_2.0_6.0_10_2_100_0.01_0_c2' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_6.0_10_2_100_0.01_0_c2' using ($1):($3) w p ls 2 t "random immunity",


#preferential
#./bin/SIRimmunity 1000 10000 2.0 5.0 20 5 200 1
set title "preferential attachment (n = 200, N_c = 20, initial clique=1)"
set output 'plots/sf_200_20_1a.png'
plot 	 'outfiles/vax_protect_sf_200_20_1a' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
	 'outfiles/vax_protect_sf_200_20_1a' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_sf_200_20_1a' using ($1):($3) w p ls 2 t "random immunity",

#./bin/SIRimmunity 1000 1000 2.0 4.0 20 5 200 1
set title "preferential attachment (n = 200, N_c = 20, initial clique=1)"
set output 'plots/sf_200_20_1b.png'
plot 	 'outfiles/vax_protect_sf_200_20_1b' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
	 'outfiles/vax_protect_sf_200_20_1b' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_sf_200_20_1b' using ($1):($3) w p ls 2 t "random immunity",

##
## Individual based networks
##

set terminal pngcairo size 600,400 enhanced font 'Verdana,14'

#individual based ring preferential (variable infectivity)
#./bin/SIRimmunity -c2 5000 50 1000 2.0 1 1 7 200 10 2 2
set title "individual ring-preferential network"
set output 'plots/rsf_2.0_1_1_7_200_10_2_2_c2.png'
plot 	 'outfiles/vax_protect_2.0_1_1_7_200_10_2_2_c2' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
	 'outfiles/vax_protect_2.0_1_1_7_200_10_2_2_c2' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_1_1_7_200_10_2_2_c2' using ($1):($3) w p ls 2 t "random immunity",

#individual based ring preferential (constant infectivity)
#./bin/SIRimmunity -c1 5000 50 1000 2.0 1 1 7 200 10 4 4
set title "individual ring-preferential network (constant infectivity)"
set output 'plots/rsf_2.0_1_1_7_200_10_4_4_c1.png'
plot 	 'outfiles/vax_protect_2.0_1_1_7_200_10_4_4_c1' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
	 'outfiles/vax_protect_2.0_1_1_7_200_10_4_4_c1' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_1_1_7_200_10_4_4_c1' using ($1):($3) w p ls 2 t "random immunity",

#individual based ring preferential (constant susceptibility)
#./bin/SIRimmunity -c3 5000 50 1000 2.0 1 1 7 200 10 6 3
set title "individual ring-preferential (constant susceptibility)"
set output 'plots/rsf_2.0_1_1_7_200_10_6_3_c3.png'
plot 	 'outfiles/vax_protect_2.0_1_1_7_200_10_6_3_c3' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
	 'outfiles/vax_protect_2.0_1_1_7_200_10_6_3_c3' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_1_1_7_200_10_6_3_c3' using ($1):($3) w p ls 2 t "random immunity",

#individual based ring preferential (constant susceptibility)
#./bin/SIRimmunity -c3 5000 50 1000 2.0 1 1 7 200 10 4 4
set title "individual ring-preferential (constant susceptibility)"
set output 'plots/rsf_2.0_1_1_7_200_10_4_4_c3.png'
plot 	 'outfiles/vax_protect_2.0_1_1_7_200_10_4_4_c3' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
	 'outfiles/vax_protect_2.0_1_1_7_200_10_4_4_c3' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_1_1_7_200_10_4_4_c3' using ($1):($3) w p ls 2 t "random immunity",

#individual random (variable infectivity)
#./bin/SIRimmunity -c2 5000 50 1000 2.0 1 1 6 200 5
set title "individual random"
set output 'plots/random_2.0_1_1_6_200_5_c2.png'
plot 	 'outfiles/vax_protect_2.0_1_1_6_200_5_c2' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
	 'outfiles/vax_protect_2.0_1_1_6_200_5_c2' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_2.0_1_1_6_200_5_c2' using ($1):($3) w p ls 2 t "random immunity",

#individual grid (no rewiring)
#./bin/SIRimmunity -c2 5000 50 1000 3.0 1 1 4 15 15 0.0
set title "individual-based grid (no rewiring)"
set output 'plots/grid_3.0_1_1_4_15_15_0.0_c2.png'
plot 	 'outfiles/vax_protect_3.0_1_1_4_15_15_0.0_c2' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
	 'outfiles/vax_protect_3.0_1_1_4_15_15_0.0_c2' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_3.0_1_1_4_15_15_0.0_c2' using ($1):($3) w p ls 2 t "random immunity",

#individual grid (5% rewiring)
#./bin/SIRimmunity -c2 5000 50 1000 3.0 1 1 4 15 15 0.05
set title "individual-based grid (5% rewiring)"
set output 'plots/grid_3.0_1_1_4_15_15_0.05_c2.png'
plot 	 'outfiles/vax_protect_3.0_1_1_4_15_15_0.05_c2' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
	 'outfiles/vax_protect_3.0_1_1_4_15_15_0.05_c2' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_3.0_1_1_4_15_15_0.05_c2' using ($1):($3) w p ls 2 t "random immunity",

#individual grid (10% rewiring)
#./bin/SIRimmunity -c2 5000 50 1000 3.0 1 1 4 15 15 0.1
set title "individual-based grid (10% rewiring)"
set output 'plots/grid_3.0_1_1_4_15_15_0.1_c2.png'
plot 	 'outfiles/vax_protect_3.0_1_1_4_15_15_0.1_c2' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
	 'outfiles/vax_protect_3.0_1_1_4_15_15_0.1_c2' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_3.0_1_1_4_15_15_0.1_c2' using ($1):($3) w p ls 2 t "random immunity",

#individual grid (20% rewiring)
#./bin/SIRimmunity -c2 5000 50 1000 3.0 1 1 4 15 15 0.2
set title "individual-based grid (20% rewiring)"
set output 'plots/grid_3.0_1_1_4_15_15_0.2_c2.png'
plot 	 'outfiles/vax_protect_3.0_1_1_4_15_15_0.2_c2' using ($1):($2):($1-$1):($3-$2) w vectors nohead lc rgb 'black' notitle, \
	 'outfiles/vax_protect_3.0_1_1_4_15_15_0.2_c2' using ($1):($2) w p ls 1 t "infection-acquired immunity", \
         'outfiles/vax_protect_3.0_1_1_4_15_15_0.2_c2' using ($1):($3) w p ls 2 t "random immunity",






