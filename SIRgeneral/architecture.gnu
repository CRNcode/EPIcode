set terminal pngcairo size 500,400 enhanced font 'Verdana,14'


# color definitions
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1 # --- blue
set style line 2 lc rgb '#ad6000' lt 1 lw 2 pt 5 ps 1 # --- red
set style line 6 lc rgb '#60ad00' lt 1 lw 2 pt 7 ps 1 # --- green
set style line 3 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 4 lc rgb '#ad6000' lt 1 lw 2 pt 7 ps 2 # --- red
set style line 5 lc rgb '#000000' lt 1 lw 2 pt 7 ps 2 # --- black

set xrange [-0.1:2.1]
set ytics 5
set lmargin 8
set bmargin 4
set xlabel "-log10(α)"
set ylabel "mean outbreak size (%)"
set style fill  transparent solid 0.2 noborder
set style circle radius 0.02

#Gamma (variable infectivity)
#./bin/SIRgamma -c2 100 1000 2.0 2.5 50 30 0
set title "Gamma networks"
set output 'plots/gamma_effect_2.0_2.5_50_30_0.png'
plot 'outfiles/gamma_effect_2.0_2.5_50_30_0' using (-log10($1)):($4) with circles lc rgb "blue" t "each grid", \
'outfiles/gamma_effect_mean_2.0_2.5_50_30_0' using (-log10($1)):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all",

#Gamma (variable infectivity, symmetric)
#./bin/SIRgamma -c2 100 1000 2.0 2.5 50 30 1
set title "symmetric Gamma networks"
set output 'plots/gamma_effect_2.0_2.5_50_30_1_c2.png'
plot 'outfiles/gamma_effect_2.0_2.5_50_30_1_c2' using (-log10($1)):($4) with circles lc rgb "blue" t "each grid", \
'outfiles/gamma_effect_mean_2.0_2.5_50_30_1_c2' using (-log10($1)):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all", 

#Gamma (constant infectivity)
#./bin/SIRgamma -c1 100 1000 2.0 2.5 50 30 0
set title "Gamma networks (constant infectivity)"
set output 'plots/gamma_effect_2.0_2.5_50_30_0_c1.png'
plot 'outfiles/gamma_effect_2.0_2.5_50_30_0_c1' using (-log10($1)):($4) with circles lc rgb "blue" t "each grid", \
'outfiles/gamma_effect_mean_2.0_2.5_50_30_0_c1' using (-log10($1)):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all",

#Gamma (constant susceptibility)
#./bin/SIRgamma -c3 100 1000 2.0 2.5 50 30 0
set title "Gamma networks (constant susceptibility)"
set output 'plots/gamma_effect_2.0_2.5_50_30_0_c3.png'
plot 'outfiles/gamma_effect_2.0_2.5_50_30_0_c3' using (-log10($1)):($4) with circles lc rgb "blue" t "each grid", \
'outfiles/gamma_effect_mean_2.0_2.5_50_30_0_c3' using (-log10($1)):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all",

set style circle radius 0.1
set key bottom right
#./bin/SIRrandom -c1 100 1000 2.0 3.0 50 100
set title "random networks (constant infectivity)"
set xrange [3.6:14.4]
set yrange [0:25]
set xlabel "mean degree"
set output 'plots/random_meandeg_2.0_3.0_50_100_c1.png'
plot 'outfiles/random_meandeg_2.0_3.0_50_100_c1' using ($1):($4) with circles lc rgb "blue" t "each grid", \
         'outfiles/random_meandeg_mean_2.0_3.0_50_100_c1' using ($1):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all",

#./bin/SIRrandom -c3 100 1000 2.0 3.0 50 100
set title "random networks (constant susceptibility)"
set xrange [3.6:14.4]
set yrange [0:25]
set xlabel "mean degree"
set output 'plots/random_meandeg_2.0_3.0_50_100_c3.png'
plot 'outfiles/random_meandeg_2.0_3.0_50_100_c3' using ($1):($4) with circles lc rgb "blue" t "each grid", \
         'outfiles/random_meandeg_mean_2.0_3.0_50_100_c3' using ($1):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all",


#./bin/SIRrandom -c2 100 1000 2.0 2.5 50 100
set title "random networks"
#unset title
#unset ylabel
set xrange [3.6:14.4]
set yrange [0:20]
set xlabel "mean degree"
set output 'plots/random_meandeg_2.0_2.5_50_100.png'
plot 'outfiles/random_meandeg_2.0_2.5_50_100' using ($1):($4) with circles lc rgb "blue" t "each grid", \
         'outfiles/random_meandeg_mean_2.0_2.5_50_100' using ($1):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all",



set yrange [0:20]
set style circle radius 0.2


#./bin/SIRsfring 100 1000 2.0 4.0 10 200
set title "ring-preferential networks"
#unset ylabel
set xrange [3.5:24.5]
set xlabel "initial ring size"
set output 'plots/ringsf_2.0_4.0_10_200.png'
plot 'outfiles/outsize_ringsf_2.0_4.0_10_200_rewire' using ($1):($2) with circles lc rgb "blue" t "each grid", \
         'outfiles/outsize_ringsf_2.0_4.0_10_200_rewire_mean' using ($1):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all"

set key top left

#./bin/SIRring 100 1000 2.0 4.0 50 100 4
set xrange [-0.5:20.5]
set title "small-world rings"
set xlabel "rewiring probability"
#unset title
set output 'plots/sw_rewire_2.0_4.0_50_100_4.png'
plot 'outfiles/sw_rewire_2.0_4.0_50_100_4' using (100*($1)):($2) with circles lc rgb "blue" t "each grid", \
         'outfiles/sw_rewire_mean_2.0_4.0_50_100_4' using (100*($1)):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all",

#./bin/SIRring -c1 100 1000 2.0 4.0 50 100 4
set xrange [-0.5:20.5]
set title "small-world rings (constant infectivity)"
set xlabel "rewiring probability"
#unset title
set output 'plots/sw_rewire_2.0_4.0_50_100_4_c1.png'
plot 'outfiles/sw_rewire_2.0_4.0_50_100_4_c1' using (100*($1)):($2) with circles lc rgb "blue" t "each grid", \
         'outfiles/sw_rewire_mean_2.0_4.0_50_100_4_c1' using (100*($1)):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all",

#./bin/SIRring -c3 100 1000 2.0 4.0 50 100 4
set xrange [-0.5:20.5]
set title "small-world rings (constant susceptibility)"
set xlabel "rewiring probability"
#unset title
set output 'plots/sw_rewire_2.0_4.0_50_100_4_c3.png'
plot 'outfiles/sw_rewire_2.0_4.0_50_100_4_c3' using (100*($1)):($2) with circles lc rgb "blue" t "each grid", \
         'outfiles/sw_rewire_mean_2.0_4.0_50_100_4_c3' using (100*($1)):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all",

#./bin/SIRgrid 100 1000 2.0 3.5 30 19 13
set title "small-world grids"
#unset title
#unset ylabel
set output 'plots/grid_rewire_2.0_3.5_30_19_13.png'
plot 'outfiles/grid_rewire_2.0_3.5_30_19_13' using (100*($1)):($2) with circles lc rgb "blue" t "each grid", \
         'outfiles/grid_rewire_mean_2.0_3.5_30_19_13' using (100*($1)):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all",


#./bin/SIRgrid -c1 100 1000 2.0 3.5 30 19 13
set title "small-world grids (constant infectivity)"
#unset title
#unset ylabel
set output 'plots/grid_rewire_2.0_3.5_30_19_13_c1.png'
plot 'outfiles/grid_rewire_2.0_3.5_30_19_13_c1' using (100*($1)):($2) with circles lc rgb "blue" t "each grid", \
         'outfiles/grid_rewire_mean_2.0_3.5_30_19_13_c1' using (100*($1)):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all",

#./bin/SIRgrid -c3 100 1000 2.0 3.5 30 19 13
set title "small-world grids (constant susceptibility)"
#unset title
#unset ylabel
set output 'plots/grid_rewire_2.0_3.5_30_19_13_c3.png'
plot 'outfiles/grid_rewire_2.0_3.5_30_19_13_c3' using (100*($1)):($2) with circles lc rgb "blue" t "each grid", \
         'outfiles/grid_rewire_mean_2.0_3.5_30_19_13_c3' using (100*($1)):($2) w circles lc rgb "red" fs solid 1.0 border lt -1 t "mean of all",
	 


