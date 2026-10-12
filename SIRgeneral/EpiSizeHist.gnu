set terminal pngcairo size 500,400 enhanced font 'Verdana,12'
set xlabel "outbreak size (%)"
set ylabel "frequency"
set style fill solid 0.5
#stats "outfiles/OutbreakHist_2.0_3.5_30_4_19_13_0.10" using 1 nooutput
#set print $mean
#print STATS_mean

set lmargin 8
set bmargin 4
unset ytics
binwidth=2
bin(x,width)=width*floor(x/width)
set boxwidth binwidth
set xrange [2:80]
set yrange[0:2600]
M=1
L=0.95
unset key
binnumber(x) = floor(x/binwidth)
rounded(x) = binwidth * (binnumber(x) + 0.5 )

set output "plots/EpiSizeHistGrid_000.png"
set title "no rewiring"
#set arrow nohead from STATS_mean,0.0 to STATS_mean,700.0 lw 2
plot  "outfiles/OutbreakHist_2.0_3.5_30_4_19_13_0.00" using (rounded($1)):(1) smooth freq with boxes fs transparent t ""

set output "plots/EpiSizeHistGrid_002.png"
set title "2% rewiring"
#set arrow nohead from STATS_mean,0.0 to STATS_mean,700.0 lw 2
plot  "outfiles/OutbreakHist_2.0_3.5_30_4_19_13_0.02" using (rounded($1)):(1) smooth freq with boxes fs transparent t ""

set output "plots/EpiSizeHistGrid_005.png"
set title "5% rewiring"
#set arrow nohead from STATS_mean,0.0 to STATS_mean,700.0 lw 2
plot  "outfiles/OutbreakHist_2.0_3.5_30_4_19_13_0.05" using (rounded($1)):(1) smooth freq with boxes fs transparent t ""

set output "plots/EpiSizeHistGrid_010.png"
set title "10% rewiring"
#set arrow nohead from STATS_mean,0.0 to STATS_mean,700.0 lw 2
plot  "outfiles/OutbreakHist_2.0_3.5_30_4_19_13_0.10" using (rounded($1)):(1) smooth freq with boxes fs transparent t ""

set output "plots/EpiSizeHistGrid_020.png"
set title "20% rewiring"
#set arrow nohead from STATS_mean,0.0 to STATS_mean,700.0 lw 2
plot  "outfiles/OutbreakHist_2.0_3.5_30_4_19_13_0.20" using (rounded($1)):(1) smooth freq with boxes fs transparent t ""
