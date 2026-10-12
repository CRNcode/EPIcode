#!/usr/bin/gnuplot -c 
set terminal pngcairo size 500,400 enhanced font 'Verdana,10'

# color definitions
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1 # --- blue
set style line 2 lc rgb '#ad6000' lt 1 lw 2 pt 7 ps 1 # --- red
set style line 6 lc rgb '#60ad00' lt 1 lw 2 pt 7 ps 1 # --- green
set style line 3 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2 # --- blue
set style line 4 lc rgb '#ad6000' lt 1 lw 2 pt 7 ps 2 # --- red
set style line 5 lc rgb '#000000' lt 1 lw 2 pt 7 ps 2 # --- black

unset key
#set key
unset xtics
unset ytics
unset cbtics
#set cblabel "density"
RepNo=ARG2
Nx=ARG3
Ny=ARG4
sampling=ARG5
leak=ARG6
Nc=ARG7
rewire=ARG8
set cbrange [0:1]


do for [i=0:ARG1] {
ii=i
unset title
set title sprintf('I, leak=%s, R_0=%s, N_c=%s, rewire=%s. Time=%04.02fs', leak, RepNo, Nc, rewire, i/sampling)
set output sprintf('pngheatI/%04.0f.png',i)

set view map
splot [-0.5:Nx-0.5] [-0.5:Ny-0.5] sprintf('dataI/m%04.0f',i) matrix using 1:2:3 with image

}

do for [i=0:ARG1] {
ii=i
unset title
set title sprintf('S, leak=%s, R_0=%s, N_c=%s, rewire=%s. Time=%04.02fs', leak, RepNo, Nc, rewire, i/sampling)
set output sprintf('pngheatS/%04.0f.png',i)

set view map
splot [-0.5:Nx-0.5] [-0.5:Ny-0.5] sprintf('dataS/m%04.0f',i) matrix using 1:2:3 with image
}

#Final susceptible state
set output sprintf('stills/sus_%s_%s_%s_%s_%s_%s.png',RepNo,leak,Nc,Nx,Ny,rewire)

do for [i=ARG1:ARG1] {
set view map
set title sprintf('final S, leak=%s, R_0=%s, N_c=%s, rewire=%s', leak, RepNo, Nc, rewire)
infile=sprintf('dataS/m%04.0f',i)
splot [-0.5:Nx-0.5] [-0.5:Ny-0.5] infile matrix using 1:2:3 with image
}

command = sprintf('cp %s finalS/sus_%s_%s_%s_%s_%s_%s', infile, RepNo,leak,Nc,Nx,Ny,rewire)
system(command)


set terminal pngcairo size 1000,400 enhanced font 'Verdana,10'
set xrange [0:ARG1/sampling]
set yrange [0:1.1]
set ytics ("max" 1.0)
set xtics
set lmargin 7
set bmargin 3

set label "infected" at -2,0.4 rotate left  font "Times,16"
set label "time" at 22,-0.1 font "Times,16"
set title sprintf('total infections')

n=0
do for [i=0:ARG1] {
n=n+1
ii=i
qq=(i<2000?0:i-2000)

    set output sprintf('pngtot/%04.0f.png',n)

    plot 'outfiles/SIRtorusmean' using ($1-qq*0.01):($2) every 1::0:0:ii:0 w l ls 1 notitle, \
         'outfiles/SIRtorusmean' using ($1-qq*0.01):($2) every ::ii:0:ii:0 w p ls 1 notitle,

}


