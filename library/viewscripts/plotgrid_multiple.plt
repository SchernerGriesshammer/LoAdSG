reset

set terminal qt size 500,500

set size ratio 1

unset key
#set key horizontal left

  set xrange [0:1]

  set yrange [0:1]
unset xtics
unset ytics
set style rect fs empty border lc rgb '#ff1c08' lw 2.0 dt 2

set multiplot layout 4,4 columns

set lmargin 0.5
set rmargin 0.5
set tmargin 0.5
set bmargin 0.5

set format ''
set multiplot next
set multiplot next
set multiplot next
set multiplot next

#set multiplot next #plot '../results/gridpoints15.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#0000ff' notitle
  #plot '../results/gridpoints14.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#0000ff' notitle
set ylabel "t_y = 3" font ",18"
plot '../results/gridpointsALL.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#bfbfbf' notitle, '../results/gridpoints13.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#0000ff' notitle
unset ylabel
set ylabel "t_y = 2" font ",18"
plot '../results/gridpointsALL.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#bfbfbf' notitle, '../results/gridpoints12.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#0000ff' notitle
unset ylabel
set xlabel "t_x = 1" font ",18"
set ylabel "t_y = 1" font ",18"
plot '../results/gridpointsALL.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#bfbfbf' notitle, '../results/gridpoints11.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#0000ff' notitle
unset xlabel
unset ylabel
set multiplot next
#set multiplot next

  #plot '../results/gridpointsALL.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#bfbfbf' notitle, '../results/gridpoints24.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#0000ff' notitle
plot '../results/gridpointsALL.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#bfbfbf' notitle, '../results/gridpoints23.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#0000ff' notitle
  
plot '../results/gridpointsALL.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#bfbfbf' notitle, '../results/gridpoints22.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#0000ff' notitle
set xlabel "t_x = 2" font ",18"
plot '../results/gridpointsALL.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#bfbfbf' notitle, '../results/gridpoints21.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#0000ff' notitle
unset xlabel
set multiplot next
#set multiplot next

 #plot '../results/gridpointsALL.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#bfbfbf' notitle, '../results/gridpoints34.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#0000ff' notitle
set object 1 rect from screen 0.7525,0.7525 to 1.05,1.05 front
plot '../results/gridpointsALL.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#bfbfbf' notitle, '../results/gridpoints33.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#0000ff' notitle
unset object 1
plot '../results/gridpointsALL.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#bfbfbf' notitle, '../results/gridpoints32.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#0000ff' notitle
set xlabel "t_x = 3" font ",18"
plot '../results/gridpointsALL.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#bfbfbf' notitle, '../results/gridpoints31.gnu' using 1:($3==1 ? $2 : 1/0) pt 7 ps 0.6 linecolor rgb '#0000ff' notitle
unset xlabel
set multiplot next


unset multiplot


set terminal epslatex color colortex
set output "grid_singularity.tex"
replot
set terminal qt
set output

pause mouse close
