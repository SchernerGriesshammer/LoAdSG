reset
print "first argument     : ", ARG1
set palette defined ( 0 "blue", 1 "red")
# splot "grid2.gnu"using 1:2:3:($4!=1?$4:1/0) palette pt 7

# splot "print_solution1.gnu" palette pt 7
unset xtics
unset ytics
String0 = 'grid'
String1 = ARG1
String2 = ARG2
String3 = ".gnu"


set key font ",14"
String = String0.String1.String2.String3




set size ratio -1

set key outside
#set key horizontal left

  set xrange [0:1]

  set yrange [0:1]

#1:2:3:($3>0.5?$2:1/0)
#set object 1 rect from 0,0 to 0.25,0.25 lw 1.2 fs empty border lc rgb 'grey'
#set object 2 rect from 0.75,0.75 to 1.0,1.0 lw 1.2 fs empty border lc rgb 'grey'

pss=0.2;


#splot String using 1:2:3:($4 == 1 ?$2:1/0)  pt 7 ps 0.3 linecolor "blue" title " active nodes",# String using 1:3:($4==0 ? $2 : 1/0) pt 7 ps 0.3 linecolor "red" title "hanging nodes"
plot 'grid6.gnu' pt 7 ps pss linecolor rgb '#0060ad' notitle, \
#'fillupgrid.gnu' i 0 using 1:($3==1 ? $2 : 1/0) pt 7 ps pss linecolor rgb '#0060ad' notitle, \
    # 'regulargrid1_addpoint.gnu'  i 0 using 1:($3==1 ? $2 : 1/0) pt 6 ps 4 lw 2 linecolor rgb 'red' notitle,\
    # 'regulargrid1_addpoint.gnu'  i 0 using 1:($3==1 ? $2 : 1/0) pt 7 ps pss linecolor rgb 'red' notitle,\
#	0/0  with points pt 7 ps 1.0 linecolor rgb '#0060ad' title "active nodes" ,\
#	0/0  with points pt 7 ps 1.0 linecolor rgb 'red' title "hanging nodes"

 # '< echo 0.0625 0.0625 0.02' w circ notitle lc "grey" lw 1.2,\
 #'< echo 0.9375 0.9375 0.02' w circ notitle lc "grey" lw 1.2,\

set terminal epslatex color colortex
set output "grid_singularity_2.tex"
replot
set terminal qt
set output

pause mouse close
