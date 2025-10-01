
print "first argument     : ", ARG1
set palette defined ( 0 "blue", 1 "red")
# splot "grid2.gnu"using 1:2:3:($4!=1?$4:1/0) palette pt 7

# splot "print_solution1.gnu" palette pt 7

String0 = 'plotting_data'
String1 = ARG1
String2 = ARG2
String3 = ".dat"


String = String0.String1.String2.String3

#set xtics nomirror
#set x2tics

#set autoscale xfix
#set autoscale x2fix
set xlabel 'DOFS'

#set x2label 'level'

#using 1:2:x2tic(1)
set multiplot layout 3,1
set logscale x
set logscale y
set format y "10^{%L}"
plot String  i 0 with linespoints  title "infty"
unset multiplot
unset format y


set terminal epslatex color colortex
set output "error_singularity1.tex"
plot String using 3:4 with linespoints notitle

set output "error_singularity2.tex"
plot String using 3:5 with linespoints notitle

set output "error_singularity3.tex"


plot String using 3:6 with linespoints notitle

set terminal qt
set output

#pause mouse close
