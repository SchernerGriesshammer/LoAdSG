# in Terminal:   gnuplot -c  plotdata.plt  uSolution2




reset
print "first argument     : ", ARG1

set ticslevel 0
# splot "grid2.gnu"using 1:2:3:($4!=1?$4:1/0) palette pt 7


# splot "print_solution1.gnu" palette pt 7

#String0 = "plotdata/"
#String1 = ARG1
#String2 = ARG2
#String3 = ".gnu"


#String = String0.String1.String2.String3

String = "error9.gnu"
unset key
unset colorbox
splot String pt 7 ps 0.2 lc  rgb '#0060ad'



pause mouse close
