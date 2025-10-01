set terminal pdfcairo enhanced
set output '6D_infty.pdf'

N_end = 10
N=1
set key font ",14"
set ytics font ", 12"

set style line 1 \
    linecolor rgb '#0060ad' \
    linetype 3 linewidth 1.5 \
    pointtype 5 ps 1\
    	dt 5\
  
set style line 2 \
    linecolor rgb '#ff7f0e' \
    linetype 1 linewidth 1.5\
     pointtype 8 ps 1
  
  
#set xrange[5:25000]
 set xlabel "DOF" font ",20"
 set ylabel '         { e_{H_1} }' font ",20"
 #set ylabel '{H_1 error, L_{/Symbol \245} error }'


datafile = "data.txt"

# Clear datafile
system(sprintf("echo '' > %s", datafile))

# Calculate and print data points
r = 0.0
N = 0.0
set print datafile append
do for [n = 1:10] {
    x = (2**n)*(n**5)
    y = ((2**(-n)))**2*(n**5)
    print x, y
}

print " "
print " "

do for [n = 1:10] {
    x = (2**n)*(n**5)
    y = ((2**(-n)))*(n**5)
    print x, y
}

print " "
print " "

do for [n = 2:10] {
    x = (2**n)**2
    y = ((2**(-n))**2)
    print x, y
}

print " "
print " "

data='data.dat'
data2='data.txt'


unset print
set logscale y 10
set logscale x 10



set format y "10^{%L}"
set format x "10^{%L}"
set key bottom left

plot data i 0  with linespoints title 'const. coeff. regular',\
data i 1  with linespoints  title 'var. coeff. regular',\
data i 2  with linespoints  title 'const. coeff. adaptive',\
#data i 3  with linespoints  title 'var. coeff. adaptive',\
#data i 5  with linespoints  title 'h log(h)^{d-1}'
#data i 4  with linespoints  title 'h^2 log(h)^{d-1}'
#data i 5  with linespoints  title 'h log(h)^{d-1}'

 



pause mouse close
