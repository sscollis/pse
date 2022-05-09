plot "space.out" using 1:2 with lines   
replot "inflow.+01" using 1:2 with lines
replot "space.out" using 1:3 with lines   
replot "inflow.+01" using 1:3 with lines

plot "space.out" using 1:4 with lines   
replot "inflow.+01" using 1:4 with lines
replot "space.out" using 1:5 with lines   
replot "inflow.+01" using 1:5 with lines

plot "space.out" using 1:6 with lines   
replot "inflow.+01" using 1:6 with lines
replot "space.out" using 1:7 with lines   
replot "inflow.+01" using 1:7 with lines

plot "space.out" using 1:8 with lines   
replot "inflow.+01" using 1:8 with lines
replot "space.out" using 1:9 with lines   
replot "inflow.+01" using 1:9 with lines

