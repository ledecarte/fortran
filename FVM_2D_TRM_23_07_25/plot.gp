# First: PNG Output
set terminal pngcairo size 800,600 enhanced font 'Arial,10'
set output 'fvm_temperature_ansys_rainbow.png'

set xlabel "x [m]"
set ylabel "y [m]"
set title "2D Temperature Distribution (FVM)"
set cblabel "Temperature (K)"

set cbrange [300.0:3500]


set palette defined (\
  0   "#0000ff", \
  0.2 "#00ffff", \
  0.4 "#00ff00", \
  0.6 "#ffff00", \
  0.8 "#ff8000", \
  1.0 "#ff0000" \
)

set view map
set size ratio -1
unset key
set pm3d map interpolate 2,2

splot 'fvm.dat' using 1:2:3 with pm3d

# Then: On-screen Display
unset output
set terminal qt persist
replot
