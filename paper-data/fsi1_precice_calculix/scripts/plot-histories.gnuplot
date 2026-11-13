# gnuplot scripts/plot-histories.gnuplot   (run from paper-data/fsi1_precice_calculix)
set terminal pngcairo size 1400,900 font ",11"
set output "fsi1_precice_histories.png"
set datafile separator ","
set multiplot layout 2,2 title "Hron-Turek FSI1, OpenFOAM + CalculiX via preCICE (dashed: Featflow 7+0)"
set xlabel "t (s)"; set grid; set key bottom right; set xrange [0:30]
L1="histories/level1-rbf.csv"; L2="histories/level2.csv"; L4="histories/level4.csv"
set ylabel "u_x(A) (mm)"; set yrange [0.015:0.03]
plot L1 u 1:($2*1e3) every ::1 w l lw 2 t "1x", L2 u 1:($2*1e3) every ::1 w l lw 2 t "2x", L4 u 1:($2*1e3) every ::1 w l lw 2 t "4x", 0.02270493 dt 2 lc rgb "black" t "Featflow"
set ylabel "u_y(A) (mm)"; set yrange [0:0.9]
plot L1 u 1:($3*1e3) every ::1 w l lw 2 t "1x", L2 u 1:($3*1e3) every ::1 w l lw 2 t "2x", L4 u 1:($3*1e3) every ::1 w l lw 2 t "4x", 0.8208773 dt 2 lc rgb "black" t "Featflow"
set ylabel "drag (N/m)"; set yrange [13.5:15.5]
plot L1 u 1:4 every ::1 w l lw 2 t "1x", L2 u 1:4 every ::1 w l lw 2 t "2x", L4 u 1:4 every ::1 w l lw 2 t "4x", 14.29426 dt 2 lc rgb "black" t "Featflow"
set ylabel "lift (N/m)"; set yrange [-0.2:1.0]
plot L1 u 1:5 every ::1 w l lw 2 t "1x", L2 u 1:5 every ::1 w l lw 2 t "2x", L4 u 1:5 every ::1 w l lw 2 t "4x", 0.763746 dt 2 lc rgb "black" t "Featflow"
