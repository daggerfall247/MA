set terminal pdfcairo


set xlabel "MC time"
set ylabel "<P>"

set yrange [0:1]

set grid

set output "plaquette_su3.pdf"
plot "su3_beta_6_heatbath.dat" using 0:1 w l title "heatbath", "su3_beta_6_metro.dat" using 0:1 w l title "metropolis"

set output "plaquette_su2.pdf"
plot "su2_beta_6_heatbath.dat" using 0:1 w l title "heatbath", "su2_beta_6_metro.dat" using 0:1 w l title "metropolis"
