reset

# build strings
cinput=sprintf("%s.CV.dat",ARG1)
sinput=sprintf("%s.CRIT.dat",ARG1)


# handle output
if ( ARG4 eq "pdf" ){
outname=sprintf("%s.trj.pdf",ARG1)
set term pdf enhanced color rounded lw 2.5 dl 2.00 size 10,10 font "Helvetica Bold,28"
set output outname
} else {
   if ( ARG4 eq "png" ){
   outname=sprintf("%s.trj.pdf",ARG1)
   print "ERROR: SECONDA currently can only be plotted as pdf."
   set term pdf enhanced color rounded lw 2.5 dl 2.00 size 10,10 font "Helvetica Bold,28"
   set output "Trajectory.pdf"
   } else {
      if ( ARG4 eq "svg" ){
      outname=sprintf("%s.trj.pdf",ARG1)
      set term pdf enhanced color rounded lw 2.5 dl 2.00 size 10,10 font "Helvetica Bold,28"
      set output "Trajectory.pdf"
      print "ERROR: Trajectory currently can only be plotted as pdf."
      }
   }
}


# enable \AA
set encoding iso_8859_1

# user functions
file_exists(file) = system("[ -f '".file."' ] && echo '1' || echo '0'") + 0


# style of chiral volume lines
if ( file_exists("config/cv_settings.gp") ) {
load 'config/cv_settings.gp'
} else {
set style line 10 lc rgb '#67001f' lw 1.5 dt 1
set style line 9  lc rgb '#b2182b' lw 2.5 dt 1
set style line 8  lc rgb '#d6604d' lw 1.5 dt 1
set style line 7  lc rgb '#f4a582' lw 1.5 dt 1
set style line 6  lc rgb '#fddbc7' lw 1.5 dt 1
set style line 5  lc rgb '#e0e0e0' lw 1.5 dt 1
set style line 4  lc rgb '#bababa' lw 1.5 dt 1
set style line 3  lc rgb '#878787' lw 1.5 dt 1
set style line 2  lc rgb '#4d4d4d' lw 1.5 dt 1
set style line 1  lc rgb '#1a1a1a' lw 1.5 dt 1
}


# format axes
set xzeroaxis ls -1
set xr[1:ARG2]
set yr[-1:1]
set zeroaxis
unset ylabel 
set format x ""


# format grid
set grid xtics ytics mxtics mytics lt 0 lw 1 lc rgb "#444444", lw 0.5 lc rgb "#bbbbbb"
set mytics 2
set mxtics 2


# format plot size
set lmargin at screen 0.15
set rmargin at screen 0.8
set multiplot layout 2,1


##
# Begin Chiral Volume plot
##

# set labels
set label 1 center at screen 0.02, 0.775 "Chiral Volume / {\\305}^3" rotate by 90 offset 1.5,0.0 font "Helvetica Bold, 36"

# format size
set tmargin at screen 0.95
set bmargin at screen 0.60

# format key
set key out width 0.8 autotitle columnhead right box 

#plot
plot for [i=2:(ARG3+1)] cinput u 1:i w lines ls (i-1) lw 2


##
# Begin Stop Criteria plot
##

#format axes
set xzeroaxis ls -1
set yr[1e-3:9.75]
set y2tics 0, 0.25
set y2r[0.0:1]
set zeroaxis
set format y "%6.0e"
set format x "%.0f"
set ytics nomirror
set y2tics tc rgb "#FF1111" nomirror
set logscale y

# set labels
set key center at screen 0.475, 0.08 horizontal box
set label 1 center at screen 0.02, 0.425 "Stop Criteria" rotate by 90 offset 1.5,0.0 font "Helvetica Bold, 36"
set label 2 center at screen 0.475, 0.15 "Iteration" font "Helvetica Bold, 36"

# format size
set tmargin at screen 0.60
set bmargin at screen 0.25

# format grid
set grid xtics ytics mxtics mytics lt 0 lw 1 lc rgb "#444444", lw 0.5 lc rgb "#bbbbbb"
set mytics 10
set mxtics

# plot
plot sinput u 1:3 with lines lc rgb "#1111FF" lw 2 t 'A',\
     ''     u 1:4 with lines lc rgb "#1111FF" lw 2 dt 3 t 'sig[A]',\
     ''     u 1:5 with lines lc rgb "#115F5F" lw 2 t 'p',\
     ''     u 1:6 with lines lc rgb "#116F5F" lw 2 dt 3 t 'sig[p]',\
     ''     u 1:7 with lines lc rgb "#111111" lw 2 t 'R^2',\
     ''     u 1:2 with lines axis x1y2 lc rgb "#FF1111" lw 2 t 'S_{overall}'

