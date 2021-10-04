#!/bin/bash
gnuplot=/usr/bin/gnuplot
output_pdf=1
output_svg=0
output_png=0
reduce_basename=0
titania_output=""
basename=""
plot_file=""
data_file=".gnu_S_data.dat"
gnu_file="plot_Sdata.gnu"
function parse_options
{
   for (( i=0; i<${#options}; ++i )); do
      opt=${options:$i:1}
      case $opt in
         p)
            output_pdf=1
            ;;
         s)
            output_svg=1
            output_pdf=0
            ;;
         g)
            output_png=1
            output_pdf=0
            ;;
         r)
            reduce_basename=1
            ;;
        \?)
            echo "Invalid Option: -" >&2
            ;;
      esac
   done
}
# Parse the arguments
for opt in $@; do
   if [[ $opt == -* ]]; then
      options=${opt:1:${#opt}}
   else
      titania_output=$opt
   fi
done
# Load options
parse_options
basename=$titania_output
if [ $reduce_basename -eq 1 ] ; then
   basename=`basename $titania_output`
fi
plot_file="$basename.dyn"
# Parse the data to plot
NUM_RDCS=$(sed -n "$((`sed -n '/Number of RDCs/=' $titania_output`))p" $titania_output | awk '{print $4}')
XSIZE=$(($(($NUM_RDCS/3)) + 1 ))
STARTPOINT=$(sed -n '/Motional information on RDC vectors:/=' $titania_output | tac | head -n 1)
SOVERALL=$(sed -n "$(($STARTPOINT+1))p" $titania_output | awk '{print $3}')
sed -n "$(($STARTPOINT+3)),$(($STARTPOINT+2+$NUM_RDCS))p" $titania_output | awk '{print $1"-"$2"	"$3"	"$4"	"100*$5}' > $data_file
# Generate gnuplot script
echo "#!$gnuplot" > $gnu_file
echo "reset" >> $gnu_file
if [ $output_pdf -eq 1 ] ; then
   echo "set term pdf enhanced color rounded lw 1.5 dl 1.25 size $XSIZE,10 font \"Helvetica Bold,20\"" >> $gnu_file
   echo "set output '$plot_file.pdf'" >> $gnu_file
elif [ $output_svg -eq 1 ] ; then 
   echo "set term svg" >> $gnu_file
   echo "set output '$plot_file.svg'" >> $gnu_file
else
   echo "set term png" >> $gnu_file
   echo "set output '$plot_file.png'" >> $gnu_file
fi
echo "set key off" >> $gnu_file
echo "set yrange [0:1]" >> $gnu_file
echo "set xtics nomirror out rotate by 60 right nooffset" >> $gnu_file
echo "set ytics nomirror out" >> $gnu_file
echo "set multiplot layout 2,1" >> $gnu_file
echo "set ylabel 'S^2'" >> $gnu_file
echo "set key below" >> $gnu_file
echo "set style histogram gap 1" >> $gnu_file
echo "set arrow 1 from graph 0, $SOVERALL to graph 1, $SOVERALL nohead front linecolor rgb \"red\"" >> $gnu_file
echo "set label 2 at graph 0.95, $SOVERALL+0.04 \"S_{overall}\" front center" >> $gnu_file
echo "plot '$data_file' using 2:xtic(1) with histogram fill pattern 7 title 'S^2_{rdc}', '' u 3 with histogram fill pattern 6 title 'S^2_{axial}'" >> $gnu_file
echo "set ylabel '{/Symbol h} / %'" >> $gnu_file
echo "set yrange [0:*]" >> $gnu_file
echo "unset arrow 1" >> $gnu_file
echo "unset label 2" >> $gnu_file
echo "set xlabel 'RDC pair" >> $gnu_file
echo "set boxwidth 2" >> $gnu_file
echo "unset key" >> $gnu_file
echo "set style fill pattern 7 border lt -1" >> $gnu_file
echo "plot '$data_file' using 4:xtic(1) with histogram" >> $gnu_file
echo "unset ylabel" >> $gnu_file
echo "unset xlabel" >> $gnu_file
$gnuplot $gnu_file
# Clean
rm $data_file
rm $gnu_file
