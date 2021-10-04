#!/bin/bash
gnuplot=/usr/bin/gnuplot
output_pdf=1
output_svg=0
output_png=0
reduce_basename=0
titania_output=""
basename=""
plot_file=""
data_file=".gnu_mcp_data.dat"
gnu_MC_file="plot_mcp_data.gnu"

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

plot_file="$basename.mcp"

# Parse the data to plot
sed -n "$((`sed -n '/Monte Carlo output:/=' $titania_output`+1)),$ p" $titania_output > $data_file
NUM_RDCS=$(sed -n "$((`sed -n '/Number of RDCs/=' $titania_output`))p" $titania_output | awk '{print $4}')

if [ $(($NUM_RDCS % 2)) ]
then
   HALF_RDCS=$((1+$NUM_RDCS/2))
else
   HALF_RDCS=$(($NUM_RDCS/2))
fi
XSIZE=$(($HALF_RDCS*5))

if [ $output_pdf -eq 1 ] ; then
   out1="set term pdf enhanced color rounded lw 1.5 dl 1.25 size 15,$XSIZE font \"Helvetica Bold,20\""
   out2="set output '$plot_file.pdf'"
elif [ $output_svg -eq 1 ] ; then 
   out1="set term svg"
   out2="set output '$plot_file.svg'"
else
   out1="set term png"
   out2="set output '$plot_file.png'"
fi

# Generate gnuplot script
echo "reset" > $gnu_MC_file
echo "$out1" >> $gnu_MC_file
echo "$out2" >> $gnu_MC_file
echo "set encoding iso_8859_1" >> $gnu_MC_file
echo "set polar" >> $gnu_MC_file
echo "set style line 10 lt 1 lc 0 lw 0.3" >> $gnu_MC_file
echo "set grid polar pi/6 front" >> $gnu_MC_file
echo "set grid ls 10" >> $gnu_MC_file
echo "unset border" >> $gnu_MC_file
echo "set xrange [-pi:pi]" >> $gnu_MC_file
echo "set yrange [-pi:pi]" >> $gnu_MC_file
echo "set xtics axis" >> $gnu_MC_file
echo "set ytics axis" >> $gnu_MC_file
echo "unset tics" >> $gnu_MC_file
echo "unset raxis" >> $gnu_MC_file
echo "set rtics (\"45\" pi/4, \"90\" pi/2, \"135\" 3*pi/4, \"180\" pi)" >> $gnu_MC_file
echo "set size square" >> $gnu_MC_file
echo "set key lmargin" >> $gnu_MC_file
echo "set_label(x, text) = sprintf(\"set label '%s' at (3.4*cos(%f)), (3.4*sin(%f)) front center\", text, x, x) #this places a label on the outside" >> $gnu_MC_file

echo "eval set_label(0, \"90\")" >> $gnu_MC_file
echo "eval set_label(pi/6, \"120\")" >> $gnu_MC_file
echo "eval set_label(pi/3, \"150\")" >> $gnu_MC_file
echo "eval set_label(pi/2, \"180\")" >> $gnu_MC_file
echo "eval set_label(2*pi/3, \"-150\")" >> $gnu_MC_file
echo "eval set_label(5*pi/6, \"-120\")" >> $gnu_MC_file
echo "eval set_label(pi, \"-90\")" >> $gnu_MC_file
echo "eval set_label(7*pi/6, \"-60\")" >> $gnu_MC_file
echo "eval set_label(8*pi/6, \"-30\")" >> $gnu_MC_file
echo "eval set_label(9*pi/6, \"0\")" >> $gnu_MC_file
echo "eval set_label(10*pi/6, \"30\")" >> $gnu_MC_file
echo "eval set_label(11*pi/6, \"60\")" >> $gnu_MC_file
echo "eval set_label(2*pi, \"90\")" >> $gnu_MC_file
echo "set xzeroaxis ls -1" >> $gnu_MC_file
echo "set xrange [-pi:pi]" >> $gnu_MC_file
echo "set yrange [-pi:pi]" >> $gnu_MC_file

echo "set multiplot layout $HALF_RDCS,2 scale 1,0.95" >> $gnu_MC_file
echo "set key autotitle columnhead outside" >> $gnu_MC_file
echo "do for[i=1:$NUM_RDCS]{" >> $gnu_MC_file
echo "   p=2*i;" >> $gnu_MC_file
echo "   t=p-1;" >> $gnu_MC_file
echo "   plot '$data_file' using (column(p)-pi/2):t" >> $gnu_MC_file
echo "}" >> $gnu_MC_file
echo "unset multiplot" >> $gnu_MC_file
echo "unset ylabel" >> $gnu_MC_file
echo "unset xlabel" >> $gnu_MC_file

$gnuplot $gnu_MC_file

# Clean
rm $data_file
rm $gnu_MC_file
