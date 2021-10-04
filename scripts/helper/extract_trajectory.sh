#!/bin/bash


function print_help {
   echo "extract_trajectory.sh needs a titania output file"
   echo "usage: extract_trajectory.sh output.out [-fgnu_format] [-cconfig]"
   echo "gnuformat: pdf (default)"
   echo "           png"
   echo "           svg"
   echo "config: file containing chirals centers for plot"
   exit
}
if [[ "$#" < "1" ]]
then
   print_help
fi

# Get the current directory
scriptdir="$(dirname $0)"
# Load the location of config file
TITANIA_GNUPLOT_PATH="$(dirname scriptdir)/gnuplot"

for arg in "$@"
do
   if [[ "${arg:0:1}" != "-" ]]
   then
      OUT_FILE=$arg
      continue
   fi
   case "${arg:1:1}" in
   "f")
      OUT_TYPE=${arg:2}
      ;;
   "i")
      OUT_FILE=${arg:2}
      ;;
   "c")
      configfile=${arg:2}
      echo "config file: $configfile"
      ;;
   *)
      echo "unknown argument: ${arg:0:2} ${arg:2}"
      ;;
   esac
done

if [ -z ${OUT_FILE+x} ]
then
   print_help
fi
TRJ_FILE=$OUT_FILE.trj

if [ -z ${OUT_TYPE+x} ]
then
   OUT_TYPE="pdf"
fi

# Parse the data to plot
# define offset between chiral volume header and actual values
OFFSET=2
HEADSET=1
OPT_STEPS=$(sed -n '/Finished iteration/h;g;$p' $TRJ_FILE | awk '{ print $5-1 }')
NOR=$(grep 'Number of RDCs:' $TRJ_FILE | awk '{print $4}')

if [ -z ${configfile+x} ]
then
   configfile=$OUT_FILE.tmp.conf
   awk -v OFFSET=$OFFSET 'BEGIN{x=0;}
   /Relative chiral volumes/{x=1}
   (x==1) {
      getline
      for ( i = 2; i <= NF; ++i ) {
         print $i
      }
      exit 0
   }
' $OUT_FILE > $configfile
fi

awk -v OFFSET=$OFFSET 'BEGIN{printf("%5s", "iter");x=0;iter=1;swaped=0}
(FNR==NR){
   center[NR]=$1
   printf(" %10s", center[NR])
   center[0]=NR
}
/Relative chiral volumes/{x=1}
(x==1) {
   getline
   if ( swaped == 0 ) {
   printf("\n")
   for ( i = 1; i <= NF; ++i ) {
         for ( j = 1; j <= center[0]; ++j ) {
            if ( center[j] == $i ) {
               center[j] = i;
            }
         }
      }
      swaped=1
   }
   getline
   printf("%5d", iter)
   for ( i = 1; i <= center[0]; ++i ) {
         printf(" %10s", $center[i])
   }
   printf("\n")
   x=0
   iter++
}
' $configfile $TRJ_FILE > $OUT_FILE.CV.dat


awk 'BEGIN{x=0;iter=1}
/MC stop criteria:/{
   i=0;
   crit[i] = $6;
   ++i;
   while ( getline && NF>1 ){
      crit[i]=$3;
      ++i;
      for ( j = 1; j < NF; ++j ) s[j]=$j
   }
   crit[i-1]=s[5]
   x=0;
   for ( j = 0; j < i; ++j ) {
      printf("%10.3e", crit[j]);
   }
   printf("\n");
   i = 0;
   ++iter;
}
/Motional information on RDC vector/{
getline
printf("%3d %10.6f ", iter, $3)
}
' $TRJ_FILE > $OUT_FILE.CRIT.dat

CENTER=$(wc -l $configfile | awk '{ print $1 }')

gnuplot -c $TITANIA_GNUPLOT_PATH/plot_trajectory.gnu $OUT_FILE $OPT_STEPS $CENTER $OUT_TYPE

rm -f $OUT_FILE.CV.dat $OUT_FILE.CRIT.dat

if [[ -f $OUT_FILE.tmp.conf ]]
then
   rm -f $OUT_FILE.tmp.conf
fi
