#!/bin/bash
MAX=2000
# if no parameters input, go to "Usage is..."
if test $# -ne 2
  then
     echo "Usage: gen_gif <Size> <Max Iterations>"
     exit
  else
     SZ=$1
     MAX=$2
     echo $SZ
fi
#generate script to generate gnuplot formated output files using wave_exc
mkdir -p ./tmp/45_dat
echo '#!/bin/bash' > ./tmp/45_dat/gen.sh
for i in `seq 1 ${MAX}`
do
echo "./wave_exc ${SZ} ${i} ./tmp/45_dat/frame${i}" >> ./tmp/45_dat/gen.sh
done
#give execute permissions to the script
#and then run it 
chmod 755 ./tmp/45_dat/gen.sh
sh ./tmp/45_dat/gen.sh
#generate scripting file to be used by gnuplot
#to create gif file
echo "clear" > ./tmp/45_dat/plot.gpi
echo "reset" >> ./tmp/45_dat/plot.gpi
echo "set terminal gif animate delay 10" >> ./tmp/45_dat/plot.gpi
echo "set output \"wave.gif\"" >> ./tmp/45_dat/plot.gpi
echo "set hidden3d" >> ./tmp/45_dat/plot.gpi
echo "set zrange[-1.0:1.0]" >> ./tmp/45_dat/plot.gpi
echo "set yrange[0:${SZ}]" >> ./tmp/45_dat/plot.gpi
echo "set xrange[0:${SZ}]" >> ./tmp/45_dat/plot.gpi
echo "set view 36,41,1,1" >> ./tmp/45_dat/plot.gpi

for i in `seq 1 ${MAX}`
do
echo "splot \"./tmp/45_dat/frame${i}.gdt\" using 1:2:3 with lines" >> ./tmp/45_dat/plot.gpi
echo "!echo -n ." >> ./tmp/45_dat/plot.gpi
done
gnuplot ./tmp/45_dat/plot.gpi
/bin/rm -r ./tmp/45_dat
