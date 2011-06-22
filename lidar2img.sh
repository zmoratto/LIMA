#!bin/sh

#image to LOLA matching 
mkdir ../data/Apollo15-PNG
for cubFile in `echo /Users/anefian/projects/LMMP/lima/data/Apollo15-CUB/*.cub`
do
   #echo "cubFile= $cubFile"
   pngFile=${cubFile/Apollo15-CUB/Apollo15-PNG}
   pngFile=${pngFile/cub/png}
   #echo "pngFile = $pnggFile" 
   #isis2std from= $cubFile to= $pngFile format= PNG
   #/Users/anefian/StereoPipeline-1.0.4-i386-OSX/bin/ipfind --num-threads 1 $pngFile -g 0.9 -l -d 
done

#new LOLA data
./lidar2img  -l ../data/Apollo15-LOLA/ApolloLOLA_22n24n_4e6e_pts_csv.csv ../data/Apollo15-CUB/*.cub





#old LOLA data
#./lidar2img  -l ../data/Apollo15-LOLA/1E8E_21N28N/RDR_1E8E_21N28NPointPerRow_csv_table.csv ../data/Apollo15-CUB/*.cub


