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

#run lidar2img
./lidar2img  -l ../data/Apollo15-LOLA/ApolloLOLA_22n24n_4e6e_pts_csv.csv -m ../data/Apollo15-MAP-CUB ../data/Apollo15-CAM-CUB/*.lev1.cub


