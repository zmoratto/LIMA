#!bin/sh

#image to LOLA matching 

#run lidar2img single processor
./lidar2img  -l ../data/Apollo15-LOLA/ApolloLOLA_22n24n_4e6e_pts_csv.csv -m ../data/Apollo15-MAP-CUB camCubList.txt 




#run lidar2img multi processor
#echo /byss/moon/lola_tracks_over_apollo/*.csv |xargs -n 1 -I {} -P 2 ./lidar2img  -l {} -m ../data/Apollo15-MAP-CUB -c ../data/Apollo15-CAM-CUB/*.lev1.cub
#echo ../data/Apollo15-LOLA/*.csv |xargs -n 1 -I {} -P 2 ./lidar2img  -l {} -m ../data/Apollo15-MAP-CUB camCubList.txt