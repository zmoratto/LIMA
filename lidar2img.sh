#!bin/sh

#image to LOLA matching 

#run lidar2img single thread
./lidar2img  -l ../data/Apollo15-LOLA/ApolloLOLA_22n24n_4e6e_pts_csv.csv -m ../data/Apollo15-MAP-CUB -d ../data/Apollo15-DRG camCubList.txt 

#./lidar2img  -l ../data/Apollo15-LOLA/ApolloLOLA_0n2n_0e2e_pts_csv.csv -m ../data/Apollo15-MAP-CUB camCubList.txt 

#run lidar2img multi-thread
echo /byss/moon/lola_tracks_over_apollo/*.csv | xargs -n 1 -I {} -P 8 ./lidar2img  -l {} -m ../data/Apollo15-MAP-CUB /byss/moon/apollo_metric/cubes/a15/sub4_cubes/*.lev1.cub
