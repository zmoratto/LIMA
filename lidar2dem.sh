#!bin/sh
#TIFF1=`ls ../data/LROC-DTMcomparison/ASU/Apollo15_comp/apollo15_dem.tif`
#TIFF2=`ls ../data/LROC-DTMcomparison/DLR/DLR_ap_dtm.tif`
#TIFF3=`ls ../data/LROC-DTMcomparison/UofA/Apollo15_DTM.tif`
#TIFF4=`ls ../data/LROC-DTMcomparison/USGS/USGS_M_111_06_16_DEM_150cm.tif`
#TIFF=$TIFF1" "$TIFF2" "$TIFF3" "$TIFF4
#TIFF=`ls ../data/LROC-DTMcomparison/USGS/USGS_M_111_06_16_DEM_150cm.tif`
#echo $TIFF
A17_DEM=`ls -d -1 /Users/anefian/projects/LMMP/LMMP_delivery/2011_0523_Apollo17_DEM/*.tif`
echo $A17_DEM 

#old lidar data, run all DEMs one by one all.
#./lidar2dem  -l ../data/Apollo15-LOLA/1E8E_21N28N/RDR_1E8E_21N28NPointPerRow_csv_table.csv ../data/LROC-DTMcomparison/ASU/Apollo15_comp/apollo15_dem.tif
#./lidar2dem  -l ../data/Apollo15-LOLA/1E8E_21N28N/RDR_1E8E_21N28NPointPerRow_csv_table.csv ../data/LROC-DTMcomparison/DLR/DLR_ap_dtm.tif
#./lidar2dem  -l ../data/Apollo15-LOLA/1E8E_21N28N/RDR_1E8E_21N28NPointPerRow_csv_table.csv ../data/LROC-DTMcomparison/UofA/Apollo15_DTM.tif
#./lidar2dem  -l ../data/Apollo15-LOLA/1E8E_21N28N/RDR_1E8E_21N28NPointPerRow_csv_table.csv ../data/LROC-DTMcomparison/USGS/USGS_M_111_06_16_DEM_150cm.tif

#old Lidar data run all DEMs at once
#./lidar2dem  -l ../data/Apollo15-LOLA/1E8E_21N28N/RDR_1E8E_21N28NPointPerRow_csv_table.csv $TIFF

#new lidar data, run all DEMs at once
#./lidar2dem  -l ../data/Apollo15-LOLA/RDR_3E4E_24N27NPointPerRow_csv_table.csv $TIFF

#run lidar2img multi processor - local
echo ../data/Apollo15-LOLA/*.csv |xargs -n 1 -I {} -P 2 ./lidar2dem  -l {} $A17_DEM

#run lidar2img multi processor - server
#echo /byss/moon/lola_tracks_over_apollo/*.csv |xargs -n 1 -I {} -P 2 ./lidar2dem  -l {} -m ../data/Apollo15-MAP-CUB -c ../data/Apollo15-CAM-CUB/*.lev1.cub
