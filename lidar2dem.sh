#!bin/sh
TIFF1=`ls ../data/LROC-DTMcomparison/ASU/Apollo15_comp/apollo15_dem.tif`
TIFF2=`ls ../data/LROC-DTMcomparison/DLR/DLR_ap_dtm.tif`
TIFF3=`ls ../data/LROC-DTMcomparison/UofA/Apollo15_DTM.tif`
TIFF4=`ls ../data/LROC-DTMcomparison/USGS/USGS_M_111_06_16_DEM_150cm.tif`
TIFF=$TIFF1" "$TIFF2" "$TIFF3" "$TIFF4
echo $TIFF

./lidar2dem  -l ../data/Apollo15-LOLA/1E8E_21N28N/RDR_1E8E_21N28NPointPerRow_csv_table.csv $TIFF


#keep this one
#./lidar2dem  -l ../data/Apollo15-LOLA/1E8E_21N28N/RDR_1E8E_21N28NPointPerRow_csv_table.csv ../data/LROC-DTMcomparison/ASU/Apollo15_comp/apollo15_dem.tif
#./lidar2dem  -l ../data/Apollo15-LOLA/1E8E_21N28N/RDR_1E8E_21N28NPointPerRow_csv_table.csv ../data/LROC-DTMcomparison/DLR/DLR_ap_dtm.tif
#./lidar2dem  -l ../data/Apollo15-LOLA/1E8E_21N28N/RDR_1E8E_21N28NPointPerRow_csv_table.csv ../data/LROC-DTMcomparison/UofA/Apollo15_DTM.tif
#./lidar2dem  -l ../data/Apollo15-LOLA/1E8E_21N28N/RDR_1E8E_21N28NPointPerRow_csv_table.csv ../data/LROC-DTMcomparison/USGS/USGS_M_111_06_16_DEM_150cm.tif


