#!bin/sh

for tileIndex in 76_73 76_74 77_73 77_74
do 
  inputFilename=../../../msl/results/Mars/MER_HIRISE/assembled_$tileIndex"_dem.tif"  
  resDirname=../../../msl/results/Mars/MER_HIRISE/$tileIndex"_DEM"
  echo "inputFilename=$inputFilename"
  echo "resDirname=$resDirname"

  #clean the output directory
  rm $resDirname/*shade.tif 
  rm $resDirname/*.vrt
  rm $resDirname/res_*.tif

  #run the georef_resample application
  ./geotif_resample -i $inputFilename -o $resDirname -c geotif_resample_DEM_settings.txt

  #check to see if the reconstructed DEM look reasonably
  #create conf mosaic for display
  gdalbuildvrt $resDirname/dem_1_mosaic.vrt $resDirname/res_1_*.tif
  gdal_translate -of GTiff -outsize 100% 100% $resDirname/dem_1_mosaic.vrt $resDirname/dem_mosaic_1.tif
  #generate shade-relief
  gdaldem hillshade $resDirname/dem_mosaic_1.tif $resDirname/dem_mosaic_1_shade.tif
   
done
