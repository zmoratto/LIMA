#!bin/sh

for tileIndex in 76_73 76_74 77_73 77_74
do 


  inputFilename=../../../msl/results/Mars/MER_HIRISE/assembled_$tileIndex"_drg.tif"
  upsampledFilename=../../../msl/results/Mars/MER_HIRISE/upsampled_$tileIndex"_drg.tif"
    
  resDirname=../../../msl/results/Mars/MER_HIRISE/$tileIndex"_DRG"
  echo "inputFilename=$inputFilename"
  echo "resDirname=$resDirname"
   
  gdal_translate -outsize 400% 400% $inputFilename $upsampledFilename

  #clean the output directory
  rm $resDirname/*.vrt
  rm $resDirname/res_*.tif

  #run the georef_resample application
  ./geotif_resample -i $upsampledFilename -o $resDirname -c geotif_resample_DRG_settings.txt

  #check to see if the reconstructed DEM look reasonably
  #create conf mosaic for display
  gdalbuildvrt $resDirname/drg_1_mosaic.vrt $resDirname/res_1_*.tif
  gdal_translate -of GTiff -outsize 100% 100% $resDirname/drg_1_mosaic.vrt $resDirname/drg_mosaic_1.tif
   
done
