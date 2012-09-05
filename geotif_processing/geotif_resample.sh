#!bin/sh


#clean the output directory
#rm ../../../msl/results/Mars/MER_HIRISE/res_*shade.tif
#rm ../../../msl/results/Mars/MER_HIRISE/res_*.tif

#run the georef_resample application
./geotif_resample -i ../../../msl/MSLData/Mars/MER_HIRISE/DTEEC_001513_1655_001777_1650_U01.tif -o ../../../msl/results/Mars/MER_HIRISE -c muja.txt

#sanity check to see if the reconstructed DEM look reasonably
#create conf mosaic for display
#gdalbuildvrt ../../../msl/results/Mars/MER_HIRISE/dem_0_mosaic.vrt ../../../msl/results/Mars/MER_HIRISE/res_0_*.tif
#gdal_translate -of GTiff -outsize 100% 100% ../../../msl/results/Mars/MER_HIRISE/dem_0_mosaic.vrt ../../../msl/results/Mars/MER_HIRISE/dem_mosaic_0_10pct.tif

gdalbuildvrt ../../../msl/results/Mars/MER_HIRISE/dem_1_mosaic.vrt ../../../msl/results/Mars/MER_HIRISE/res_1_*.tif
gdal_translate -of GTiff -outsize 100% 100% ../../../msl/results/Mars/MER_HIRISE/dem_1_mosaic.vrt ../../../msl/results/Mars/MER_HIRISE/dem_mosaic_1_10pct.tif

gdalbuildvrt ../../../msl/results/Mars/MER_HIRISE/dem_2_mosaic.vrt ../../../msl/results/Mars/MER_HIRISE/res_2_*.tif
gdal_translate -of GTiff -outsize 100% 100% ../../../msl/results/Mars/MER_HIRISE/dem_2_mosaic.vrt ../../../msl/results/Mars/MER_HIRISE/dem_mosaic_2_10pct.tif

gdalbuildvrt ../../../msl/results/Mars/MER_HIRISE/dem_3_mosaic.vrt ../../../msl/results/Mars/MER_HIRISE/res_3_*.tif
gdal_translate -of GTiff -outsize 100% 100% ../../../msl/results/Mars/MER_HIRISE/dem_3_mosaic.vrt ../../../msl/results/Mars/MER_HIRISE/dem_mosaic_3_10pct.tif

 #generate shade-relief
  for demFile in `echo ../../../msl/results/Mars/MER_HIRISE/dem_mosaic*.tif` 
 
  do
      shadeFile=${demFile/.tif/_shade.tif}
      echo "shadeFile = $shadeFile"
      
      #this is good - START
      gdaldem hillshade $demFile $shadeFile
      #this is good - END
  done

