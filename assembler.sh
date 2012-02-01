#!bin/sh
  
#run the assembler in DEM mode with ICP
#can be used with assembler_no_ICP_settings.txt or assembler_settings.txt
#./assembler  -m DEM -b ../../../msl/MSLData/Mars/MER_HIRISE/DTEEC_001513_1655_001777_1650_U01.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-Sol-855.tif

  #generate shade-relief
  for demFile in `echo ../../../msl/results/Mars/MER_HIRISE/*dem.tif` 
 
  do
      shadeFile=${demFile/.tif/_shade.tif}
      colorShadeFile=${demFile/.tif/_clrshade.tif}
      echo "shadeFile = $shadeFile"
      echo "demFile = $demFile"
 
      #VW hillshade works well for the final DEM but does not show anything on the initDEM. Needs a fix here!      
      #hillshade -o $shadeFile -a 315 -s 0 --nodata-value 0 $demFile

      #VW colormap shows basicaly the same color on both init and final assembled DEM. (bad colormap, or a bug?)
      #colormap  --lut-file LMMP_color_medium.lut -o $colorShadeFile -s $shadeFile --mars --legend --nodata-value 0 $demFile
      
      #this is good - START
      gdaldem hillshade $demFile $shadeFile
      #this is good - END


  done

#run the assembler in DRG mode

./assembler  -m DRG -b ../../../msl/MSLData/Mars/MER_HIRISE/DT1EA_001513_1655_001777_1650_U01.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Photo-Sol-855.tif
