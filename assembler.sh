#!bin/sh
  
#run the assembler in DEM mode with ICP
#can be used with assembler_no_ICP_settings.txt or assembler_settings.txt
#./assembler  -m DEM -b ../../../msl/MSLData/Mars/MER_HIRISE/DTEEC_001513_1655_001777_1650_U01.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-Sol-855.tif
#./assembler  -m DEM -b ../../../msl/MSLData/Mars/MER_HIRISE/DTEEC_001513_1655_001777_1650_U01.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_no_icp_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-Sol-855.tif
  #generate shade-relief
  for demFile in `echo ../../../msl/results/Mars/MER_HIRISE/*dem.tif` 
 
  do
      shadeFile=${demFile/.tif/_shade.tif}
      colorShadeFile=${demFile/.tif/_clrshade.tif}
      echo "shadeFile = $shadeFile"
      echo "demFile = $demFile"
   
      #this is good - START
      #gdaldem hillshade $demFile $shadeFile
      #this is good - END


  done

#run the assembler in DRG mode
./assembler  -m DRG -b ../../../msl/MSLData/Mars/MER_HIRISE/DT1EA_001513_1655_001777_1650_U01.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Photo-Sol-855.tif
