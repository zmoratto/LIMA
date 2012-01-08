#!bin/sh
  
#refDEMFilename = "DTEEC_001513_1655_001777_1650_U01.tif"

#run the assembler in DEM mode wthout ICP
#./assembler  -m DEM -b ../../../msl/MSLData/Mars/MER_HIRISE/DTEEC_001513_1655_001777_1650_U01.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_no_icp_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-mod.tif
#./assembler  -m DEM -b ../../../msl/MSLData/Mars/MER_HIRISE/DTEEC_001513_1655_001777_1650_U01.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_no_icp_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-Sol-855.tif
#mv ../../../msl/results/Mars/MER_HIRISE/assembled_dem.tif ../../../msl/results/Mars/MER_HIRISE/assembled_init_dem.tif 

#./assembler  -m DEM_DRG -b ../../../msl/MSLData/Mars/MER_HIRISE/DEM_1m_ColumbiaHills-flat-crop-geo.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_no_icp_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-mod.tif
# mv ../../../msl/results/Mars/MER_HIRISE/assembled_drg.tif ../../../msl/results/Mars/MER_HIRISE/assembled_init_drg.tif 

#run the assembler in DEM mode with ICP
#this is good - START
#can be used with assembler_no_ICP_settings.txt or assembler_settings.txt
./assembler  -m DEM -b ../../../msl/MSLData/Mars/MER_HIRISE/DTEEC_001513_1655_001777_1650_U01.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_no_ICP_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-Sol-855.tif
#this is good END

#./assembler  -m DEM_DRG -b ../../../msl/MSLData/Mars/MER_HIRISE/DEM_1m_ColumbiaHills-flat-crop-geo.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-mod.tif

#./assembler  -m DEM_DRG -b  ../../../msl/MSLData/Mars/MER_HIRISE/DTEEC_001513_1655_001777_1650_U01.tif -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-Sol-855.tif

  #generate shade-relief
  for demFile in `echo ../../../msl/results/Mars/MER_HIRISE/*dem.tif` 
  #for demFile in `echo ../MSLData/Mars/MER_HIRISE/DEM*.tif` 
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
#create tiles in the DRG tif images and prepare to read in VW. not working since it is not supported by the installed GDAL and libtiff versions. 
#gdal_translate -co BIGTIFF=NO -co TILED=YES -co BLOCKXSIZE=128 -co BLOCKYSIZE=128 ../../../msl/MSLData/Mars/MER_HIRISE/DT1EA_001513_1655_001777_1650_U01.tif ../../../msl/MSLData/Mars/MER_HIRISE/DT1EA_001513_1655_001777_1650_U01_block.tif

#running the assembler with HiRISE images already cropped 
#this will be replaced when we will be able to read HiRISE DRGs
#./assembler  -m DRG -b ../../../msl/MSLData/Mars/MER_HIRISE/PSP_001777_1650_1m_o-crop-geo.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Photo-Sol-855.tif


./assembler  -m DRG -b ../../../msl/MSLData/Mars/MER_HIRISE/res.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Photo-Sol-855.tif
#./assembler  -m DRG -b ../../../msl/MSLData/Mars/MER_HIRISE/res.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Photo-mod.tif