#!bin/sh
  
#refDEMFilename = "DTEEC_001513_1655_001777_1650_U01.tif"

#run the assembler in DEM mode wthout ICP
#./assembler  -m DEM -b ../../../msl/MSLData/Mars/MER_HIRISE/DTEEC_001513_1655_001777_1650_U01.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_no_icp_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-mod.tif
#mv ../../../msl/results/Mars/MER_HIRISE/assembled_dem.tif ../../../msl/results/Mars/MER_HIRISE/assembled_init_dem.tif 

# ./assembler  -m DEM_DRG -b ../../../msl/MSLData/Mars/MER_HIRISE/DEM_1m_ColumbiaHills-flat-crop-geo.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_no_icp_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-mod.tif
# mv ../../../msl/results/Mars/MER_HIRISE/assembled_drg.tif ../../../msl/results/Mars/MER_HIRISE/assembled_init_drg.tif 

#run the assembler in DEM mode with ICP
#./assembler  -m DEM -b ../../../msl/MSLData/Mars/MER_HIRISE/DTEEC_001513_1655_001777_1650_U01.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-mod.tif

# ./assembler  -m DEM_DRG -b ../../../msl/MSLData/Mars/MER_HIRISE/DEM_1m_ColumbiaHills-flat-crop-geo.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-mod.tif


  #generate shade-relief
  for demFile in `echo ../../../msl/results/Mars/MER_HIRISE/*dem.tif` 
  #for demFile in `echo ../MSLData/Mars/MER_HIRISE/DEM*.tif` 
  do
      shadeFile=${demFile/.tif/_shade.tif}
      colorShadeFile=${demFile/.tif/_clrshade.tif}
      echo "shadeFile = $shadeFile"
      echo "demFile = $demFile"
      #/Users/anefian/projects/visionworkbench/build/bin/hillshade -o $shadeFile -a 315 -s 0 --nodata-value 0 $demFile
      #/Users/anefian/projects/visionworkbench/build/bin/colormap  --lut-file LMMP_color_medium.lut -o $colorShadeFile -s $shadeFile --mars --legend --nodata-value 0 $demFile
      
      /opt/local/var/macports/software/gdal/1.8.0_0+expat+universal/opt/local/bin/gdaldem hillshade $demFile $shadeFile
      #/opt/local/var/macports/software/gdal/1.8.0_0+expat+universal/opt/local/bin/gdaldem color-relief $demFile LMMP_color_medium.lut $colorShadeFile
      #/opt/local/var/macports/software/gdal/1.8.0_0+expat+universal/opt/local/bin/gdaldem hillshade ../../../msl/MSLData/Mars/MER_HIRISE/DTEEC_001513_1655_001777_1650_U01.tif ../../../msl/MSLData/Mars/MER_HIRISE/DTEEC_001513_1655_001777_1650_U01_shade.tif
  done

  #run the assembler in DRG mode
#create tiles in the DRG tif images and prepare to read in VW 
gdal_translate -co TILED=YES -co BLOCKXSIZE=128 -co BLOCKYSIZE=128 ../../../msl/MSLData/Mars/MER_HIRISE/DT1EA_001513_1655_001777_1650_U01.tif ../../../msl/MSLData/Mars/MER_HIRISE/DT1EA_001513_1655_001777_1650_U01_block.tif
  #./assembler  -m DRG -b ../../../msl/MSLData/Mars/MER_HIRISE/PSP_001777_1650_1m_o-crop-geo.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Photo-mod.tif
 
