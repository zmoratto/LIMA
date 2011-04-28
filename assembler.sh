#!bin/sh
  
  #run the assembler in DEM_DRG mode wthout ICP
  ./assembler  -m DEM_DRG -b ../../../msl/MSLData/Mars/MER_HIRISE/DEM_1m_ColumbiaHills-flat-crop-geo.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_no_icp_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-mod.tif
  mv ../../../msl/results/Mars/MER_HIRISE/assembled_dem.tif ../../../msl/results/Mars/MER_HIRISE/assembled_init_dem.tif 
  mv ../../../msl/results/Mars/MER_HIRISE/assembled_drg.tif ../../../msl/results/Mars/MER_HIRISE/assembled_init_drg.tif 
  #run the assembler in DEM_DRG mode with ICP
  ./assembler  -m DEM_DRG -b ../../../msl/MSLData/Mars/MER_HIRISE/DEM_1m_ColumbiaHills-flat-crop-geo.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Height-mod.tif


  #generate shade-relief
  for demFile in `echo ../../../msl/results/Mars/MER_HIRISE/*dem.tif` 
  #for demFile in `echo ../MSLData/Mars/MER_HIRISE/DEM*.tif` 
  do
      shadeFile=${demFile/.tif/shade.tif}
      echo "shadeFile = $shadeFile"
      echo "demFile = $demFile"
   
      /opt/local/var/macports/software/gdal/1.8.0_0+expat+universal/opt/local/bin/gdaldem hillshade $demFile $shadeFile
  
  done

  #run the assembler in DRG mode
  #./assembler  -m DRG -b ../../../msl/MSLData/Mars/MER_HIRISE/PSP_001777_1650_1m_o-crop-geo.tif  -r ../../../msl/results/Mars/MER_HIRISE -c assembler_settings.txt ../../../msl/MSLData/Mars/MER_HIRISE/Photo-mod.tif
 