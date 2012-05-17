#!bin/sh

#test script to use openCV stereo matching for stereo image pairs in /irg/data/KRex/current/navigation/ 
#the disparity images and point images are saved ../results directory
#written by Ara

for leftImageFilename in /Users/anefian/projects/navigation/data/2012_0510_k10_marscape/*left*.pgm
do
 rightImageFilename=${leftImageFilename/-left-/-right-} 
 #echo $leftImageFilename
 #echo $rightImageFilename
 #./stereo -a $leftImageFilename -b $rightImageFilename
done


#point clouds to dem/drg/vrml
for stereoPair in `echo results/*`
do
    pc_txt=`echo $stereoPair/point*.txt`
    echo $pc_txt
    pc_tif=${pc_txt/.txt/-PC.tif}
    echo $pc_tif
    dem=${pc_txt/.txt/-DEM.tif} 
    echo $dem
    img_tif=${pc_txt/.txt/.tif}
    echo $img_tif 
    ../../../msl/scripts_msl/stereo $pc_txt $pc_tif
    #../../../StereoPipeline/src/asp/Tools/point2mesh $pc_tif
    ../../../StereoPipeline/src/asp/Tools/point2dem --mercator  -r mars $pc_tif -n --orthoimage $img_tif
    hillshade $dem
done

