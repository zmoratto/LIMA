#!bin/sh

#test script to use openCV stereo matching for stereo image pairs in /irg/data/KRex/current/navigation/ 
#the disparity images and point images are saved ../results directory
#written by Ara

for leftImageFilename in /Users/anefian/projects/navigation/data/2012_0510_k10_marscape/*left*.pgm
do
 rightImageFilename=${leftImageFilename/-left-/-right-} 
 echo $leftImageFilename
 echo $rightImageFilename
 ./stereo -a $leftImageFilename -b $rightImageFilename
done

