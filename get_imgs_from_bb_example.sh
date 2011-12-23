#!bin/sh


Apollo_DEM=`ls -d -1 /Users/anefian/projects/LMMP/LMMP_delivery/2011_1102_Apollo_DEM/*.tif`
#run for all geotiifs at once
./get_imgs_from_bb -n 30 -s -30 -e -50 -w 50 $Apollo_DEM

