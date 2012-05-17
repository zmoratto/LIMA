#the calibration code for retangular patterns with -w 9 horizontal and -h 7 vertical inner squares.
#images are stored in subList.txt
#scale of rectangles on the checkerboard are 4 inches=101.6mm (-s 101.6)

#these data sets are useless, the checkerboard does not cover the image surface entirely - START
#./stereo_calibration -w 9 -h 7 2012_0405_marscape_list_short.txt -s 101.6 --show-rectified --show-corners
#these data sets are useless, the checkerboard does not cover the image surface entirely - - END

dirname=2012_0423_marscape

#create the list of stereo calibration images - START
leftFilename=$dirname'_left.txt'
rightFilename=$dirname'_right.txt'

ls $dirname/*left-*.pgm>$leftFilename
ls $dirname/*right-*.pgm>$rightFilename

paste -d\\n $leftFilename $rightFilename>$dirname'_list.txt'

rm $leftFilename
rm $rightFilename
#create the list of stereo calibration images - START

#good data sets - START
#./stereo_calibration -w 9 -h 7 2012_0213_marscape_list.txt -s 101.6 --show-rectified --show-corners
./stereo_calibration -w 9 -h 7 $dirname_list.txt -s 101.6 --show-rectified --show-corners
#good data sets - END