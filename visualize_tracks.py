#!/usr/bin/python
# Takes a track data file output by save_track_data which contains 
# alternating rows of distance data, distances between the points,
# the synthetic image, and the color of the points on the actual image (after adjustment).
# Graphs the depths of the points and shows the synthetic and original image
# colors using matplotlib.

import sys
from pylab import rand

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse

if len(sys.argv) < 2:
	print 'Please pass the track data file as an input.'
	sys.exit(-1)

f = open(sys.argv[1], 'r')

while True:
	s1 = f.readline()
	if s1 == "":
		break
	s2 = f.readline()
	if s2 == "":
		break
	s3 = f.readline()
	if s3 == "":
		break
	s4 = f.readline()
	if s4 == "":
		break

	depths = map(float, s1.split())
	x = map(float, s2.split())
	reflectance = map(float, s3.split())
	image = map(float, s4.split())
	if len(reflectance) != len(image):
		print >> sys.stderr, "Number of reflectance and image points should be equal."
	small = reduce(min, depths, depths[0])
	depths = map(lambda f: f - small, depths)
	start = x[0]
	x = map(lambda f: f - start, x)

	fig = plt.figure()
	g = fig.add_subplot(3, 1, 1)
	g.plot(x, depths)
	g.set_xlim(x[0], x[len(x)-1])

	ref = fig.add_subplot(3, 1, 2 )
	for i in xrange(len(reflectance)):
		e = Rectangle(xy=(i,0), width=1, height=1, edgecolor='none')
		r = min(reflectance[i], 1)
		e.set_facecolor([r, r, r])
		ref.add_artist(e)

	ref.set_frame_on(False)
	ref.get_xaxis().set_visible(False)
	ref.get_yaxis().set_visible(False)
	ref.set_xlim(0, len(reflectance))
	ref.set_ylim(0, 1)
	
	ref = fig.add_subplot(3, 1, 3 )
	for i in xrange(len(image)):
		e = Rectangle(xy=(i,0), width=1, height=1, edgecolor='none')
		e.set_facecolor([image[i], image[i], image[i]])
		ref.add_artist(e)

	ref.set_frame_on(False)
	ref.get_xaxis().set_visible(False)
	ref.get_yaxis().set_visible(False)
	ref.set_xlim(0, len(image))
	ref.set_ylim(0, 1)

	plt.show()

