#!/usr/bin/python
# takes a track CSV file and an image cub file as input, aligns tracks with the
# image using brute force search on multiple resolutions of the image pyramid

import sys
import os
import subprocess
import argparse
import tempfile
import math
import time

import numpy

#libraries needed for numpy conflict with those for ISIS. we put:
# alias python='OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH LD_LIBRARY_PATH="" python'
# in .bashrc to remove these libraries when loading python, here we restore them
os.environ['LD_LIBRARY_PATH'] = os.environ['OLD_LD_LIBRARY_PATH']

if not os.environ.has_key('ISISROOT'):
	print >> sys.stderr, 'The environment variable ISISROOT is undefined.'
	sys.exit(1)
ISISROOT = os.environ['ISISROOT']
ISISBIN = ISISROOT + '/bin'

# returns [minlon, minlat, maxlon, maxlat] of box enclosing LOLA tracks
def get_tracks_bbox(tracks_file):
	f = open(tracks_file, 'r')
	if f == None:
		print >> sys.stderr, 'File %s could not be opened.' % (tracks_file)
		return None
	minlon = 180; maxlon = -180; minlat = 90; maxlat = -90
	f.readline() # skip the first line of headers
	for l in f:
		line = l.split(',')
		try:
			lon = float(line[1])
			lat = float(line[2])
		except:
			print >> sys.stderr, 'Could not parse %s and %s as longitude and latitude.' % (line[1], line[2])
		if lon < minlon:
			minlon = lon
		if lon > maxlon:
			maxlon = lon
		if lat < minlat:
			minlat = lat
		if lat > maxlat:
			maxlat = lat
	f.close()
	return [minlon, minlat, maxlon, maxlat]

# get longitude / latitude borders of image
def get_image_bbox(image_file):
	proc = subprocess.Popen([ISISBIN + '/camrange', 'from=' + image_file], stdout=subprocess.PIPE)
	proc.wait()
	os.remove('print.prt')
	if proc.returncode != 0:
		print >> sys.stderr, "Error in camrange from=%s, aborting." % (image_file)
		return None
	for l in proc.stdout:
		if l.startswith('Group = UniversalGroundRange'):
			break
		if l == '':
			print >> sys.stderr, "UniversalGroundRange not found in camrange output."
	minlon = None; maxlon = None; minlat = None; maxlat = None
	for l in proc.stdout:
		if l.startswith('End_Group'):
			break
		s = l.split()
		try:
			if s[0] == 'MaximumLatitude':
				maxlat = float(s[2])
			if s[0] == 'MinimumLatitude':
				minlat = float(s[2])
		except ValueError:
			print >> sys.stderr, 'Could not convert %s to float in camrange output.' % (s[2])
			return None
	for l in proc.stdout:
		if l.startswith('Group = PositiveEast180'):
			break
		if l == '':
			print >> sys.stderr, "PositiveEast180 not found in camrange output."
	for l in proc.stdout:
		if l.startswith('End_Group'):
			break
		s = l.split()
		try:
			if s[0] == 'MaximumLongitude':
				maxlon = float(s[2])
			if s[0] == 'MinimumLongitude':
				minlon = float(s[2])
		except ValueError:
			print >> sys.stderr, 'Could not convert %s to float in camrange output.' % (s[2])
			return None
	
	if minlon == None or maxlon == None or minlat == None or maxlat == None:
		print >> sys.stderr, "Not all boundaries found in camrange output."
		return None

	return [minlon, minlat, maxlon, maxlat]

def get_image_dimensions(image_file):
	tfile = tempfile.NamedTemporaryFile(delete=False)
	out_file = tfile.name
	tfile.close()
	ret = os.system(ISISBIN + '/caminfo from=%s to=%s format=flat' % (image_file, out_file))
	os.remove('print.prt')
	if ret != 0:
		return (None, None)
	f = open(out_file, 'r')
	f.readline() # header
	s = f.readline().split(',')
	os.remove(out_file)
	return (int(s[4]), int(s[5]))

def get_point_projection(image_file, lon, lat):
	tfile = tempfile.NamedTemporaryFile(delete=False)
	out_file = tfile.name
	tfile.close()
	if lon < 0.0:
		lon = lon + 360.0
	ret = os.system(ISISBIN + '/campt from=%s to=%s format=flat type=ground latitude=%g longitude=%g > /dev/null' % (image_file, out_file, lat, lon))
	os.remove('print.prt')
	if ret != 0:
		return (None, None)
	f = open(out_file, 'r')
	s = f.readline().split(',')
	os.remove(out_file)
	return (int(float(s[1])), int(float(s[2])))

# create bounds for a trimmed image including as many of the LOLA points as possible,
# with some wiggle room around the edges for alignment
def get_bounding_box(tracks_file, image_file):
	lolabbox = get_tracks_bbox(tracks_file)
	imagebbox = get_image_bbox(image_file)
	if lolabbox == None or imagebbox == None:
		return None
	(width, height) = get_image_dimensions(image_file)
	if width == None or height == None:
		return None

	bbout = [0, 0, 0, 0]
	bbout[0] = max(lolabbox[0], imagebbox[0])
	bbout[1] = max(lolabbox[1], imagebbox[1])
	bbout[2] = min(lolabbox[2], imagebbox[2])
	bbout[3] = min(lolabbox[3], imagebbox[3])

	if bbout[0] >= bbout[2] or bbout[1] >= bbout[3]:
		print >> sys.stderr, "No overlap between LOLA points and the image."
		return None

	horiz_expand = (bbout[2] - bbout[0]) * 0.2
	bbout[0] = max(bbout[0] - horiz_expand, imagebbox[0])
	bbout[2] = min(bbout[2] + horiz_expand, imagebbox[2])
	
	vert_expand = (bbout[3] - bbout[1]) * 0.2
	bbout[1] = max(bbout[1] - horiz_expand, imagebbox[1])
	bbout[3] = min(bbout[3] + horiz_expand, imagebbox[3])

	(x1, y1) = get_point_projection(image_file, bbout[0], bbout[1])
	(x2, y2) = get_point_projection(image_file, bbout[2], bbout[3])
	x1 = min(max(x1, 1), width)
	x2 = min(max(x2, 1), width)
	y1 = min(max(y1, 1), height)
	y2 = min(max(y2, 1), height)
	return [min(x1,x2), min(y1,y2), max(x1,x2), max(y1,y2)]

def crop_image(image_file, view_bbox):
	tfile = tempfile.NamedTemporaryFile(delete=False, suffix='.cub')
	trimmed_file = tfile.name
	tfile.close()
	command = ISISBIN + '/crop from=%s to=%s sample=%d nsamples=%d line=%d nlines=%d > /dev/null' % \
		(image_file, trimmed_file, view_bbox[0] + 1, view_bbox[2] - view_bbox[0], view_bbox[1] + 1, view_bbox[3] - view_bbox[1])
	ret = os.system(command)
	os.remove('print.prt')
	if ret != 0:
		print >> sys.stderr, "Failed to crop image."
		return None
	return trimmed_file

def reduce_image(in_file, scale_factor):
	tfile = tempfile.NamedTemporaryFile(delete=False, suffix='.cub')
	reduced_file = tfile.name
	tfile.close()
	command = ISISBIN + '/reduce from=%s to=%s sscale=%g lscale=%g > /dev/null' % \
		(in_file, reduced_file, scale_factor, scale_factor)
	ret = os.system(command)
	os.remove('print.prt')
	if ret != 0:
		print >> sys.stderr, "Failed to downscale image."
		return None
	return reduced_file

# if alignments is specified, use those as starting track alignments
# if window is specified, do brute force search over a window
def align_image(tracks, image, alignments=None, window=None, output_image=None):
	if alignments != None:
		tfile = tempfile.NamedTemporaryFile(delete=False)
		alignment_file = tfile.name
		tfile.close()
		f = open(alignment_file, "w")
		for m in alignments:
			f.write('%g %g %g %g %g %g\n' % (m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2]))
		f.close()
	
	if output_image != None:
		output_image = '--outputImage ' + output_image
	else:
		output_image = ""

	tfile = tempfile.NamedTemporaryFile(delete=False)
	output_file = tfile.name
	tfile.close()
	command = '../bin/bruteforcealign -l %s -i %s -o %s %s ' % \
		(tracks, image, output_file, output_image)
	if alignments != None:
		command += '--startMatrices ' + alignment_file + ' '
	if window != None:
		command += '--transSearchWindow %d --transSearchStep %d --thetaSearchWindow %g --thetaSearchStep %g' % \
			(window[0], window[1], window[2], window[3])
	ret = os.system(command)
	if alignments != None:
		os.remove(alignment_file)
	if ret != 0:
		print >> sys.stderr, "Failed to align image."
		return None

	out = []
	f = open(output_file, 'r')
	for l in f:
		if l == "":
			continue
		temp = numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
		s = map(float, l.split())
		temp[0] = s[0:3]
		temp[1] = s[3:6]
		out.append(temp)
	f.close()
	os.remove(output_file)
	return out

def brute_force_pyramid_align(image_file, tracks_file):
	print "Computing bounding box..."
	view_bbox = get_bounding_box(tracks_file, image_file)
	if view_bbox == None:
		return None
	
	print "Cropping image file to tracks..."
	trimmed_image_file = crop_image(image_file, view_bbox)
	if trimmed_image_file == None:
		return None

	# first do a big search over the smallest image
	old_factor = 8
	print 'Downscaling image by scale factor %d...' % (old_factor)
	reduced = reduce_image(trimmed_image_file, old_factor)
	print 'Aligning image...'
	alignments = align_image(tracks_file, reduced, window=[20, 5, math.pi/20, math.pi/20], output_image="reflectance.tif")
	
	# now do smaller searches over more detailed images, using previous alignment as starting point
	scale_reductions = [4, 2, 1]
	search_window = [2, 1, math.pi / 40, math.pi / 80]
	for factor in scale_reductions:
		# scale matrix for new resolution
		for i in range(len(alignments)):
			r = float(old_factor) / factor
			t1 = numpy.array([[r, 0, 0], [0, r, 0], [0, 0, 1]])
			t2 = numpy.array([[1.0/r, 0, 0], [0, 1.0/r, 0], [0, 0, 1]])
			alignments[i] = t1.dot(alignments[i]).dot(t2)
		old_factor = factor
		if factor == 1:
			reduced = trimmed_image_file
		else:
			print "Reducing image by scale factor %d." % (factor)
			reduced = reduce_image(trimmed_image_file, factor)
		print 'Aligning image...'
		alignments = align_image(tracks_file, reduced, alignments=alignments, output_image='reflectance.tif')
		# also use finer theta search window
		search_window[2] = search_window[2] / 2
		search_window[3] = search_window[3] / 2
	return alignments

parser = argparse.ArgumentParser(description='Align LOLA tracks to an image.')
parser.add_argument('-t', '--tracks', type=argparse.FileType('r'), required=True, metavar='tracks', help='A CSV file of LOLA tracks.')
parser.add_argument('-i', '--image', type=argparse.FileType('r'), required=True, metavar='image', help='A cub satelite image.')
res = parser.parse_args()

tracks = res.tracks
tracks_file = tracks.name
tracks.close()
image = res.image
image_file = image.name
image.close()

start = time.time()
brute_force_pyramid_align(image_file, tracks_file)
print 'Time taken: ' + str(time.time() - start)

