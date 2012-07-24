#!/usr/bin/python
# takes a track CSV file and an image cub file as input, aligns tracks with the
# image using brute force search on multiple resolutions of the image pyramid

import sys
import os
import os.path
import shutil
import argparse
import tempfile
import math
import time
import glob
import numpy

from isis import *

CACHE_DIR = '../cache/'

#libraries needed for numpy conflict with those for ISIS. we put:
# alias python='OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH LD_LIBRARY_PATH="" python'
# in .bashrc to remove these libraries when loading python, here we restore them
os.environ['LD_LIBRARY_PATH'] = os.environ['OLD_LD_LIBRARY_PATH']

def bboxes_overlap(b1, b2):
	if b1[0] > b2[2]:
		return False
	if b2[0] > b1[2]:
		return False
	if b1[1] > b2[3]:
		return False
	if b2[1] > b1[3]:
		return False
	return True

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

def get_overlapping_tracks(tracks_dir, bbox):
	overlap = []
	for f in glob.glob(tracks_dir + '/*.csv'):
		b = get_tracks_bbox(f)
		if bboxes_overlap(bbox, b):
			overlap.append(f)
	return overlap


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

# scale the images, write a file containing the scale factors and file paths
def construct_image_pyramid(image_file):
	scale_reductions = [8, 4, 2]
	base = os.path.splitext(os.path.basename(image_file))[0]
	pyramid_file = os.path.join(CACHE_DIR, base + '.pyr')
	f = open(pyramid_file, 'w')
	for s in scale_reductions:
		fn = os.path.join(CACHE_DIR, base + '_' + str(s) + '.cub')
		if not os.path.exists(fn):
			reduced = reduce_image(image_file, s)
			if fn == None:
				f.close()
				os.remove(pyramid_file)
				return None
			shutil.move(reduced, fn)
		f.write('%d %s\n' % (s, fn))
	f.write('1 %s\n' % (image_file))
	f.close()
	return pyramid_file


# if alignments is specified, use those as starting track alignments
# if window is specified, do brute force search over a window
def align_image(tracks=None, tracks_list=None, image=None, image_pyramid=None, alignments=None, global_alignment=False, window=None, output_image=None, gcp_file=None):
	if tracks == None and tracks_list == None:
		return None
	if tracks_list != None:
		tracks_file_list = tempfile.NamedTemporaryFile(delete=False)
		for t in tracks_list:
			tracks_file_list.write(t + '\n')
		tracks_file_list.close()
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
	if tracks != None:
		command = '../bin/bruteforcealign -l %s -o %s %s -d tracks.dat ' % \
			(tracks, output_file, output_image)
	else:
		command = '../bin/bruteforcealign -t %s -o %s %s -d tracks.dat ' % \
			(tracks_file_list.name, output_file, output_image)
	if image != None:
		command += '-i ' + image + ' '
	if image_pyramid != None:
		command += '-p ' + image_pyramid + ' '
	if gcp_file != None:
		command += '-g ' + gcp_file + ' '
	if alignments != None:
		command += '--startMatrices ' + alignment_file + ' '
	if global_alignment:
		command += '--globalAlignment true '
	if window != None:
		command += '--transSearchWindow %d --transSearchStep %d --thetaSearchWindow %g --thetaSearchStep %g' % \
			(window[0], window[1], window[2], window[3])
	print command
	ret = os.system(command)
	if alignments != None:
		os.remove(alignment_file)
	if tracks_list != None:
		os.remove(tracks_file_list.name)
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
	
	pyramid_file = construct_image_pyramid(image_file)
	return align_image(tracks=tracks_file, image_pyramid=pyramid_file, global_alignment=True, output_image="reflectance.tif", gcp_file="gcp")
	
	#trimmed_image_file = os.path.join(CACHE_DIR, os.path.splitext(os.path.basename(image_file))[0] + '_%d_%d_%d_%d.cub' % (view_bbox[0], view_bbox[1], view_bbox[2], view_bbox[3]))
	#if not os.path.exists(trimmed_image_file):
	#	print "Cropping image file to tracks..."
		#trimmed_image_temp = crop_image(image_file, view_bbox)
		#shutil.move(trimmed_image_temp, trimmed_image_file)
	#if trimmed_image_file == None:
	#	return None

	#pyramid_file = construct_image_pyramid(trimmed_image_file)
	#return align_image(tracks=tracks_file, image_pyramid=pyramid_file, global_alignment=True, output_image="reflectance.tif", gcp_file="gcp")

def align_tracks(image_file, tracks_directory):
	print "Computing bounding box..."
	image_bbox = get_image_bbox(image_file)
	print "Determining overlapping tracks..."
	tracks = get_overlapping_tracks(tracks_directory, image_bbox)
	print "Constructing image pyramid..."
	pyramid_file = construct_image_pyramid(image_file)
	print "Aligning tracks..."
	return align_image(tracks_list=tracks, image_pyramid=pyramid_file, window=[20, 5, 0.0, math.pi/20], output_image="reflectance.tif")

parser = argparse.ArgumentParser(description='Align LOLA tracks to an image.')
parser.add_argument('-t', '--tracks', type=argparse.FileType('r'), required=False, metavar='tracks', help='A CSV file of LOLA tracks.')
parser.add_argument('-i', '--image', type=argparse.FileType('r'), required=True, metavar='image', help='A cub satelite image.')
parser.add_argument('-d', '--trackDirectory', type=str, required=False, metavar='trackDirectory', help='Directory to look for track files.')
res = parser.parse_args()

if res.tracks != None:
	tracks = res.tracks
	tracks_file = tracks.name
	tracks.close()
image = res.image
image_file = image.name
image.close()

start = time.time()
if res.tracks != None:
	brute_force_pyramid_align(image_file, tracks_file)
elif res.trackDirectory != None:
	align_tracks(image_file, res.trackDirectory)
else:
	print "Must specify either track file or track directory."
	sys.exit(1)
print 'Time taken: ' + str(time.time() - start)

