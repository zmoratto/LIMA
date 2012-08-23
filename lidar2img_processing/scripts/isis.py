# Call ISIS utilities from python.

import sys
import os
import subprocess
import tempfile

if not os.environ.has_key('ISISROOT'):
	print >> sys.stderr, 'The environment variable ISISROOT is undefined.'
	sys.exit(1)
ISISROOT = os.environ['ISISROOT']
ISISBIN = ISISROOT + '/bin'

def get_image_coords(image_file, lat, lon):
	proc = subprocess.Popen([ISISBIN + '/campt', 'from=' + image_file, 'type=ground', 'latitude=' + str(lat), 'longitude=' + str(lon)], stdout=subprocess.PIPE)
	proc.wait()
	os.remove('print.prt')
	if proc.returncode != 0:
		print >> sys.stderr, "Error in campt, aborting."
		return None
	x = None
	y = None
	for l in proc.stdout:
		s = l.split()
		try:
			if s[0] == 'Sample':
				x = float(s[2])
				if y != None:
					break
			if s[0] == 'Line':
				y = float(s[2])
				if x != None:
					break
		except ValueError:
			print >> sys.stderr, 'Could not convert %s to float in campt output.' % (s[2])
			return None
	if x == None or y == None:
		print >> sys.stderr, 'Did not find latitude and longitude in campt output.'
		return None
	return (x, y)

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

def normalize_image(image_file):
	tfile = tempfile.NamedTemporaryFile(delete=False, suffix='.cub')
	normalized_file = tfile.name
	tfile.close()
	command = ISISBIN + '/bandnorm from=%s to=%s average=cube > /dev/null' % (image_file, normalized_file)
	ret = os.system(command)
	os.remove('print.prt')
	if ret != 0:
		print >> sys.stderr, "Failed to normalize image."
		return None
	return normalized_file

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

