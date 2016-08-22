#!/usr/bin/env python
# =======================================================================================
# Copyright (C) 2016 Kondrak Matthias
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# =======================================================================================
# Title:		sigma.py
# Author:		Kondrak Matthias
# Contact:		matthias.kondrak@uibk.ac.at
# Address:		University of Innsbruck, Institute of Astro and Particle Physics
#				Technikerstrasse 25, 08/22
#				A-6020 Innsbruck
# Date:			22/08/2016
# Version:		3.0
# 
# Description:	This script uses the "scipy.stats" library and the sigmaclip
#				function to clip the data. Input is the data file (ASCII) 
#				and a selectable clipping factor for a "threshold" ->
#				mean - std * factor < data < mean + std * factor
#				Clipping is applied to each orbit seperatly, if there are more
#				than 2 datapoints within such an orbit.
#				This script overwrites the input file!
#
# Parameters:	two;
#				-----------------------------------------------------------------
#				$one = (string) Input file
#				$two = (float) Factor of clipping
#				-----------------------------------------------------------------
#
# Input:		ASCII-file with comma-separator ',' and factor for clipping
#
# Output:   	A string to stdout with the entries of the first column provided
#				by the input (unique identifier!) and a 'nan'-entry for deleted/removed
#				values.
#
# Exit:   		0 - everything went fine
#				1 - not enough input-parameters
#				2 - file does not exist
#
# =======================================================================================
# Import libraries
import sys
import os.path
import numpy as np
from scipy.stats import sigmaclip

# =======================================================================================
# Checks if there are enough input parameters; else exits
if (len(sys.argv) <= 2):
	print "Usage", str(sys.argv[0]), "<FILE> <SIGMA-VALUE>"
	sys.exit(1)

# Checks if file exists
if not (os.path.isfile(sys.argv[1])):
	print "File does not exist!"
	sys.exit(2)

# Checks if the sigma value is too small (which leads to an error)
if float(sys.argv[2]) < 1.3:
	# Sets a default value of 1.3
	clip = 1.3
else:
	clip = float(sys.argv[2])

# Defines the array for clipping
hjd = list()
mag = list()
tmp = list()
out = list()
orbit = list()

# ---------------------------------------------------------------------------------------
# Opens the file - reading
f = open(sys.argv[1], 'r')

# Reads line by line of the file
for line in f:
	# Gets the line up to the line-seperator
	line = line.strip()
	# Gets the column from the current line, seperated by ','
	column = line.split(',')
	# Adds the value to the list
	hjd.append(column[0])
	mag.append(float(column[1]))
	orbit.append(float(column[5]))
	
# Close file
f.close()

# Gets the first orbital number
n = orbit[0]
# Loops through the whole dataset
for i in range(0, len(mag)):
	# Appends the next value to list
	tmp.append(mag[i])
	# Checks if the orbital number has changed
	if not orbit[i] == n:
		# Checks if there are enough values for clipping procedure
		if len(tmp) > 2:
			# Clipps and gets back an array with clipped values
			c, down, up = sigmaclip(tmp, clip, clip)
			# Sets looping variable to 0
			k = 0
			# Loops through the current orbital values
			for j in range(0, len(tmp)):
				# Checks if the value was not clipped
				if tmp[j] == c[k]:
					# Appends the non-clipped value to output
					out.append(tmp[j])
					# Increases variable, if it does not exceed
					# length of list
					if not (k + 1) == len(c):
						k = k + 1
				# If value was clipped - appends a NaN (not a number) entry
				else:
					out.append(np.nan)
		# If there are to few values for clipping
		else:
			# Just appends NAN entries
			for j in range(0, len(tmp)):
				out.append(np.nan)
		# Flushes/resets the list
		del tmp[:]
		# Sets new orbital number
		n = orbit[i + 1]
	
# Same procedure for the last orbit in file - because loop will not handle it
if len(tmp) > 2:
	c, down, up = sigmaclip(tmp, clip, clip)
	k = 0
	for j in range(0, len(tmp)):
		if tmp[j] == c[k]:
			out.append(tmp[j])
			if not (k + 1) == len(c):
				k = k + 1
		else:
			out.append(np.nan)
else:
	for j in range(0, len(tmp)):
		out.append(np.nan)
		
# ---------------------------------------------------------------------------------------
# Writes the modified data to std-out (only first column)
sys.stdout.write(str(out))

# Exit
sys.exit(0)
# ---------------------------------------------------------------------------------------
