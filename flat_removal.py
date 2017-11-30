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
# Title:		flat_removal.py
# Author:		Kondrak Matthias
# Contact:		matthias.kondrak@uibk.ac.at
# Address:		University of Innsbruck, Institute of Astro and Particle Physics
#				Technikerstrasse 25, 08/22
#				A-6020 Innsbruck
# Date:			22/08/2016
# Version:		3.0
# 
# Description:	This script is a graphical assistent for a 'flat_removal'. 
#				The magnitude vs. HJD is plotted. The user can draw two borders
#				with the mouse to exclude datapoints with a large scattering.
#				A string with the magnitude of the two borders is returned, which
#				can be removed in a next step.
#
# Parameters:	one;
#				-----------------------------------------------------------------
#				$one = (string) Input file
#				-----------------------------------------------------------------
#
# Input:		ASCII-file with comma-seperator ',' and a header in the first
#				line. Mouse-input by drawing inside the plot.
#				File has a minimum of 2 columns (HJD, mag)
#
# Output:   	String to STDOUT with the y-coordinates of the rectangles
#				Y0 Y1
#				Separated by space ' ' and sorted.
#				Default value is -50 50
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
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# =======================================================================================
# Defines the class for drawing rectangulars and selecting the aperature correction
class draw_rect(object):
    # Initialises all variables
    def __init__(self, axis1):
        # Defines the axis used
        self.ax1 = axis1
        # Defines the default y-coordinate of the rectangle
        self.y0 = -50
        self.y1 = 50
        # Defines variables for drawing enabled/disabled and default output
        self.draw = 0
        self.output = "-50 50"
        # Defines a rectangle for the plot (mag/hjd)
        # Sets the rectangular width to 10000 - should be enough
        self.rect1 = Rectangle((0, 0), 10000, 0, alpha=0.3)
        # Adds the rectangle to the corresponding axis
        self.ax1.add_patch(self.rect1)

    # Defines what should happen on pressing the mousebutton
    def on_press(self, event):
        # At first gets the toolbar-button-pressed-info and checks if
        # no toolbar-button is pressed (zooming or shifting/moving)
        # So zooming/shifting can be done without drawing a new rectangle
        toolbar = plt.get_current_fig_manager().toolbar
        if toolbar.mode == '':
            # If the mouse is INSIDE the plot
            if event.inaxes is self.ax1:
                # Enables drawing
                self.draw = 1
                # Gets the y-coordinate, where the mouse it located
                self.y0 = event.ydata

    # On mouse button release
    def on_release(self, event):
        # Only if it is INSIDE the first plot
        if event.inaxes is self.ax1:
            # Sorts the x- and y-values and sets the output
            if self.y0 < self.y1:
                self.output = str(self.y0) + " " + str(self.y1)
            else:
                self.output = str(self.y1) + " " + str(self.y0)

            # Redraws the rectangle (if last drawing went wrong)
            self.rect1.set_height(self.y1 - self.y0)
            self.rect1.set_xy((0, self.y0))
            self.ax1.figure.canvas.draw()

        # Sets drawing to 0 - disable drawing
        self.draw = 0

    # On mouse motion
    def on_motion(self, event):
        # If drawing is enabled
        if (self.draw == 1):
            # Checks if inside first plot
            if event.inaxes is self.ax1:
                # Gets the new (because moving) coordinates, calculates
                # the rectangle and draws it during moving, - this can take
                # a few seconds due to calulating/redrawing - depending on
                # size of data
                self.y1 = event.ydata
                self.rect1.set_height(self.y1 - self.y0)
                self.rect1.set_xy((0, self.y0))
                self.ax1.figure.canvas.draw()
            else:
                # If not inside the plot or leaving it - disables drawing
                self.draw = 0

    # Returns the output
    def ret(self):
        return self.output


# =======================================================================================
# Checks if there are enough input parameters; else exits
if (len(sys.argv) < 1):
    print
    "Usage", str(sys.argv[0]), "<FILE>"
    sys.exit(1)

# Checks if file exists
if not (os.path.isfile(sys.argv[1])):
    print
    "File does not exist!"
    sys.exit(2)

# Defines variables
hjd = list()
mag = list()

# ---------------------------------------------------------------------------------------
# Opens the file
f = open(sys.argv[1], 'r')

# Gets the header and splits it up, so that every entry of array has the right name 
header = f.readline()
header = header.strip()
header = header.split(',')

# Reads line by line of the file
for line in f:
    # Gets the line up to the line-seperator
    line = line.strip()
    # Gets the column from the current line, seperated by ','
    column = line.split(',')
    # Adds the value to the list
    hjd.append(column[0])
    mag.append(column[1])

# Close file
f.close()

# ---------------------------------------------------------------------------------------
# Generates a figure for plotting with corresponding axis
fig, ax1 = plt.subplots(1, 1, figsize=(16, 9))
fig.canvas.set_window_title('Flat removal')

# Plots the data and labes the axis
ax1.plot(hjd, mag, 'xr')
ax1.set_xlabel(header[0])
ax1.set_ylabel(header[1])

# ---------------------------------------------------------------------------------------
# Generates the drawing-objects - sets the corresponding axes to the objects
rect = draw_rect(ax1)

# Enables/connects the event-handlers to the objects created above
# So that clicking/releasing/moving mouse is enabled
fig.canvas.mpl_connect('button_press_event', rect.on_press)
fig.canvas.mpl_connect('motion_notify_event', rect.on_motion)
fig.canvas.mpl_connect('button_release_event', rect.on_release)

# Changes the layout a little bit - plots are bigger/less grey space between
fig.tight_layout()

# Shows the plot - script is 'halted' until the plot is closed
plt.show()

# Writes the output to the stdout
sys.stdout.write(rect.ret())

# Exit
sys.exit(0)
# ---------------------------------------------------------------------------------------
