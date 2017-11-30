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
# Title:		aperture_correction1.py
# Author:		Kondrak Matthias
# Contact:		matthias.kondrak@uibk.ac.at
# Address:		University of Innsbruck, Institute of Astro and Particle Physics
#				Technikerstrasse 25, 08/22
#				A-6020 Innsbruck
# Date:			22/08/2016
# Version:		3.0
# 
# Description:	This script is the graphical help for choosing the aperture
#				correction. The data is plotted in 2x3 graphs and supports
#				chopping-mode data. Each line contains the following plots
#				for one position of the chopping:
#				y-center (y-axis) vs. x-center (x-axis),
#				mag vs. x-center,
#				mag vs. y-center.
#				The user can draw rectangles in the left plots, choosing data with
#				specific x- and y-center positions for further calculations and 
#				processing. These values can also be adjusted drawing 'borders' 
#				in the middle and/or right plot.
#				Closing the plot will return the choosen values to stdout.
#
# Parameters:	one;
#				-----------------------------------------------------------------
#				$one = (string) Input file
#				-----------------------------------------------------------------
#
# Input:		ASCII-file with comma-seperator ',' and a header in the first
#				line. Mouse-input by drawing inside the plot.
#				File needs a minimum of 4 columns (HJD, mag, x-center, y-center)
#
# Output:	 	String to STDOUT with the borders of the rectangles
#				X0 X1 Y0 Y1 X2 X3 Y2 Y3
#				Separated by space ' ' and sorted.
#				Default value is 0 100 0 100 0 100 0 100
#
# Exit:   		0 - everything went fine
#				1 - not enough input-parameters
#				2 - file does not exist
#
# =======================================================================================
# Import libraries
import sys
import os.path
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# =======================================================================================
# Defines the class for drawing rectangulars and selecting the aperture correction
class draw_rect(object):
    # Initialises all variables
    def __init__(self, axis1, axis2, axis3):
        # Defines the axes used by the specific rectangular for each plot
        self.ax1 = axis1
        self.ax2 = axis2
        self.ax3 = axis3
        # Defines the default x- and y-coordinates of the rectangle
        self.x0 = 0
        self.x1 = 100
        self.y0 = 0
        self.y1 = 100
        # Defines variables for drawing enabled/disabled and default output
        self.draw = 0
        self.output = "0 100 0 100"
        # Defines a rectangule for each plot (y/x, mag/x, mag/y)
        # Sets the rectangular y-position for the 2nd and 3rd plot to -50 and the
        # height to 100 to cover a default range from -50 to 50 (mag)
        self.rect1 = Rectangle((0, 0), 0, 0, alpha=0.3)
        self.rect2 = Rectangle((0, -50), 0, 100, alpha=0.3)
        self.rect3 = Rectangle((0, -50), 0, 100, alpha=0.3)
        # Connects the rectangle to the corresponding axis
        self.ax1.add_patch(self.rect1)
        self.ax2.add_patch(self.rect2)
        self.ax3.add_patch(self.rect3)

    # Defines what should happen on pressing the mousebutton
    def on_press(self, event):
        # If not the left mouse-button is clicked
        if event.button != 1:
            return
        # At first gets the toolbar-button-pressed-info and checks if
        # no toolbar-button is pressed (zooming or shifting/moving)
        # So zooming/shifting can be done without drawing a new rectangle
        toolbar = plt.get_current_fig_manager().toolbar
        if toolbar.mode == '':

            # If the mouse is INSIDE the first plot
            if event.inaxes is self.ax1:
                # Enables drawing
                self.draw = 1
                # Gets the coordinates, where the mouse it located
                self.x0 = event.xdata
                self.y0 = event.ydata

            # If mouse in second plot
            if event.inaxes is self.ax2:
                # Enables drawing just of the x-coord (y=mag is static)
                self.draw = 1
                self.x0 = event.xdata

            if event.inaxes is self.ax3:
                # Enables drawing just of the x-coord (y=mag is static)
                # This x-value of this graph is the actual y-center-value
                self.draw = 1
                self.y0 = event.xdata

    # On mouse button release
    def on_release(self, event):
        # If not the left mouse-button is clicked
        if event.button != 1:
            return

        toolbar = plt.get_current_fig_manager().toolbar
        if toolbar.mode == '':
            # Only if it is INSIDE the first plot
            if event.inaxes is self.ax1:
                # Sorts the x- and y-values and sets the output variable
                if self.x0 < self.x1:
                    if self.y0 < self.y1:
                        self.output = str(self.x0) + " " + str(self.x1) + " " + str(self.y0) + " " + str(self.y1)
                    else:
                        self.output = str(self.x0) + " " + str(self.x1) + " " + str(self.y1) + " " + str(self.y0)
                else:
                    if self.y0 < self.y1:
                        self.output = str(self.x1) + " " + str(self.x0) + " " + str(self.y0) + " " + str(self.y1)
                    else:
                        self.output = str(self.x1) + " " + str(self.x0) + " " + str(self.y1) + " " + str(self.y0)

                # Redraws the rectangle (if last drawing went wrong)
                # Width and height is just difference of x- and y-values
                self.rect1.set_width(self.x1 - self.x0)
                self.rect1.set_height(self.y1 - self.y0)
                self.rect1.set_xy((self.x0, self.y0))
                self.rect2.set_width(self.x1 - self.x0)
                self.rect2.set_xy((self.x0, 0))
                self.rect3.set_width(self.y1 - self.y0)
                self.rect3.set_xy((self.y0, 0))
                self.ax1.figure.canvas.draw()
                self.ax2.figure.canvas.draw()
                self.ax3.figure.canvas.draw()

            # Same for 2nd plot, but just with x:
            if event.inaxes is self.ax2:
                if self.x0 < self.x1:
                    self.output = str(self.x0) + " " + str(self.x1) + " " + str(self.y0) + " " + str(self.y1)
                else:
                    self.output = str(self.x1) + " " + str(self.x0) + " " + str(self.y0) + " " + str(self.y1)

                # Just needs to redraw the first and second graph, because
                # the 3rd was not modified if drawing in the second graph
                self.rect1.set_width(self.x1 - self.x0)
                self.rect1.set_xy((self.x0, self.y0))
                self.rect2.set_width(self.x1 - self.x0)
                self.rect2.set_xy((self.x0, -50))
                self.ax1.figure.canvas.draw()
                self.ax2.figure.canvas.draw()

            # Same for 3rd plot, but not x-value == y-center-value
            if event.inaxes is self.ax3:
                if self.y0 < self.y1:
                    self.output = str(self.x0) + " " + str(self.x1) + " " + str(self.y0) + " " + str(self.y1)
                else:
                    self.output = str(self.x0) + " " + str(self.x1) + " " + str(self.y1) + " " + str(self.y0)

                # Redraws first and 3rd graph, because 2nd was not modified
                self.rect1.set_height(self.y1 - self.y0)
                self.rect1.set_xy((self.x0, self.y0))
                self.rect3.set_width(self.y1 - self.y0)
                self.rect3.set_xy((self.y0, -50))
                self.ax1.figure.canvas.draw()
                self.ax3.figure.canvas.draw()

            # Sets drawing to 0 - disables drawing after mouse was released
            self.draw = 0

    # On mouse motion
    def on_motion(self, event):
        # If not the left mouse-button is clicked
        if event.button != 1:
            return
        # If drawing is enabled
        if (self.draw == 1):
            # Checks if inside first plot
            if event.inaxes is self.ax1:
                # Gets the new (because moving) coordinates, calculates
                # the rectangle and draws it during moving/motion
                # Also redraws other plots - this can take a few seconds
                # due to calulating/redrawing - depends on data
                self.x1 = event.xdata
                self.y1 = event.ydata
                self.rect1.set_width(self.x1 - self.x0)
                self.rect1.set_height(self.y1 - self.y0)
                self.rect1.set_xy((self.x0, self.y0))
                self.rect2.set_width(self.x1 - self.x0)
                self.rect2.set_xy((self.x0, -50))
                self.rect3.set_width(self.y1 - self.y0)
                self.rect3.set_xy((self.y0, -50))
                self.ax1.figure.canvas.draw()
                self.ax2.figure.canvas.draw()
                self.ax3.figure.canvas.draw()

            # Same for second plot - but just redraws first and 2nd graph
            elif event.inaxes is self.ax2:
                self.x1 = event.xdata
                self.rect1.set_width(self.x1 - self.x0)
                self.rect1.set_xy((self.x0, self.y0))
                self.rect2.set_width(self.x1 - self.x0)
                self.rect2.set_xy((self.x0, -50))
                self.ax1.figure.canvas.draw()
                self.ax2.figure.canvas.draw()

            # Same for third plot
            elif event.inaxes is self.ax3:
                self.y1 = event.xdata
                self.rect1.set_height(self.y1 - self.y0)
                self.rect1.set_xy((self.x0, self.y0))
                self.rect3.set_width(self.y1 - self.y0)
                self.rect3.set_xy((self.y0, -50))
                self.ax1.figure.canvas.draw()
                self.ax3.figure.canvas.draw()
            # If mouse moves outside a graph - disables drawing
            else:
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
xy = list()
x1 = list()
y1 = list()
x2 = list()
y2 = list()
mag = list()
mag1 = list()
mag2 = list()

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
    # Adds the values to the corresponding list
    # Gets the xy values in ONE list/array as a tuple
    xy.append((float(column[2]), float(column[3])))
    mag.append(column[1])

# Close file
f.close()

# ---------------------------------------------------------------------------------------
# Makes an array out of that list and transposes it, so that x values are in column 0
# and y values are in column 1
xy = np.array(xy)
xy = xy.T

# Needs to calculate a border to distinguish the chopped data points
# Chopping may not always be aligned in x- or y-direction only, but in a combindation 
# of both. So chopped data may be 'diagonal'
# To distinguish them properly, an algorithm for detection is written here
# Calculates the mean value of the x- and y-center
xmean = np.mean(xy[0])
ymean = np.mean(xy[1])

# Calculates the best (least squared) fit of a linear regression of the whole data
# and gets the slope k (0th entry)
k = np.polyfit(xy[0], xy[1], 1)[0]

# Checks if k is 0 -> if so, changes its value, because of division through zero might
if k == 0:
    k = 0.0001

# Generates an orthogonal function using ordinary linear equation: y = k * x + d
# Deriving (orthogonal k') using k * k' = -1 -> k' = -1/k
# Using xmean and ymean for deriving a 'footpoint' of the new curve
# d = y - k' * x = y + 1/k * x
d = ymean + 1 / k * xmean

# Generates the polynomial which describes the curve for splitting/chopping
f = np.poly1d([-1 / k, d])

# Splits the data due to chopping mode - according to the curve derived above
for i in range(0, len(xy[0])):
    # For values below that curve - adds to first list
    if xy[1][i] < f(xy[0][i]):
        x1.append(xy[0][i])
        y1.append(xy[1][i])
        mag1.append(mag[i])
    # Else - adds to second list
    else:
        x2.append(xy[0][i])
        y2.append(xy[1][i])
        mag2.append(mag[i])

# ---------------------------------------------------------------------------------------
# Generates a figure for plotting with corresponding axis and title
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(16, 9))
fig.canvas.set_window_title('aperture correction')

# Plots the data and labes the axis
ax1.plot(x1, y1, 'xr')
ax1.set_xlabel(header[2])
ax1.set_ylabel(header[3])

# Also sets a title at the second graph
ax2.plot(x1, mag1, 'xr')
ax2.set_xlabel(header[2])
ax2.set_ylabel(header[1])
ax2.set_title("aperture correction of first chopping data")

ax3.plot(y1, mag1, 'xr')
ax3.set_xlabel(header[3])
ax3.set_ylabel(header[1])

ax4.plot(x2, y2, 'xr')
ax4.set_xlabel(header[2])
ax4.set_ylabel(header[3])

ax5.plot(x2, mag2, 'xr')
ax5.set_xlabel(header[2])
ax5.set_ylabel(header[1])
ax5.set_title("aperture correction of second chopping data")

ax6.plot(y2, mag2, 'xr')
ax6.set_xlabel(header[3])
ax6.set_ylabel(header[1])

# ---------------------------------------------------------------------------------------
# Generates the drawing-objects - sets the corresponding axes to the objects
upper = draw_rect(ax1, ax2, ax3)
lower = draw_rect(ax4, ax5, ax6)

# Enables/connects the event-handlers to the objects created above so that
# clicking/releasing/moving mouse is enabled
fig.canvas.mpl_connect('button_press_event', upper.on_press)
fig.canvas.mpl_connect('motion_notify_event', upper.on_motion)
fig.canvas.mpl_connect('button_release_event', upper.on_release)
fig.canvas.mpl_connect('button_press_event', lower.on_press)
fig.canvas.mpl_connect('motion_notify_event', lower.on_motion)
fig.canvas.mpl_connect('button_release_event', lower.on_release)

# Changes the layout a little bit - plots are bigger/less grey space between
fig.tight_layout()

# Shows the plot - script is 'halted' until the plot is closed
plt.show()

# Gets the output string and writes it to the stdout
output = upper.ret() + ' ' + lower.ret()
sys.stdout.write(output)

# Exit
sys.exit(0)
# ---------------------------------------------------------------------------------------
