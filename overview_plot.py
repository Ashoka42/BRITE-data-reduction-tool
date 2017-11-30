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
# Title:		overview_plot.py
# Author:		Kondrak Matthias
# Contact:		matthias.kondrak@uibk.ac.at
# Address:		University of Innsbruck, Institute of Astro and Particle Physics
#				Technikerstrasse 25, 08/22
#				A-6020 Innsbruck
# Date:			22/08/2016
# Version:		3.0
# 
# Description:	This script plots an overview of the data in a 2x2 plot
#				1st graph is column 3 vs. column 2
#				2nd graph is column 1 vs. column 0
#				3rd graph is column 1 vs. column 2
#				4th graph is column 1 vs. column 3
#				The user can select data points with a click; properties of the
#				data are displayed; the selected data point is highlighted with
#				a circle within ALL 4 plots.
#				The user can go through the datapoints using the keyboard with
#				'n' (next) and 'p' (previous).
#				The user can delete/remove datapoints with 'CTRL' + click.
#				The corresponding values are replaced by 'nan' (not a number)
#
# Parameters:	one;
#				-----------------------------------------------------------------
#				$one   = (string) Input file
#				-----------------------------------------------------------------
#
# Input:		ASCII-file with comma-seperator ',' and a header in the first
#				line. File has a minimum of 4 columns (HJD, mag, x- and y-center)
#
# Output:		An overview plot. Clicking the data points in the plot displays
#				additional information about the data points and marks them with
#				a blue circle.
#				A string to stdout with the entries of the first column provided
#				by the input (unique identifier!) and a 'nan'-entry for deleted/removed
#				values.
#			
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


# =======================================================================================
# Defines a class for clicking
class PointBrowser(object):
    # Initialize the object
    def __init__(self, a, b, c, d):
        # Sets the variables for removal, lastindex and CTRL-key pressed
        self.remi = -1
        self.lastind = 0
        self.ctrl_on = False

        # Gets the x- and y-limit for the plots -> to not change them when redrawing
        self.ax1xlim = ax1.get_xlim()
        self.ax1ylim = ax1.get_ylim()
        self.ax2xlim = ax2.get_xlim()
        self.ax2ylim = ax2.get_ylim()
        self.ax3xlim = ax3.get_xlim()
        self.ax3ylim = ax3.get_ylim()
        self.ax4xlim = ax4.get_xlim()
        self.ax4ylim = ax4.get_ylim()

        self.ax1xlabel = ax1.get_xlabel()
        self.ax1ylabel = ax1.get_ylabel()
        self.ax2xlabel = ax2.get_xlabel()
        self.ax2ylabel = ax2.get_ylabel()
        self.ax3xlabel = ax3.get_xlabel()
        self.ax3ylabel = ax3.get_ylabel()
        self.ax4xlabel = ax4.get_xlabel()
        self.ax4ylabel = ax4.get_ylabel()

        # Also needs the specific data arrays/lists for removing single data points
        self.a = a
        self.b = b
        self.c = c
        self.d = d

        # Initializes the axis
        self.init_axis()

    # Init of axis
    def init_axis(self):

        # Sets the axis-limits
        ax1.set_xlim(self.ax1xlim)
        ax1.set_ylim(self.ax1ylim)
        ax2.set_xlim(self.ax2xlim)
        ax2.set_ylim(self.ax2ylim)
        ax3.set_xlim(self.ax3xlim)
        ax3.set_ylim(self.ax3ylim)
        ax4.set_xlim(self.ax4xlim)
        ax4.set_ylim(self.ax4ylim)

        # Sets the labels
        ax1.set_xlabel(self.ax1xlabel)
        ax1.set_ylabel(self.ax1ylabel)
        ax2.set_xlabel(self.ax2xlabel)
        ax2.set_ylabel(self.ax2ylabel)
        ax3.set_xlabel(self.ax3xlabel)
        ax3.set_ylabel(self.ax3ylabel)
        ax4.set_xlabel(self.ax4xlabel)
        ax4.set_ylabel(self.ax4ylabel)

        # For every sublot - generates a text object and sets the postion of it
        self.text1 = ax1.text(0.05, 0.95, 'selected: none', transform=ax1.transAxes, va='top')
        self.text2 = ax2.text(0.05, 0.95, 'selected: none', transform=ax2.transAxes, va='top')
        self.text3 = ax3.text(0.05, 0.95, 'selected: none', transform=ax3.transAxes, va='top')
        self.text4 = ax4.text(0.05, 0.95, 'selected: none', transform=ax4.transAxes, va='top')

        # Generates a selection object (blue circle around datapoint) and sets its
        # visibility to False
        self.selected1, = ax1.plot([c[0]], [d[0]], 'o', ms=12, alpha=0.4, color='blue', visible=False)
        self.selected2, = ax2.plot([a[0]], [b[0]], 'o', ms=12, alpha=0.4, color='blue', visible=False)
        self.selected3, = ax3.plot([c[0]], [b[0]], 'o', ms=12, alpha=0.4, color='blue', visible=False)
        self.selected4, = ax4.plot([d[0]], [b[0]], 'o', ms=12, alpha=0.4, color='blue', visible=False)

        # If the user presses a key on the keyboard

    def on_press(self, event):
        # If there was not an object selected last - leaves this function
        if self.lastind is None:
            return
        # If the key what was pressed was not 'n', 'p' or CTRL - do nothing
        if event.key not in ('n', 'p', 'control'):
            return
        # If the key pressed was 'n' - increment/next datapoint
        if event.key == 'n':
            inc = 1
        # If it was 'p' - decrement/previous datapoint
        if event.key == 'p':
            inc = -1

        if event.key == 'control':
            self.ctrl_on = True
        else:

            # The new index is the old one plus or minus 1
            self.lastind += inc
            # Clips the values, so that the index is always in the range of the array.
            # E.g. no negative values or values greater than the maximum number of elements
            # in the array are allowed!
            # Updates the variable and calls the update function
            self.lastind = np.clip(self.lastind, 0, len(self.c) - 1)
            self.update()

    # If the user releases the CTRL key
    def on_release(self, event):
        # Checks if the CTRL-key was released and sets the variable
        if event.key == 'control':
            self.ctrl_on = False
        else:
            return

            # If the user picks an object/point per click

    def on_pick(self, event):
        # Checks the lenght of the event.ind (how much points are selected)
        N = len(event.ind)
        # If nothing was selected - exits function
        if not N:
            return True
        # If multiple points were selected (overlapping of the circles) - informs
        # the user and exits function - needs to exit function, because not a
        # unique datapoints can be (or was) selected - so no informations can
        # be displayed - also throws an exception/error!
        if N > 1:
            print
            "Multiple objects lie within selection range. Zoome in to select a single object!"
            return True

            # Gets the location of that click
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata

        # Calculates the distance - hypothenuse
        distances = np.hypot(x - self.c[event.ind[0]], y - self.d[event.ind[0]])
        # Returns the indices of the minimum values along an axis
        indmin = distances.argmin()
        # Gets the data index of the selected data
        dataind = event.ind[indmin]

        # If CTRL is pressed
        if self.ctrl_on == True:
            # Sets the index of the data which should be removed
            self.remi = dataind
        else:
            # Updates the variable if no change was made
            self.lastind = dataind

        # Calls the update function and sets the removal variable afterwards
        self.update()
        self.remi = -1

    # Updating the plot
    def update(self):
        # Checks if the lastindex variable exist - if not exits function
        if self.lastind is None:
            return

        # If a point has to be removed
        if self.remi != -1:

            # Replaces the corresponding point with a NaN (Not-a-Number) entry
            self.a[self.remi] = np.nan
            self.b[self.remi] = np.nan
            self.c[self.remi] = np.nan
            self.d[self.remi] = np.nan

            # Clears all axis for new plotting
            ax1.cla()
            ax2.cla()
            ax3.cla()
            ax4.cla()

            # Plots every point
            ax1.plot(self.c, self.d, 'xr', picker=5)
            ax2.plot(self.a, self.b, 'xr', picker=5)
            ax3.plot(self.c, self.b, 'xr', picker=5)
            ax4.plot(self.d, self.b, 'xr', picker=5)
            # Also initializes the text and so on...
            self.init_axis()

        else:
            # Sets the local variable to the last index
            dataind = self.lastind

            # Sets the circles of the selected objects to True
            self.selected1.set_visible(True)
            # Sets the circle to the selected dataindex (x- and y-position of the dataindex)
            self.selected1.set_data(c[dataind], d[dataind])
            self.selected2.set_visible(True)
            self.selected2.set_data(a[dataind], b[dataind])
            self.selected3.set_visible(True)
            self.selected3.set_data(c[dataind], b[dataind])
            self.selected4.set_visible(True)
            self.selected4.set_data(d[dataind], b[dataind])

            # Generates a string for the text
            name1 = header[2] + ' : ' + str(c[dataind]) + ' ; ' + header[3] + ' : ' + str(d[dataind])
            name2 = header[0] + ' : ' + str(a[dataind]) + ' ; ' + header[1] + ' : ' + str(b[dataind])
            name3 = header[2] + ' : ' + str(c[dataind]) + ' ; ' + header[1] + ' : ' + str(b[dataind])
            name4 = header[3] + ' : ' + str(d[dataind]) + ' ; ' + header[1] + ' : ' + str(b[dataind])

            # Sets the textname
            self.text1.set_text(name1)
            self.text2.set_text(name2)
            self.text3.set_text(name3)
            self.text4.set_text(name4)
        # Redraws the figure
        fig.canvas.draw()


# =======================================================================================
# Checks if there are enough input parameters; else exits
if (len(sys.argv) < 1):
    print
    "Usage", str(sys.argv[0]), "<FILE>"
    print
    "E.g.: /str(sys.argv[0]) file.dat"
    sys.exit(1)

# Checks if file exists
if not (os.path.isfile(sys.argv[1])):
    print
    "File does not exist!"
    sys.exit(2)

# Defines the lists for plotting
a = list()
b = list()
c = list()
d = list()

# ---------------------------------------------------------------------------------------
# Opens the file - reading
f = open(sys.argv[1], 'r')

# Gets the header and splits it up, so that every entry of array header[i] has the
# column name of the corresponding column
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
    a.append(float(column[0]))
    b.append(float(column[1]))
    c.append(float(column[2]))
    d.append(float(column[3]))

# Close file
f.close()

# ---------------------------------------------------------------------------------------
# Defines a figure for plotting and adjusts some parameters
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 9))
fig.canvas.set_window_title('Overview')
plt.suptitle(
    "Overview - click on a data point to get more information - press 'n' or 'p' to go through the data points")

# Generates the 'line' element - used for selection function - with a picker of 5
# Picker: range of pixels around the datapoint for selecting this specific point
ax1.plot(c, d, 'xr', picker=5)
# Labels the axes
ax1.set_xlabel(header[2])
ax1.set_ylabel(header[3])

# Same for the other plots
ax2.plot(a, b, 'xr', picker=5)
ax2.set_xlabel(header[0])
ax2.set_ylabel(header[1])

ax3.plot(c, b, 'xr', picker=5)
ax3.set_xlabel(header[2])
ax3.set_ylabel(header[1])

ax4.plot(d, b, 'xr', picker=5)
ax4.set_xlabel(header[3])
ax4.set_ylabel(header[1])

# Generates the PointBrowser object
browser = PointBrowser(a, b, c, d)
# Connects the functions (key-pressing and mouse-picking) to the figure
fig.canvas.mpl_connect('pick_event', browser.on_pick)
fig.canvas.mpl_connect('key_press_event', browser.on_press)
fig.canvas.mpl_connect('key_release_event', browser.on_release)

# Shows the plot
plt.show()

# Writes the modified data to std-out (only first column)
sys.stdout.write(str(a))

# Exit
sys.exit(0)
# ---------------------------------------------------------------------------------------
