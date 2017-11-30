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
# Title:		decor2.py
# Author:		Kondrak Matthias
# Contact:		matthias.kondrak@uibk.ac.at
# Address:		University of Innsbruck, Institute of Astro and Particle Physics
#				Technikerstrasse 25, 08/22
#				A-6020 Innsbruck
# Date:			22/08/2016
# Version:		3.0
# 
# Description:	This script decorelates BRITE-data in chopping mode. The
#				user can choose between decorrelation of temperature, x- and y-
#				center and x1-, x2-, y1-, y2- center for chopping mode.
#				The modified magnitude will be written (incl. HJD for unique
#				identification) to the stdout.
#
# Parameters:	one;
#				-----------------------------------------------------------------
#				$one   = (string) Input file
#				-----------------------------------------------------------------
#
# Input:		ASCII-file with comma-seperator ','. The input needs to have at least
#				five columns. (HJD, magnitude, x, y and temperature)
#
# Output:		The output is the HJD with the corresponding magnitude written to 
#				stderr! This has to be done to get not into conflict with the 'print'
#				command, which directs its output to stdout
#
# Exit:			0 - everything went fine
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
# Defines a decorrelation function
def decor(x, y):
    # Gets the slope of the linear polynomial, which fits best
    k = fun(x, y)[1]
    # Generates a new array, which multiplies k to the whole array x
    TMP = map(lambda x: x * k, x)
    # Returns the decorrelated array (y - TMP) (subtraction)
    return map(lambda x, y: x - y, y, TMP)


# Defines a function for calculating and returning a ploynomial
def fun(x, y):
    z = np.polyfit(x, y, 1)
    return np.poly1d(z)


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

# Defines the lists for calculation
t = list()
xy = list()
x1 = list()
x2 = list()
y1 = list()
y2 = list()
arr = list()
hjd = list()
mag = list()
hjd1 = list()
hjd2 = list()
mag1 = list()
mag2 = list()

# Defines a variable for endless loop and another for writing the output file
var = True
count = 0

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
    t.append(float(column[4]))
    xy.append((float(column[2]), float(column[3])))

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

# ---------------------------------------------------------------------------------------
# Loops through until user stops - result is okay
while var:

    # Splits the data due to chopping mode - according to the curve derived above
    for i in range(0, len(xy[0])):
        # For values below that curve - adds to first list
        if xy[1][i] < f(xy[0][i]):
            hjd1.append(hjd[i])
            mag1.append(mag[i])
            x1.append(xy[0][i])
            y1.append(xy[1][i])
        # Else - adds to second list
        else:
            hjd2.append(hjd[i])
            mag2.append(mag[i])
            x2.append(xy[0][i])
            y2.append(xy[1][i])

    # Calculates the correlation coefficient matrix and gets the first line
    tmp = np.corrcoef(np.array([mag, t, xy[0], xy[1]], dtype=float))[0]
    # Calculates the r2 values in percentage
    tcorr = tmp[1] * tmp[1] * 100
    xcorr = tmp[2] * tmp[2] * 100
    ycorr = tmp[3] * tmp[3] * 100
    # Also for the chopped data
    tmp = np.corrcoef(np.array([mag1, x1, y1], dtype=float))[0]
    x1corr = tmp[1] * tmp[1] * 100
    y1corr = tmp[2] * tmp[2] * 100
    tmp = np.corrcoef(np.array([mag2, x2, y2], dtype=float))[0]
    x2corr = tmp[1] * tmp[1] * 100
    y2corr = tmp[2] * tmp[2] * 100

    # Prints the result of the correlation coefficient value
    print
    "Correlation coefficient - higher value means stronger correlation (in percent)"
    print
    "T-corr:  %7.5f" % (tcorr)
    print
    "y-corr:  %7.5f" % (ycorr)
    print
    "y1-corr: %7.5f" % (y1corr)
    print
    "y2-corr: %7.5f" % (y2corr)
    print
    "x-corr:  %7.5f" % (xcorr)
    print
    "x1-corr: %7.5f" % (x1corr)
    print
    "x2-corr: %7.5f" % (x2corr)

    # Defines a figure for plotting and adjusts some parameters
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 9))
    plt.suptitle("Correlations")

    # Plots the correlation with a linear correlation
    ax1.plot(t, mag, 'xr', t, fun(t, mag)(t), '-b')
    # Labels the axes
    ax1.set_xlabel("Temperature[C]")
    ax1.set_ylabel("Magnitude[mag]")

    # Same for the other plots
    ax2.plot(hjd, mag, 'xr')
    ax2.set_xlabel("HJD")
    ax2.set_ylabel("Magnitude[mag]")

    ax3.plot(y1, mag1, 'xr', y1, fun(y1, mag1)(y1), '-b')
    ax3.plot(y2, mag2, 'xb', y2, fun(y2, mag2)(y2), '-r')
    ax3.set_xlabel("YCEN[pixels]")
    ax3.set_ylabel("Magnitude[mag]")

    ax4.plot(x1, mag1, 'xr', x1, fun(x1, mag1)(x1), '-b')
    ax4.plot(x2, mag2, 'xb', x2, fun(x2, mag2)(x2), '-r')
    ax4.set_xlabel("XCEN[pixels]")
    ax4.set_ylabel("Magnitude[mag]")

    # Shows the plot
    plt.show()

    # Asks the user what to do...
    print
    "What do you want to do?"
    print
    "Enter: '1' - temperature decorrelation, '2' - y decorrelation, '3' y1 - decorrelation, '4' y2 - decorrelation, '5' x - decorrelation '6' x1 - decorrelation, '7' x2 - decorrelation, '8' QUIT"

    # Gets the user input
    decision = raw_input("> ")
    # Checks the input and does the decorrelation
    if decision == "1":
        mag = decor(t, mag)

    elif decision == "2":
        mag = decor(xy[1], mag)

    # If a partial decorrelation is done
    elif decision == "3":
        mag1 = decor(y1, mag1)
        # Calculates an offset, which is caused due to splitting the dataset
        offset = np.mean(mag1) - np.mean(mag2)
        j = 0
        # Loops through all data
        for i in range(0, len(mag)):
            # Checks if the dates from the two datasets matches
            if hjd[i] == hjd1[j]:
                # If they match, replaces the magnitude from the
                # whole file with the decorrelated for a specific
                # variable
                # Also subtracts the offset to compensate for unequal
                # mean of chopping-mode
                mag[i] = mag1[j] - offset
                # Increases counter, if it does not exceed length of array
                if not (j + 1) == len(hjd1):
                    j = j + 1

    # Same for decorrelation of y2, x, x1 and x2
    elif decision == "4":
        mag2 = decor(y2, mag2)
        offset = np.mean(mag1) - np.mean(mag2)
        j = 0
        for i in range(0, len(mag)):
            if hjd[i] == hjd2[j]:
                mag[i] = mag2[j] + offset
                if not (j + 1) == len(hjd2):
                    j = j + 1

    elif decision == "5":
        mag = decor(xy[0], mag)

    elif decision == "6":
        mag1 = decor(x1, mag1)
        offset = np.mean(mag1) - np.mean(mag2)
        j = 0
        for i in range(0, len(mag)):
            if hjd[i] == hjd1[j]:
                mag[i] = mag1[j] - offset
                if not (j + 1) == len(hjd1):
                    j = j + 1

    elif decision == "7":
        mag2 = decor(x2, mag2)
        offset = np.mean(mag1) - np.mean(mag2)
        j = 0
        for i in range(0, len(mag)):
            if hjd[i] == hjd2[j]:
                mag[i] = mag2[j] + offset
                if not (j + 1) == len(hjd2):
                    j = j + 1

    # Exits loop if user selects '8'
    elif decision == "8":
        var = False

        # Prepares the output
        for j in range(0, len(hjd1)):
            arr.append(str(hjd1[j]) + ";" + str(mag1[j]))

        for j in range(0, len(hjd2)):
            arr.append(str(hjd2[j]) + ";" + str(mag2[j]))

        # Writes the output to stderr to get not into conflict with the 'print'
        # command, which is directed to stdout
        sys.stderr.write(str(arr))

    # Informs the user that the input was invalid
    else:
        print
        "Input was not valid! Please try again."

    # Flushes/resets the arrays
    del x1[:]
    del y1[:]
    del hjd1[:]
    del mag1[:]
    del x2[:]
    del y2[:]
    del hjd2[:]
    del mag2[:]

# ---------------------------------------------------------------------------------------
# Exit
sys.exit(0)
# ---------------------------------------------------------------------------------------
