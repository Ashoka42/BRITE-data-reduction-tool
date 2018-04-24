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
# Title:		decor1.py
# Author:		Kondrak Matthias
# Contact:		matthias.kondrak@uibk.ac.at
# Address:		University of Innsbruck, Institute of Astro and Particle Physics
#				Technikerstrasse 25, 08/22
#				A-6020 Innsbruck
# Date:			22/08/2016
# Version:		3.0
#
# Description:	This script decorelates BRITE-data. The	user can choose between
#				decorrelation of temperature, x- and y-center.
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
    return list(map(lambda x, y: x - y, y, TMP))


# Defines a function for calculating and returning a ploynomial
def fun(x, y):
    z = np.polyfit(x, y, 1)
    return np.poly1d(z)


# =======================================================================================
# Checks if there are enough input parameters; else exits
if (len(sys.argv) < 1):
    print("Usage", str(sys.argv[0]), "<FILE>")
    print("E.g.: /str(sys.argv[0]) file.dat")
    sys.exit(1)

# Checks if file exists
if not (os.path.isfile(sys.argv[1])):
    print("File does not exist!")
    sys.exit(2)

# Defines the lists for calculation
t = list()
x = list()
y = list()
arr = list()
hjd = list()
mag = list()

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
    x.append(float(column[2]))
    y.append(float(column[3]))

# Close file
f.close()

# Loops through until user stops - result is okay
while var:

    # Calculates the correlation coefficient matrix and gets the first line
    tmp = np.corrcoef(np.array([mag, t, x, y], dtype=float))[0]
    # Calculates the r2 values in percentage
    tcorr = tmp[1] * tmp[1] * 100
    xcorr = tmp[2] * tmp[2] * 100
    ycorr = tmp[3] * tmp[3] * 100

    # Prints the result of the correlation coefficient value
    print("Correlation coefficient - higher value means stronger correlation (in percent)")
    print("T-corr:  %7.5f" % (tcorr))
    print("y-corr:  %7.5f" % (ycorr))
    print("x-corr:  %7.5f" % (xcorr))

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

    ax3.plot(y, mag, 'xr', y, fun(y, mag)(y), '-b')
    ax3.set_xlabel("YCEN[pixels]")
    ax3.set_ylabel("Magnitude[mag]")

    ax4.plot(x, mag, 'xr', x, fun(x, mag)(x), '-b')
    ax4.set_xlabel("XCEN[pixels]")
    ax4.set_ylabel("Magnitude[mag]")

    # Shows the plot
    plt.show()

    # Asks the user what to do...
    print("What do you want to do?")
    print("Enter: '1' - temperature decorrelation, '2' - y decorrelation, '3' x - decorrelation, '4' QUIT")

    # Gets the user input
    decision = raw_input("> ")
    # Checks the input and does the decorrelation
    if decision == "1":
        mag = decor(t, mag)

    elif decision == "2":
        mag = decor(y, mag)

    elif decision == "3":
        mag = decor(x, mag)

    # Exits loop if user selects '8'
    elif decision == "4":
        var = False

        # Prepares the output
        for j in range(0, len(hjd)):
            arr.append(str(hjd[j]) + ";" + str(mag[j]))

        # Writes the output to stderr!
        sys.stderr.write(str(arr))

    # Informs the user that the input was invalid
    else:
        print("Input was not valid! Please try again.")

# ---------------------------------------------------------------------------------------
# Exit
sys.exit(0)
# ---------------------------------------------------------------------------------------
