#!/bin/bash
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
# Title:		BRITE-DATA-REDUCTION-V3
# Author:		Kondrak Matthias
# Contact:		matthias.kondrak@uibk.ac.at
# Address:		University of Innsbruck, Institute of Astro and Particle Physics
#				Technikerstrasse 25, 08/22
#				A-6020 Innsbruck
# Date:			22/08/2016
# Version:		3.0
#
# Description:	This script guides the user through a data reduction process
#				of BRITE-satellite data. This script contains:
#				* aparature correction
#				* flat removal
#				* decorrelation
#				* sigma clipping
#
#				This script is backwords compatible using BRITE-DATA without
#				chopping-mode.
#				For more informations see README.
#
# Parameters:	one;
#				-----------------------------------------------------------------
#				$one = (string) Input file
#				-----------------------------------------------------------------
#
# Output:   	INPUT_dd.mm.yyyy_hh:mm:ss.data
#				INPUT_dd.mm.yyyy_hh:mm:ss.log
#
#
# Exit:   		0 - Success
#				1 - No input-parameter
#				2 - File does not exist
#				3 - User quit
#				4 - All datapoints removed
#
# Bash version:	4.3.42(1)
# AWK version:  4.1.3
# =======================================================================================
# ---------------------------------------------------------------------------------------
# 								VARIABLES - TO BE SET BY USER
# ---------------------------------------------------------------------------------------
# Sets the HJD, which will be subtracted from the data for easier reading
HJDsub=2456600.0

# Threshold to identify different orbits - TO BE SET BY USER
thres=0.02
# Threshold to remove whole orbit
# Orbit is removed, if more than (1-x)*datapoints within one orbit are removed
orb_thres=0.2

# ---------------------------------------------------------------------------------------
# 								System-variables
# ---------------------------------------------------------------------------------------
# Sets the default encoding to en_US.UTF-8 (US-english with UFT-8 format)
# Using '.' as decimal separator and US-timeformat.
export LC_ALL=en_US.UTF-8
# Sets the path for the temporary files to $TMPDIR or default /tmp
# If $TMPDIR is set to RAM - may improve speed, also reduces writing (important for SSD)
if [ $(env | grep -q TMPDIR) ]; then
	TD=${TMPDIR}/
else
	TD=/tmp/
fi

# Sets the name of a logfile, where all applied changes (shape of aperture correction,
# sigma clipping, decorrelation, etc.) are documented and the name of the output-file
# Filename is name of the inputfile + date (dd.mm.yyyy_hh:mm:ss) + .log/.data
strdate=$(date +"%F_%T")
LOG=$1_${strdate}.log
OUT=$1_${strdate}.data

# =======================================================================================
#									Orbital remover
# =======================================================================================
# Function to remove orbits with large scattering - input: 1st file is larger datafile
# 2nd file is file with removed data points
remove_orbit () {
	# Generates temporary files
	touch $tmpdir/TMP_OR1
	touch $tmpdir/TMP_OR2
	touch $tmpdir/TMP_OR3

	TMP_OR1=$tmpdir/TMP_OR1
	TMP_OR2=$tmpdir/TMP_OR2
	TMP_OR3=$tmpdir/TMP_OR3

	# Gets the orbits from the file and counts them (unique)
	awk -F ',' '{print $6}' $1 | uniq -c > ${TMP_OR1}
	awk -F ',' '{print $6}' $2 | uniq -c > ${TMP_OR2}

	# Gets the number of specific orbits to an array with key = orbit from the 1st dataset
	# Loops through the other dataset and checks, if the number of orbits, which were not
	# removed is larger than the original number times a threshold. If so, prints the
	# orbit
	awk -v othres=${orb_thres} 'NR == FNR {a[$2] = $1; next} {if (a[$2] * othres < $1) {print $2}}' ${TMP_OR1} ${TMP_OR2} > ${TMP_OR3}

	# Loops through the reduced dataset and only prints the orbits, which were not
	# removed in the step before.
	awk -F',' 'NR == FNR {a[$1]; next} {if ($6 in a) {print $0}}' ${TMP_OR3} $2 > ${TMP_OR1}

	# Copies result to input-file (overwriting)
	cp ${TMP_OR1} $2

	# Removes temporary files
	rm ${TMP_OR1}
	rm ${TMP_OR2}
	rm ${TMP_OR3}
}
# =======================================================================================
# ---------------------------------------------------------------------------------------
# 								Initialization
# ---------------------------------------------------------------------------------------
# Checks if AWK and Python are installed, otherwise exits;
type awk >/dev/null 2>&1 || { echo >&2 "ERROR: Awk is required but not installed. Exit."; exit 1; }
type python3 >/dev/null 2>&1 || { echo >&2 "ERROR: Python3 is required but not installed. Exit."; exit 1; }
type mktemp >/dev/null 2>&1 || { echo >&2 "ERROR: mktemp is required but not installed. Exit."; exit 1; }

# Checks if system uses 'real awk' per default
# Script not working if awk is linked to mawk (mawk-interpretor used instead of awk/gawk)
# Therefore checks version of awk
if awk --version 2>&1 | grep -q "GNU Awk"
then
	# Awk is used per default -> no problem arise
	rawk=awk
elif gawk --version 2>&1 | grep -q "GNU Awk"
then
		# Awk is used per default -> no problem arise
		echo "Using gnu awk on mac‚"
		rawk=gawk
elif awk -Wv 2>&1 | grep -q "mawk"
then
	# Mawk is used per default as 'awk'
	# Informs user about the problem
	echo "WARNING!: Your default awk-interpreter is mawk. This causes unpredictable errors. Using gawk instead."

	# Checks if gawk is present (try to use gawk instead)
	type gawk >/dev/null 2>&1 || { echo >&2 "ERROR: Gawk is also not installed. Exit."; exit 1; }

	# Prints location of gawk and uses gawk as default awk
	echo "Using gawk located in: "
	which gawk
	echo ""
	rawk=gawk
fi


# Checks the input parameter - exits the script if no file was selected
if [ $# -lt 1 ] ; then
	printf "Usage $0 \<file.dat\>\n"
	exit 1;
fi

# Checks if the input file exists - otherwise exit
if ! [ -e $1 ]; then
	echo "File does not exist!"
	exit 2;
fi

# Generates a temporary directory - first command for linux, second for OS-X
tmpdir=$(echo $(mktemp -d 2>/dev/null || mktemp -d -t 'tmp'))

# Generates temporary files and links them to variables
touch $tmpdir/TMP1
touch $tmpdir/TMP2
touch $tmpdir/TMP3
touch $tmpdir/TMP4

TMP=$tmpdir/TMP1
TEM=$tmpdir/TMP2
DATA=$tmpdir/TMP3
INPUT=$tmpdir/TMP4

# ---------------------------------------------------------------------------------------
#								Preparing data
# ---------------------------------------------------------------------------------------
# Greps the version number of the file/reduction process (R1, R2, R3, ....)
version=$(grep -e "c ReleaseV" $1 | awk -F'=' '{print $2}' | awk -F'/' '{print $1}')

# Greps the 'data-part' and also ignores 'nan'-entries from the original file
# The header starts with 'c', and 'nan' lines can not be processed in further steps
grep -v "c" $1 | grep -i -v "nan" > $INPUT

# Substracts a specified date from the HJD for easier reading/plotting
# Converts the data to magnitude, adding an arbitrary constant of 14.5
# Using CSV format with ',' as a separator for easier processing with python
$rawk -v HJD=$HJDsub '{printf("%10.6lf,%9.7lf,%5.2lf,%5.2lf,%5.2lf\n"), $1-HJD, -2.5*log($2)/log(10.0)+14.5, $3, $4, $5}' $INPUT > $TEM

# Identifies orbits/larger timespawns between shots.
# Input value is a threshold, which specifies the time between two orbits
# In the beginning sets the orbital value to 1 (n = 1)
# For the first line (NR == 1) it sets the current value to the variable. Then it
# goes on (next). If the next (current) value minus the old one is larger than the
# threshold, it adds 1 to the orbital number. Then it prints the whole line and
# adds a ',' (separator) and the current orbital number. Afterwards it sets the current
# value to the variable.
$rawk -F',' -v t=$thres 'BEGIN{n = 1} NR == 1 {old = $1; print $0","1; next} {if (($1 - old) > t) {n += 1}; {print $0","n}; old = $1} ' $TEM > $DATA

# Now the file contains various columns:
# column1 = HJD  / Heliocentric Julian Date at start of exposure (subtracted an
#			arbitrary date) [day]
# column2 = FLUX / measured flux [mag]
# column3 = XCEN / Profile centre of gravity with respect to raster origin [pixel]
# column4 = YCEN / Profile centre of gravity with respect to raster origin [pixel]
# column5 = CCDT / CCD Temperature [C]
# column6 = Orbit / Number of orbit - needed for special treatments of different orbits
#	Sigma clipping per orbit, removing whole orbit if to many datapoints are outliers

# ---------------------------------------------------------------------------------------
#							Generates logfile and informs user
# ---------------------------------------------------------------------------------------
# Generates header for the logfile
printf "! Data reduction for BRITE-data
! Version: 3.0
! Date: 22.08.2016
! Author: Matthias Kondrak
! Version 3.0 - contains:
!
!  * Aperture correction: Select data points around some x- and y-center in order to
!	remove datapoints with a large spatial scattering:
!	- Use graphical assistance/selection
!	- Use keyboard to enter values
!	- Supports chopping mode
!
!
!  * Flat removal: Remove data points above/below a certain flux/magnitude threshold
!	in order to remove large flux/magnitude scattering:
!	- Use graphical assistance/selection
!	- Use keyboard to enter values
!
!  * Decorrelation:
!	- Temperature decorrelation
!	- x- and y-center decorrelation
!	- supports x1-, x2-, y1-, y2-center decorrelation for chopping mode
!
!  * Sigma clipping: Remove data points with large orbital scattering
!	- Use keyboard to enter a clipping value
!
!
! FORMAT: ! starts comment lines
!	  * name of the variable (variable is located in the next line)
!         Contains all choices made by the user and changes of the original data
!         If user exits before the script finished: QUIT BY USER
!         File ends with END
!\n\n" > $LOG

# Generates the first entries of the logfile with the name of the file, the date and time,
# the number of lines in that file, the subtracted HJD.
printf "! General data:
* Input
%s
* Version of the input file
%s
* Date
%s
* Time
%s
* Number of lines
%i
* HJD subtraction
%s\n\n" $1 $version $(date +"%F") $(date +"%T") $(wc -l $1 | awk '{print $1}') $HJDsub >> $LOG

# Informs the user about the version and some other things
printf "BRITE-DATA-REDUCTION-V3  Copyright (C) 2016 Kondrak Matthias
This program comes with ABSOLUTELY NO WARRANTY;
This is free software, and you are welcome to redistribute it under certain conditions;
For additional informations see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------
* WELCOME *
You are using BRITE-DATA-REDUCTION v3.0
The shell will guide you through the reduction process. See the README for further information.
The data is converted to magnitudes and a HJD of '$HJDsub' is subtracted. \n"

# ---------------------------------------------------------------------------------------
#									Overview
# ---------------------------------------------------------------------------------------
# Informs the user of the next step
printf "Plotting data ...
This may take one moment.
If you are finished viewing the data, please close the plot to continue...\n"

# Adds a header to the data and only plots relevant columns for the overview
$rawk -F',' 'BEGIN{print "HJD,Magnitude[mag],XCEN[pixels],YCEN[pixels]"}{print $1","$2","$3","$4}' $DATA > $TMP

# Plots the overview - gets an array via stdout back for removign specific datapoints
# This 'expert-mode' is used in a later step, so output is not required
grommorg=$(python3 overview_plot.py $TMP)

# ---------------------------------------------------------------------------------------
#								Aparature correction
# ---------------------------------------------------------------------------------------
# =======================================================================================
# May generate a aberration correction with the data provided - using PSFC and RTSC?!?
# =======================================================================================
# Informs user and writes logfile
printf "\nShould aparature correction be applied?\n"
echo "! Data reduction - aperture correction" >> $LOG

# Checks the version (chopping mode or not)
if [ ${version} == "R1" ] || [ ${version} == "R2" ] ;
then
	# Give some choices
	select decision in "No." "Using default values." "Manually enter data with keyboard." "Using graphical assistance." "Quit."
	do
		# Select
		case $decision in
			"No.")
				# Writes to logfile, informs the user and sets default values
				echo "! No aperture correction" >> $LOG
				echo "No aperture correction"

				temparr[0]=0
				temparr[1]=100
				temparr[2]=0
				temparr[3]=100
				break ;;

			"Using default values.")
				# Writes to logfile, informs the user and sets default values
				echo "! Default values" >> $LOG
				echo "Using default values."

				# Gets the rectangular shape for the automatic aberration correction,
				# which is preselected in the file (flag-column). At first just uses
				# flagged values, then sorts the data, then extracts the first and
				# last entry. Then pastes the two values into an array
				xmin=$(awk -F "," '{if ($6 == 1) print $3}' $DATA | sort -n | head -n1)
				xmax=$(awk -F "," '{if ($6 == 1) print $3}' $DATA | sort -n | tail -n1)
				# Same for y-values
				ymin=$(awk -F "," '{if ($6 == 1) print $4}' $DATA | sort -n | head -n1)
				ymax=$(awk -F "," '{if ($6 == 1) print $4}' $DATA | sort -n | tail -n1)

				temparr[0]=$xmin
				temparr[1]=$xmax
				temparr[2]=$ymin
				temparr[3]=$ymax
				break ;;

			# Input by keyboard
			"Manually enter data with keyboard.")
				# Informs user and writes logfile
				echo "Please enter your aperture values for x: x_min x_max y_min y_max"
				echo "e.g.: 0 100 0 100"
				echo "! Using keyboard data" >> $LOG
				# Reads the input
				read user_input
				# Sets the user input to an array
				temparr=(${user_input})
				break ;;

			# For graphical assistance
			"Using graphical assistance.")
				# Informs user and writes logfile
				echo "Please use your mouse to draw a rectangular around the data, which you want to use. When you are finish, close the window to continue."
				echo "! Using graphical data" >> $LOG

				# Gets the rectangular shape for the automatic aberration correction,
				# which is preselected in the file (flag-column). At first just uses
				# flagged values, then sorts the data, then extracts the first and
				# last entry. Then pastes the two values into an array
				xmin=$(awk -F "," '{if ($6 == 1) print $3}' $DATA | sort -n | head -n1)
				xmax=$(awk -F "," '{if ($6 == 1) print $3}' $DATA | sort -n | tail -n1)
				# Same for y-values
				ymin=$(awk -F "," '{if ($6 == 1) print $4}' $DATA | sort -n | head -n1)
				ymax=$(awk -F "," '{if ($6 == 1) print $4}' $DATA | sort -n | tail -n1)

				# Calls the aperture_correction script, if the script is closed
				# the values of the rectangle are returned to an array
				temparr=($(python3 aperture_correction1.py $TMP $xmin $xmax $ymin $ymax))
				break ;;

			# Quit
			"Quit.")
				# Informs user and writes logfile
				echo "Exit"
				echo "QUIT BY USER" >> $LOG
				exit 3 ;;
		esac
	done
else
	# Give some choices
	select decision in "No." "Manually enter data with keyboard." "Using graphical assistance." "Quit."
	do
		# Select
		case $decision in
			"No.")
				# Writes to logfile, informs the user and sets default values
				echo "! Default values" >> $LOG
				echo "Using default values."
				temparr[0]=0
				temparr[1]=100
				temparr[2]=0
				temparr[3]=100
				temparr[4]=0
				temparr[5]=100
				temparr[6]=0
				temparr[7]=100
				break ;;

			# Input by keyboard
			"Manually enter data with keyboard.")
				# Same for chopping mode
				echo "Please enter your aperture values for x: x1_min x1_max y1_min y1_max x2_min x2_max y2_min y2_max"
				echo "e.g.: 0 100 0 100 0 100 0 100"
				echo "! Using keyboard data" >> $LOG
				read user_input
				temparr=(${user_input})
				break ;;

			# For graphical assistance
			"Using graphical assistance.")
				# Informs user and writes logfile
				echo "Please use your mouse to draw a rectangular around the data, which you want to use. When you are finish, close the window to continue."
				echo "! Using graphical data" >> $LOG
				# Does the aperture_correction
				temparr=($(python3 aperture_correction2.py $TMP))
				break ;;

			# Quit
			"Quit.")
				# Informs user and writes logfile
				echo "Exit"
				echo "QUIT BY USER" >> $LOG
				exit 3 ;;
		esac
done
fi

# Informs the user and writes logfile
echo ""
echo "Applying aperture corrections with:" ${temparr[@]}
# Checks for version
if [ ${version} == "R1" ] || [ ${version} == "R2" ] ;
then
	echo "* x_min x_max y_min y_max" >> $LOG
	# Applies the aperture correction - excludes data points, which are not within the
	# rectangular shape(s) selected above - needs echo, because of the 'default'
	# option in R1/R2
	$rawk -F "," -v x0=$(echo ${temparr[0]}) -v x1=$(echo ${temparr[1]}) -v y0=$(echo ${temparr[2]}) -v y1=$(echo ${temparr[3]}) '{if ($3 > x0 && $3 < x1 && $4 > y0 && $4 < y1) print $0}' OFS=',' $DATA > $TEM
else
	echo "* x1_min x1_max y1_min y1_max x2_min x2_max y2_min y2_max" >> $LOG
	$rawk -F "," -v x0=${temparr[0]} -v x1=${temparr[1]} -v y0=${temparr[2]} -v y1=${temparr[3]} -v x2=${temparr[4]} -v x3=${temparr[5]} -v y2=${temparr[6]} -v y3=${temparr[7]} '{if (($3 > x0 && $3 < x1 && $4 > y0 && $4 < y1) || ($3 > x2 && $3 < x3 && $4 > y2 && $4 < y3)) print $0}' OFS=',' $DATA > $TEM
fi
echo ${temparr[@]} >> $LOG

# Removes orbits, where to much data points were removed
remove_orbit $DATA $TEM

# Checks if the number of lines of the file is 0
if [ $(wc -l $TEM | awk '{print $1}') == 0 ] ;
then
	# Informs the user and writes logfile
	echo ""
	echo "No values are contained in the ouput. Please try again with other parameters."
	echo "!No values are contained using these parameters." >> $LOG
	echo "END" >> $LOG
	# Exits
	exit 4;
fi

# ---------------------------------------------------------------------------------------
#									FLAT REMOVAL
# ---------------------------------------------------------------------------------------
# Informs user and writes logfile
echo "Would you like to do a flat removal?"
printf "\n\n! Data reduction - Flat removal\n" >> $LOG

# Give some choices
select decision in "No." "Manually enter data with keyboard." "Using graphical assistance." "Quit."
do
	case $decision in
		"No.")
			# Writes logfile
			echo "! No removal" >> $LOG
			# Sets default values
			flat[0]=-50
			flat[1]=50
			break ;;

		"Manually enter data with keyboard.")
			# Informs user and writes logfile
			echo "Please enter a lower and upper threshold."
			echo "! Using keyboard input" >> $LOG
			# Reads the input and sets it to an array
			read user_input
			flat=(${user_input})
			break ;;

		"Using graphical assistance.")
			# Informs user and writes logfile
			echo "Please use your mouse select data, close plot to continue"
			echo "! Using graphical data" >> $LOG
			# Prepares the data - header and just uses needed columns
			$rawk -F',' 'BEGIN{print "HJD,Magnitude[mag]"}{print $1","$2}' $TEM > $TMP
			# Starts flat_removal script; Returns lower and upper limit
			flat=($(python3 flat_removal.py $TMP))
			break ;;

		"Quit.")
			# Informs user and writes logfile
			echo "Exit"
			echo "QUIT BY USER" >> $LOG
			exit 3 ;;
	esac
done

# Informing user and writing logfile
echo "Applying flat removal with:" ${flat[@]}
echo "* Flat removal:" >> $LOG
echo ${flat[@]} >> $LOG

# Removes all data points, which are below/above the lower/upper threshold
$rawk -F "," -v m0=${flat[0]} -v m1=${flat[1]} '{if ($2 > m0 && $2 < m1) print $0}' OFS=',' $TEM > $TMP

# Removes orbits, where to much datapoints are thrown away
remove_orbit $DATA $TMP

# Checks if the number of lines of the file is 0
if [ $(wc -l $TMP | awk '{print $1}') == 0 ] ;
then
	# Informs the user and writes logfile
	echo ""
	echo "No values are contained in the ouput. Please try again with other parameters."
	echo "!No values are contained using these parameters." >> $LOG
	echo "END" >> $LOG
	# Exits
	exit 4;
fi

# ---------------------------------------------------------------------------------------
#								Decorrelation
# ---------------------------------------------------------------------------------------
# Informs user
echo ""
echo "Starting decorrelation ..."
# Starts the decorrelation script depending on the version of the file
if [ ${version} == "R1" ] || [ ${version} == "R2" ] ;
then
	# The output of the decorrelation is redirected to the stderr, so here the stderr
	# from the decorrelation is redirected to the temporary/working file
	python3 decor1.py $TMP 2> $TEM
else
	python3 decor2.py $TMP 2> $TEM
fi

# To get the required syntax the ',' is replaced by a new line, the '[' and ']' denoting
# the array are removed, '>' from the output mechanism and ' ' are removed and the
# ';' is replaced by a ',' for using the same separator in all files
sed -i -e "s/,/\n/g; s/>//g; s/\[//; s/\]//; s/ //g; s/'//g; s/;/,/g" $TEM
sort -o $TEM $TEM

# Now the 2 files (before and after the decorrelation) are compared. The magnitude of
# the decorrelated file is put to a variable with index of the HJD. Then the second
# (original) file is read and the 2nd column (magnitude) is replaced by the
# value from the first file, if the key/index (HJD) is present. Then the whole column
# is printed.
$rawk -F ',' 'FNR == NR {a[$1] = $2; next}{if (a[$1] != ""){$2 = a[$1]; print $0}}' OFS=',' $TEM $DATA > $TMP

# ---------------------------------------------------------------------------------------
#								Sigma clipping
# ---------------------------------------------------------------------------------------
# Informs user and writes logfile
printf "\nShould sigma-clipping be performed?\n"
echo "! Data reduction - sigma clipping" >> $LOG

# Again choices to make
select decision in "No." "Yes." "Quit."
do
	case $decision in
		"No.")
			# Informs user and writes logfile
			echo "! No clipping applied" >> $LOG
			break ;;

		"Yes.")
			# Informs user and writes logfile
			echo "Please enter the threshold value for clipping."
			echo "! Sigma clipping applied" >> $LOG
			# Gets user input and sets variable
			read user_input
			sig=(${user_input})
			# Writes data to logfile
			echo "* sigma" >> $LOG
			echo $sig >> $LOG
			cp $TMP WORKING.dat
			# Applies sigma clipping - NAN entries
			python3 sigma.py $TMP $sig > $TEM
			# Removes all "NAN" entries
			sed -i -e "s/'//g; s/,/\n/g; s/\[//g; s/\]//g; s/ //g" $TEM
			paste -d',' $TMP $TEM > $TEM.tmp
			grep -i -v 'nan' $TEM.tmp > $TMP
			break ;;

		"Quit.")
			# Informs user and writes logfile
			echo "Exit"
			echo "QUIT BY USER" >> $LOG
			exit 3 ;;
	esac
done

# Removes orbits, where to much datapoints are thrown away
remove_orbit $DATA $TEM

# Checks if the number of lines of the file is 0
if [ $(wc -l $TEM | awk '{print $1}') == 0 ] ;
then
	# Informs the user and writes logfile
	echo ""
	echo "No values are contained in the ouput. Please try again with other parameters."
	echo "!No values are contained using these parameters." >> $LOG
	echo "END" >> $LOG
	# Exits
	exit 4;
fi

# ---------------------------------------------------------------------------------------
#									Final overview
# ---------------------------------------------------------------------------------------
# Informs the user
echo ""
echo "Plotting the output. Individual points may be removed now."

# Adds a header to the data for plotting
$rawk -F',' 'BEGIN{print "HJD,Magnitude[mag],XCEN[pixels],YCEN[pixels]"}{print $1","$2","$3","$4}' $TMP > $TEM

# Plots the overview - gets an array via stdout back for removign specific datapoints
# The user can now remove/delete single points using 'CTRL' + click
python3 overview_plot.py $TEM > $TEM.tmp
# Since overview_plot.py returns an array of the HJD incl. removed points denoted by
# 'nan' (not a number), the array has to be trimed (removing of '[' and ']' and
# replacing ',' by a new line
sed -i -e "s/\[//g; s/\]//g; s/,/\n/g" $TEM.tmp
# Then the array is pasted along the original one and lines containing a 'nan' entry
# are removed
paste -d',' $TMP $TEM.tmp > $TEM
grep -i -v 'nan' $TEM > $TMP

# Checks if the number of lines of the file is 0
if [ $(wc -l $TMP | awk '{print $1}') == 0 ] ;
then
	# Informs the user and writes logfile
	echo ""
	echo "No values are contained in the ouput. Please try again with other parameters."
	echo "!No values are contained using these parameters." >> $LOG
	echo "END" >> $LOG
	# Exits
	exit 4;
fi

# ---------------------------------------------------------------------------------------
#									Output
# ---------------------------------------------------------------------------------------
# For each version: appends the 'missing columns' like flags or blurring coefficient
# And also generates a header for the column descriptions
if [ ${version} == "R1" ] || [ ${version} == "R2" ] ;
then
	# Generates the header
	header=$(echo "HJD,Magnitude[mag],XCEN[pixels],YCEN[pixels],Temperature[C],Orbit,Aperture_flag")
	# Gets the missing columns from the input file and formats them
	# Using a ',' to separate the first column from the rest with ';'
	$rawk -v HJD=$HJDsub '{printf("%10.6lf,%i\n"), $1-HJD, $7}' $INPUT > $TEM
elif [ ${version} == "R3" ];
then
	header=$(echo "HJD,Magnitude[mag],XCEN[pixels],YCEN[pixels],Temperature[C],Orbit,PSCF,RTSC")
	$rawk -v HJD=$HJDsub '{printf("%10.6lf,%8.6lf;%6.2lf\n"), $1-HJD, $7, $8}' $INPUT > $TEM
elif [ ${version} == "R4" ];
then
	header=$(echo "HJD,Magnitude[mag],XCEN[pixels],YCEN[pixels],Temperature[C],Orbit,PSFC1,PSFC2,RTSC")
	$rawk -v HJD=$HJDsub '{printf("%10.6lf,%8.6lf;%8.6lf;%6.2lf\n"), $1-HJD, $7, $8, $9}' $INPUT > $TEM
fi

# Prints the header of the output file
echo $header > $OUT
# Uses the first column from the TEM file as an index for the array 'a'. The values
# of the array are given by the other entries in the corresponding line of TEM.
# Loops through the second file (TMP) and appends the missing columns, if the
# index of the array is correct
$rawk -F ',' 'FNR == NR {a[$1] = $2; next}{if (a[$1] != ""){print $1,$2,$3,$4,$5,$6,a[$1]}}' OFS=',' $TEM $TMP >> $OUT

# Replaces the ';' by a ',' for an equal separator and removes all ' '.
sed -i -e 's/;/,/g; s/ //g' $OUT

# Informs the user of the ouput.
echo ""
echo "Output:" $OUT

# Logs the END of this script
printf "\nEND" >> $LOG
# Delets temporary directory
rm -r $tmpdir

# Informs the user that the script is finished
echo "Script finished"
exit 0;
# ---------------------------------------------------------------------------------------
