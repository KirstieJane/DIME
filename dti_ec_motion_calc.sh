#!/bin/bash

# Script written by Mark Jenkinson of FMRIB; posted in reply to Charlotte's question on the FSL forum 11/01/2011
# https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;4f8352a3.1101

# Commented and edited by Kirstie Whitaker (kw401@cam.ac.uk) on 3rd February 2014
# Specifically added the ability to not calculate the difference
# between b0 images and diffusion weighted images

# Define usage function
if [ $# -lt 2 ] ; then 
  echo "Usage: `basename $0` <eddy current ecclog file> <bvals file>"
  exit 1;
fi

# Read in ecclog file and the bvals file
logfile=$1
bvalsfile=$2

# Define the directory in which the ecclog file resides
dwi_dir=`dirname ${logfile}`

# Remove files that are created by this script
rm -f ${dwi_dir}/ec_*.txt

# Get the basename (which is also the basename of the eddy corrected nifti file)
basenm=`basename $logfile .ecclog`

# Find the line numbers for the line before each of the registration matrices
# by searching for the word "Final", which precedes 4 rows that look like:
#    Final result:                                  line number: n
#    1.000000 0.000000 0.000000 0.000000            line number: n1
#    0.000000 1.000000 0.000000 0.000000 
#    0.000000 0.000000 1.000000 0.000000 
#    0.000000 0.000000 0.000000 1.000000            line number: n2
# or
#     Final result:                                 line number: n
#     1.011887 -0.002517 0.010715 -1.205730         line number: n1
#     0.002287 1.007846 0.002932 -1.780132
#     -0.012578 -0.000272 1.002718 1.177598 
#     0.000000 0.000000 0.000000 1.000000           line number: n2
#
# The variable nums contains the line numbers of the word "Final"
nums=(`grep -n 'Final' $logfile | sed 's/:.*//'`)

# The variable bvals contains the bvalue for each of these transformations
bvals=(`cat ${bvalsfile}`)

# Set a variable saying that we're starting so we know to keep the next file as
# the first volume for calculating the absolute movement
firsttime=yes

# Set two variables saying whether we've found the first diffusion weighted image
# or the first b0 image
foundfirstnotb0=no
foundfirstb0=no

# Set a counter going for each registration matrix
i=0

# Loop through each of the registration volumes
while [[ $i -lt  ${#bvals[@]} ]]; do
    n=${nums[i]} 
    b=${bvals[i]}

    # Figure out if you're looking at a b0 volume or not
    if [[ `echo "${b} < 50" | bc` == 1 ]]; then
        b0=yes
    else
        b0=no
    fi

    # Calculate the line numbers of the first (n1) and last (n2) lines
    # for the registration matrix (see above for example)
    n1=`echo $n + 1 | bc` ; 
    n2=`echo $n + 5 | bc` ;
    # Write the matrix into the current.mat file
    # note that this overwrites the file each time
    sed -n  "$n1,${n2}p" $logfile > ${dwi_dir}/current.mat ; 

    # If this is the first time you're running this loop
    # then save the matrix as first.mat, change the firsttime marker
    # and copy over the matrix to previous.mat
    # The first.mat will be the same for all the rest of the comparisons
    # and the previous.mat will be the one directly before the current matrix
    if [[ $firsttime = yes ]]; then
        firsttime=no
        cp ${dwi_dir}/current.mat ${dwi_dir}/first.mat
        cp ${dwi_dir}/current.mat ${dwi_dir}/previous.mat
    fi

    # If you haven't yet found the first non b0
    # then save this one as the firstnotb0.mat
    if [[ $foundfirstnotb0 == no && ${b0} == 'no' ]]; then
        cp ${dwi_dir}/current.mat ${dwi_dir}/firstnotb0.mat
        cp ${dwi_dir}/current.mat ${dwi_dir}/previousnotb0.mat
        foundfirstnotb0=yes
    fi

    # If you haven't yet found the first b0
    # then save this one as the firstb0.mat
    if [[ $foundfirstb0 == no && ${b0} == 'yes' ]]; then
        cp ${dwi_dir}/current.mat ${dwi_dir}/firstb0.mat
        cp ${dwi_dir}/current.mat ${dwi_dir}/previousb0.mat
        foundfirstb0=yes
    fi



    # Now calculate the root mean square difference between the
    # current matrix and the four possibilities:
    ### ABSOLUTE DIFFERENCE COMPARED TO FIRST
    absval=`$FSLDIR/bin/rmsdiff ${dwi_dir}/current.mat ${dwi_dir}/first.mat ${dwi_dir}/$basenm`

    ### RELATIVE DIFFERENCE COMPARED TO PREVIOUS
    relval=`$FSLDIR/bin/rmsdiff ${dwi_dir}/current.mat ${dwi_dir}/previous.mat ${dwi_dir}/$basenm`

    # If you've found the first b0, then we can compare the 
    # current volume to the first and previous not b0 volumes
    # Note that this is only meaningful for the non-b0 volumes
    # but we'll run all of the comparisons for ease.
    # Just don't consider them when you calculate the means etc!
    if [[ ${b0} == 'no' ]]; then
        
        ### ABSOLUTE DIFFERENCE COMPARED TO FIRST NON-B0
        absvalnotb0=`$FSLDIR/bin/rmsdiff ${dwi_dir}/current.mat ${dwi_dir}/firstnotb0.mat ${dwi_dir}/$basenm`

        ### RELATIVE DIFFERENCE COMPARED TO PREVIOUS NON-B0

        relvalnotb0=`$FSLDIR/bin/rmsdiff ${dwi_dir}/current.mat ${dwi_dir}/previousnotb0.mat ${dwi_dir}/$basenm`

        ### ABSOLUTE DIFFERENCE COMPARED TO FIRST B0
        # doesn't mean anything (will be subject to artefact)
        absvalb0='.'

        ### RELATIVE DIFFERENCE COMPARED TO PREVIOUS B0
        # doesn't mean anything (will be subject to artefact)
        relvalb0='.'

    elif [[ ${b0} == 'yes' ]]; then
        ### ABSOLUTE DIFFERENCE COMPARED TO FIRST NON-B0
        # doesn't mean anything (will be subject to artefact)
        absvalnotb0='.'

        ### RELATIVE DIFFERENCE COMPARED TO PREVIOUS NON-B0
        # doesn't mean anything (will be subject to artefact)
        relvalnotb0='.'

        ### ABSOLUTE DIFFERENCE COMPARED TO FIRST B0
        absvalb0=`$FSLDIR/bin/rmsdiff ${dwi_dir}/current.mat ${dwi_dir}/firstb0.mat ${dwi_dir}/$basenm`

        ### RELATIVE DIFFERENCE COMPARED TO PREVIOUS B0
        relvalb0=`$FSLDIR/bin/rmsdiff ${dwi_dir}/current.mat ${dwi_dir}/previousb0.mat ${dwi_dir}/$basenm`
    else
        ### None of the values will mean anything so don't fill them in
        absvalnotb0='.'
        relvalnotb0='.'
        absvalb0='.'
        relvalb0='.'         
    fi
    
    # Write the absolute and relative rms values into the ec_disp.txt file
    echo $absval $relval >> ${dwi_dir}/ec_disp.txt
    echo $absvalnotb0 $relvalnotb0 >> ${dwi_dir}/ec_disp_notb0.txt
    echo $absvalb0 $relvalb0 >> ${dwi_dir}/ec_disp_b0.txt

    # Copy over the current matrix to previous.mat ready for the next loop
    cp ${dwi_dir}/current.mat ${dwi_dir}/previous.mat

    # If the current matrix does not represent a b0 acquisition
    # then copy this to previousnotb0.mat
    if [[ $b0 == no ]]; then
        cp ${dwi_dir}/current.mat ${dwi_dir}/previousnotb0.mat
    else
        cp ${dwi_dir}/current.mat ${dwi_dir}/previousb0.mat
    fi

    # Now find all the rotations and translations from the current matrix
    # and save them in ec_rot.txt and ec_trans.txt
    $FSLDIR/bin/avscale --allparams ${dwi_dir}/current.mat ${dwi_dir}/$basenm | grep 'Rotation Angles' | sed 's/.* = //' >> ${dwi_dir}/ec_rot.txt ;
    $FSLDIR/bin/avscale --allparams ${dwi_dir}/current.mat ${dwi_dir}/$basenm | grep 'Translations' | sed 's/.* = //' >> ${dwi_dir}/ec_trans.txt ;

    # Finally, increase the counter and carry on the loop
    i=`echo $i + 1 | bc`
done

# clean up temp files
/bin/rm -f ${dwi_dir}/grot_labels.txt ${dwi_dir}/grot.oldmat ${dwi_dir}/grot.refmat ${dwi_dir}/grot.mat 
