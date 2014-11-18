#!/usr/bin/env python

'''
DIME_report.py

This code provides a succinct report of your participant's
motion during their diffusion tensor imaging scan so that
you can easily compare quality across participants.

You can access the help by typing:
    DIME_report.py --help
    
Inputs:
    dwi_file      4D diffusion weighted nifti file in nii.gz format
    bvals_file    text file containing b values in FSL format
                      ie: one entry per volume in one row
    bvecs_file    text file containing b vectors in FSL format
                      ie: one entry for each volume's x, y and z components
                          on three rows
    sub_id        text string containing the subject id
    
Created by Kirstie Whitaker on 17th November 2014
Contact: kw401@cam.ac.uk or HappyPenguin on Github
'''

#==============================================================================
# Import the modules you need
#==============================================================================
import argparse
from glob import glob
import os
import itertools as it
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pylab as plt
import nibabel as nib
import numpy as np
import pandas as pd
import shutil
import sys

#==============================================================================
# Define some useful functions
#==============================================================================
def setup_argparser():
    '''
    Setup argparser is the help function.
    You can access this by typing DIME_report.py --help
    '''
    
    #------------------------------------------------------
    # Build a basic parser.
    #------------------------------------------------------
    help_text = ('Create a quality control report for a diffusion weighted acquisition')
    
    sign_off = 'Author: Kirstie Whitaker <kw401@cam.ac.uk>'
    
    parser = argparse.ArgumentParser(description=help_text, epilog=sign_off)
    
    #------------------------------------------------------
    # Now add the arguments
    #------------------------------------------------------
    
    # Required argument: dwi_file
    parser.add_argument(dest='dwi_file', 
                            type=str,
                            metavar='dwi_file',
                            help='Full path to 4D diffusion weighted nifti image')
    
    # Required argument: bvals_file
    parser.add_argument(dest='bvals_file', 
                            type=str,
                            metavar='bvals_file',
                            help='Text file containing b values in "FSL" format - one entry per volume in one row')

    # Required argument: bvecs_file
    parser.add_argument(dest='bvecs_file', 
                            type=str,
                            metavar='bvecs_file',
                            help='Text file containing b vectors in "FSL" format - three entries per volume arranged on three rows')

    # Required argument: subid
    parser.add_argument(dest='sub_id', 
                            type=str,
                            metavar='sub_id',
                            help='Subject id')
                            
    arguments = parser.parse_args()
    
    return arguments, parser

#=============================================================================
def add_header(fig, grid, sub_id):
    '''
    The usefulness of the DIME report is that you can get lots of information
    in one place. So let's make a nice pretty header that contains the
    subject id, a space for the date, and check boxes for "PASS" and "FAIL"
    '''
    ax = plt.Subplot(fig, grid[0])
    fig.add_subplot(ax)

    # The header simply says:
    header_text = "Diffusion Imaging Motion Evaluation\n\nSubID: {}  Date:________________".format(sub_id[:15].ljust(15,'_'))
    
    ax.text(0.05, 0.5, header_text, transform=ax.transAxes, fontsize=14,
                   horizontalalignment='left', verticalalignment='center')
    
    # On the right we'll add two options:
    quality_text = "Pass"

    ax.text(0.82, 0.55, quality_text, transform=ax.transAxes, fontsize=18,
                   horizontalalignment='right', verticalalignment='bottom')

    ax.add_patch(patches.Rectangle((0.84,0.55),0.1,0.35, color='black', fill=False))

    quality_text = "Fail"

    ax.text(0.82, 0.15, quality_text, transform=ax.transAxes, fontsize=18,
                   horizontalalignment='right', verticalalignment='bottom')
    
    ax.add_patch(patches.Rectangle((0.84,0.15),0.1,0.35, color='black', fill=False))
    
    # Turn off axis labels
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_frame_on(False)

    return fig

#=============================================================================
def add_background(fig, grid):
    '''
    Make the background of the grid black so it looks nice :)
    '''
    ax = plt.Subplot(fig, grid[0])
    fig.add_subplot(ax)

    # Add a black background
    black = ax.imshow(np.ones([100,100]),
                            interpolation='none',
                            cmap='gray',
                            aspect='auto')
    # Turn off axis labels
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_frame_on(False)

    return fig
    
#==============================================================================
def plot_dti_slices(background_file, overlay_file, fig, grid, ax_name_list, cmap='jet'):
    '''
    plot_dti_slices places an image of an overlay file on top of the given
    background file and then returns it in a particular matplotlib grid
    location within a figure
    '''
    #------------------------------------------------------
    # Load the data (background and overlay files)
    #------------------------------------------------------
    bg_img = nib.load(background_file)
    bg = bg_img.get_data()
    
    # If the background file is 4D then only keep the first volume
    if len(bg.shape) == 4:
        bg = bg[:,:,:,0]
    
    overlay_img = nib.load(overlay_file)
    overlay = overlay_img.get_data()

    #------------------------------------------------------
    # Ensure the data is in float format
    #------------------------------------------------------
    bg = bg/1.
    overlay = overlay/1.
    
    #------------------------------------------------------
    # Scale the data by its maximum
    #------------------------------------------------------
    bg = bg / bg.max()
    overlay = overlay / overlay.max()    
        
    # Now we're going to loop through the different slice orientations
    # There is almost certainly a better way to do this. If you know 
    # of any tools that do this in just a few lines (probably in nipy?)
    # please do let Kirstie (HappyPenguin) know on GitHub
    # www.github.com/HappyPenguin/DIME

    for i, axis_name in enumerate(ax_name_list):
        if axis_name == 'axial':
            # Align so that right is right
            overlay_plot = np.rot90(overlay)
            overlay_plot = np.fliplr(overlay_plot)
            bg_plot = np.rot90(bg)
            bg_plot = np.fliplr(bg_plot)
    
        elif axis_name == 'coronal':
            overlay_plot = np.rot90(overlay)
            bg_plot = np.rot90(bg)
            overlay_plot = np.flipud(np.swapaxes(overlay_plot, 0, 2))
            bg_plot = np.flipud(np.swapaxes(bg_plot, 0, 2))

        elif axis_name == 'sagittal':
            overlay_plot = np.flipud(np.swapaxes(overlay, 0, 2))
            bg_plot = np.flipud(np.swapaxes(bg, 0, 2))

        # Calculate the number of images that you can fit on the page
        # n = height/width for bg_plot / height/0.15*width for figure
        # This gives a slight overlap of 0.15*width so you don't
        # waste too much space
        n = (( np.float(bg_plot.shape[1])/bg_plot.shape[2] ) 
                * (np.float(figsize[0])/ (0.15 * figsize[1])))

        # Round the number of images you can fit on the page to 
        # the lowest integer
        n_floor = np.int(np.floor(n))
        
        # Create a grid that contains spaces for the n images
        inner_grid = gridspec.GridSpecFromSubplotSpec(1, n_floor,
                         subplot_spec=grid[i], wspace=0.0, hspace=0.0)
        
        # Loop through the volume plotting n evenly spaced images:
        for j, slice_id in enumerate(np.linspace(0 , bg_plot.shape[2], n_floor+2)[1:-1]):
        
            # Define the slices you want to show
            bg_slice = bg_plot[:,:,slice_id]
            overlay_slice = overlay_plot[:,:,slice_id]
            
            ax = plt.Subplot(fig, inner_grid[j])
            fig.add_subplot(ax)

            # Add a black background
            black = ax.imshow(np.ones_like(bg_slice),
                                    interpolation='none',
                                    cmap='gray')
            
            # Mask the overlay data
            m_overlay_slice = np.ma.masked_where(overlay_slice==0, overlay_slice)

            # First show the background slice
            im1 = ax.imshow(bg_slice,
                                interpolation='none',
                                cmap='gray',
                                vmin = 0,
                                vmax = 1)
    
            # Then overlay the overlay_slice
            im2 = ax.imshow(m_overlay_slice,
                                interpolation='none',
                                cmap=cmap,
                                vmin = 0,
                                vmax = 1,
                                alpha = 0.3)
               
            # Turn off axis labels
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_frame_on(False)
            
    return fig

#=============================================================================
def plot_movement_params(dti_dir, fig, grid):
    measures = ['abs', 'rel']
    measure_suffixes = [ '', '_notb0' ]
    
    # Loop through the measure suffixes ('', '_notb0')
    # which you can think of as the groups of volumes that are being considered
    # and find the data from each of thoese files
    for i, suffix in enumerate(measure_suffixes):
        
        ax = plt.Subplot(fig, grid[i])
        fig.add_subplot(ax)

        # Read in the files
        disp = pd.read_csv(os.path.join(dti_dir, 'ec_disp{}.txt'.format(suffix)),
                            delimiter=' ', header=None,
                            names=['abs'+suffix, 'rel'+suffix], na_values='.')

        # Loop through the three different values that you want to know
        # for the two different measures (abs and rel)
        for measure in measures:

            ax.plot(disp[measure+suffix][disp[measure+suffix].notnull()], label=measure)
                
        # Label the x axis according to which plot this is:
        if suffix == '':
            ax.set_xlabel('Volume Number')
        else:
            ax.set_xlabel('Diff weighted Volume Number')
        
        # Set the y axis to always between 0 and 3
        ax.set_ylim(0,3)
        
        # Only label the first y axis
        if i == 0:
            # And label the yaxis
            ax.set_ylabel('Displacement (mm)')
            
        # Add a legend
        leg = ax.legend(loc=2, fontsize=8)
        leg.get_frame().set_alpha(0.5)
        
    return fig

#=============================================================================
def report_movement_params(dmri_motion_file, fig, movement_grid):
    '''
    Fill in a little text box with the TRACULA summary measures
    '''
    ax = plt.Subplot(fig, grid[2])
    fig.add_subplot(ax)
    
    motion = pd.read_csv(dmri_motion_file, delimiter=' ')

    # Insert some explanatory text
    header_text = 'TRACULA summary movement measures'

    ax.text(0.5, 0.92, header_text, transform=ax.transAxes, fontsize=18,
               horizontalalignment='center', verticalalignment='center')
    
    reference_text = [ u'For more details see:', 
    u'  Yendiki, A., Koldewyn, K., Kakunoori, S., Kanwisher, N., & Fischl, B. (2013).',
    u'  Spurious group differences due to head motion in a diffusion MRI study.',
    u'  NeuroImage, 88, 79–90. doi:10.1016/j.neuroimage.2013.11.027' ]
    
    ax.text(0.03, 0.15, '\n'.join(reference_text), transform=ax.transAxes, fontsize=11,
               horizontalalignment='left', verticalalignment='center')
               
    # Average Translation
    report_text = '{}: {:2.3f} mm'.format('Average Translation', np.float(motion['AvgTranslation']))
    ax.text(0.5, 0.78, report_text, transform=ax.transAxes, fontsize=14,
                            horizontalalignment='center', verticalalignment='center')

    # Average Rotation
    report_text = u'{}: {:2.3f}°'.format('Average Rotation', np.float(motion['AvgRotation']))
    ax.text(0.5, 0.63, report_text, transform=ax.transAxes, fontsize=14,
                            horizontalalignment='center', verticalalignment='center')
    
    # Portion of slices with dropout
    report_text = '{}: {:2.1f}%'.format('Proportion of slices with dropout'.ljust(35), np.float(motion['PercentBadSlices']))
    ax.text(0.5, 0.48, report_text, transform=ax.transAxes, fontsize=14,
                            horizontalalignment='center', verticalalignment='center')
                            
    # Dropout score
    report_text = '{}: {:2.1f}'.format('Dropout score', np.float(motion['AvgDropoutScore']))
    ax.text(0.5, 0.33, report_text, transform=ax.transAxes, fontsize=14,
                            horizontalalignment='center', verticalalignment='center')
    
    
    return fig

#=============================================================================
def tensor_histogram(fa_file, mo_file, sse_file, wm_mask_file, fig, grid):
    
    # Load in the data
    fa_img = nib.load(fa_file)
    fa = fa_img.get_data()
    mo_img = nib.load(mo_file)
    mo = mo_img.get_data()
    sse_img = nib.load(sse_file)
    sse = sse_img.get_data()
    
    wm_mask_img = nib.load(wm_mask_file)
    wm_mask = wm_mask_img.get_data()
    wm_mask[wm_mask>0] = 1
    
    # Mask the fa data with the white matter mask
    # so we're only looking inside the mask
    fa = fa * wm_mask
    
    # Add a subplot to the first space in the grid
    # and enter a histogram of FA values
    ax = plt.Subplot(fig, grid[0])
    fig.add_subplot(ax)    
    ax.hist(fa[fa>0].reshape(-1), bins=np.linspace(0,1,100), color='green',histtype='stepfilled')
    # Label the x axis:
    ax.set_xlabel('Fractional Anisotropy')
    # Set the y axis to always between 0 and 2500
    ax.set_ylim(0,2500)
    # Adjust the power limits so that you use scientific notation on the y axis
    ax.ticklabel_format(style='sci', axis='y')
    ax.yaxis.major.formatter.set_powerlimits((-3,3))

    # Only label this first y axis as they're all the same
    ax.set_ylabel('Number of voxels')
    
    # Add a subplot to the second space in the grid
    # and plot a histogram of mode values
    ax = plt.Subplot(fig, grid[1])
    fig.add_subplot(ax)    
    ax.hist(mo[fa>0].reshape(-1), bins=np.linspace(-1,1,100), color='orange', histtype='stepfilled')
    # Label the x axes:
    ax.set_xlabel('Mode of Anisotropy')
    # Set the y axis to always between 0 and 3500
    ax.set_ylim(0,3500)
    # Adjust the power limits so that you use scientific notation on the y axis
    #plt.ticklabel_format(style='sci', axis='y')
    ax.yaxis.major.formatter.set_powerlimits((-3,3))

    
    # Add a subplot to the third space in the grid
    # and plot a histogram of sum of square errors
    # Note that low values are very good - they indicate voxels
    # that have a good fit to the diffusion tensor model.
    # The y-axis is therefore limited so that the histogram highlights
    # "bad" fit voxels.

    ax = plt.Subplot(fig, grid[2])
    fig.add_subplot(ax)    
    ax.hist(sse[fa>0].reshape(-1), bins=np.linspace(0,5,100), color='red', histtype='stepfilled')
    # Label the x axis:
    ax.set_xlabel('Sum of Square Errors')
    # Set the y axis to always between 0 and 3500
    ax.set_ylim(0,100)
    
    return fig







#=============================================================================
# Define the files you're going to need or create
#=============================================================================
# Read in the arguments from argparse
arguments, parser = setup_argparser()
dwi_file = arguments.dwi_file
bvals_file = arguments.bvals_file
bvecs_orig_file = arguments.bvecs_file
sub_id = arguments.sub_id

# Figure out the name of the directory in which the dwi file is saved
dwi_dir = os.path.dirname(arguments.dwi_file)

# Figure out the name of the DIME scripts directory
dime_scripts_dir = os.path.dirname(os.path.realpath(sys.argv[0]))

# Define all the filenames that you're going to create
dwi_ec_file = dwi_file.replace('.nii.gz', '_ec.nii.gz')
dwi_ec_brain_file = dwi_ec_file.replace('.nii.gz', '_brain.nii.gz')
dwi_ec_brain_mask_file = dwi_ec_file.replace('.nii.gz', '_brain_mask.nii.gz')
dmri_motion_file = os.path.join(dwi_dir, 'dmri_motion_output.txt')
ecclog_file = dwi_ec_file.replace('.nii.gz', '.ecclog')
shutil.move(bvecs_orig_file, os.path.join(dwi_dir, 'bvecs_orig'))
bvecs_rot_file = os.path.join(dwi_dir, 'bvecs_rotated')

fdt_root = os.path.join(dwi_dir, 'FDT', sub_id)
fa_file = fdt_root + '_FA.nii.gz'
fa_ero_file = fdt_root + '_FA_ero.nii.gz'
mo_file = fdt_root + '_MO.nii.gz'
sse_file = fdt_root + '_sse.nii.gz'

# Set up the log files
output_log_file = os.path.join(dwi_dir, 'DIME_report_log_output')
error_log_file = os.path.join(dwi_dir, 'DIME_report_log_error')

# For now delete them every time....
# I might change this in the future..
if os.path.exists(output_log_file):
    os.remove(output_log_file)
if os.path.exists(error_log_file):
    os.remove(error_log_file)

#=============================================================================
# Run FSL's eddy correct
#=============================================================================
command = 'eddy_correct {} {} 0 > {} 2> {}'.format(dwi_file,
                                                    dwi_ec_file,
                                                    output_log_file, 
                                                    error_log_file)

if not os.path.isfile(dwi_ec_file):
    print '    Running eddy_correct....'
    os.system(command)

print '      Eddy_correct complete'

#=============================================================================
# Rotate the bvecs file as a result of the rotations applied by eddy_correct
#=============================================================================
command = 'xfmrot {} {} {}'.format(ecclog_file, 
                                    bvecs_orig_file, 
                                    bvecs_rot_file)

if not os.path.isfile(bvecs_rot_file):
    print '    Rotating bvecs....'
    os.system(command)

print '      Bvecs rotated'

#=============================================================================
# Extract motion values as calculated by eddy_correct
#=============================================================================
# NOTE That this can probably be wrapped into the dmri_motion
# command - there's just a hidden flag that I don't know about....
motion_calc_script = os.path.join(dime_scripts_dir, 'dti_ec_motion_calc.sh')
command = '{} {} {} > {} 2> {}'.format(motion_calc_script,
                                            ecclog_file,
                                            bvals_file,
                                            output_log_file, 
                                            error_log_file)
                                            
if not os.path.isfile(os.path.join(dwi_dir, 'ec_rot.txt')):
    print '    Extracting FSL motion measures....'
    os.system(command)

print '      FSL motion measures extracted'

#=============================================================================
# Run the dmri_motion command that is part of Freesurfer's TRACULA
#=============================================================================
command = ( 'dmri_motion --dwi {} '
            '--bval {} --mat {} '
            '--out {} > {} 2> {}'.format(dwi_file, 
                                            bvals_file, 
                                            ecclog_file, 
                                            dmri_motion_file,
                                            output_log_file, 
                                            error_log_file) )

if not os.path.isfile(dmri_motion_file):
    print '    Extracting TRACULA motion parameters....'
    os.system(command)

print '      TRACULA motion parameters extracted'

#=============================================================================
# Run FSL's brain extraction command
#=============================================================================
command = 'bet {} {} -f 0.15 -m > {} 2> {}'.format(dwi_ec_file, 
                                                    dwi_ec_brain_file,
                                                    output_log_file, 
                                                    error_log_file)


if not os.path.isfile(dwi_ec_brain_file):
    print '    Running brain extraction....'
    os.system(command)

print '      Brain extraction complete'

#=============================================================================
# Fit the diffusion tensor by running FSL's FDT command
#=============================================================================
command = ( 'dtifit -k {} -m {} '
            '-b {} -r {} --sse '
            '--save_tensor -o {} '
            ' > {} 2> {}'.format(dwi_ec_file, 
                                    dwi_ec_brain_mask_file,
                                    bvals_file, 
                                    bvecs_rot_file,
                                    fdt_root,
                                    output_log_file, 
                                    error_log_file) )

if not os.path.isfile(fa_file):
    print '    Fitting tensor....'
    try:
        os.makedirs(os.path.dirname(fdt_root))
    except OSError as exception:
        pass
        
    os.system(command)

print '      Tensor fit complete'

#=============================================================================
# Create a rough white matter mask based on the FA image
#=============================================================================
command = ( 'fslmaths {} -ero -ero -thr 0.1 -bin {}'.format(fa_file,
                                                             fa_ero_file) )

if not os.path.isfile(fa_ero_file):
    print '    Creating white matter mask by eroding & thresholding FA image'
    os.system(command)
    
print '      White matter mask created'

#=============================================================================
# Create a figure that's the same size as an A4 piece of paper
# (sorry to folk who use letter paper as standard!)
#=============================================================================

figsize = (8.3,11.6)
fig = plt.figure(figsize = figsize)

#=============================================================================
# Set up the various locations within this figure that you're going to fill
#=============================================================================
# Header grid - to contain the header at the top of the page
header_grid = gridspec.GridSpec(1,1)
header_grid.update(left=0.05, right=0.95, top = 0.98, bottom = 0.9)

# Background A - to go behind brain_grid A
bgA_grid = gridspec.GridSpec(1, 1)
bgA_grid.update(left=0.05, right=0.95, top = 0.9, bottom = 0.65)

# Brain grid A - non-diffusion weighted image and brain mask
brainA_grid = gridspec.GridSpec(3, 1)
brainA_grid.update(left=0.05, right=0.95, top = 0.9, bottom = 0.65)

# Movement grid - movement and eddy_correct realignment parameters 
movement_grid = gridspec.GridSpec(1, 3)
movement_grid.update(left=0.1, right=0.95, top = 0.63, bottom = 0.5, wspace=0.2)

# Background B - to go behind brain_grid B
bgB_grid = gridspec.GridSpec(1, 1)
bgB_grid.update(left=0.05, right=0.95, top = 0.45, bottom = 0.2)

# Brain grid B - FA image and white matter mask
brainB_grid = gridspec.GridSpec(3, 1)
brainB_grid.update(left=0.05, right=0.95, top = 0.45, bottom = 0.2)

# Histogram grid - histograms of FA, MO, and sum of square errors
hist_grid = gridspec.GridSpec(1, 3)
hist_grid.update(left=0.1, right=0.95, top = 0.18, bottom = 0.05, wspace=0.2)


#=============================================================================
# Fill in these plotting areas using the functions defined above
#=============================================================================

fig = add_header(fig, header_grid, sub_id)
fig = add_background(fig, bgA_grid)
fig = plot_dti_slices(dwi_ec_brain_file, dwi_ec_brain_mask_file, fig, brainA_grid, ['sagittal', 'coronal', 'axial'], cmap='cool_r')
fig = plot_movement_params(dwi_dir, fig, movement_grid)
fig = report_movement_params(dmri_motion_file, fig, movement_grid)
fig = add_background(fig, bgB_grid)
fig = plot_dti_slices(fa_file, fa_ero_file, fig, brainB_grid, ['sagittal', 'coronal', 'axial'], cmap='cool')
fig = tensor_histogram(fa_file, mo_file, sse_file, fa_ero_file, fig, hist_grid)


#=============================================================================
# Finally, save the figure
#=============================================================================

report_filename = os.path.join(dwi_dir, 'DIME_report.jpg')
fig.savefig(report_filename, bbox_inches=0, dpi=300)

#=============================================================================
# **That's it! You're done :)**
#=============================================================================

#=============================================================================
# If you have any feedback please add it to the GitHub page
# www.github.com/HappyPenguin/DIME

# I wish you low motion in your data for years to come
# Kx
#=============================================================================
