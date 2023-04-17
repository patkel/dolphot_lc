"""Test dlc_params.py | dolphot_scripts.py | image_sub.py | process_ims.py modules."""

import dolphot_lc as dlc
import os
from astropy.io import fits
from astropy.wcs import WCS


######################################################################
##### example required to run tests #######
imroot = os.getcwd()	#example folder / base directory

orig_im_loc = f'{imroot}/image_backup/ims'          # image (WFC3 instrume)
orig_temp_loc = f'{imroot}/image_backup/template'   # template
im_loc = f'{imroot}/Images'				# path to working directory
ref_image = f'{imroot}/ref/coadd_CLASH_sci.fits'       	#registration image
dolphot_path = f'{imroot}/dolphot2.0/bin'

sn_ra_me, sn_dec_me = '13:47:31.8180', '-11:45:51.914'

sexpath = 'source-extractor'

dolphot_params = {
    'UseWCS': 1,
    'raper': 5,
    'rchi': 1.5,
    'rsky0': 8,
    'rsky1': 3,
    'rpsf': 15,
    'WFC3IRpsfType': 1,      #wide field camera 3 IR
    }

objCoords = {'S1': [206.8826, -11.7644],
             'S2': [206.8826, -11.7644],
             'S3': [206.8826, -11.7644],
             'S4': [206.8826, -11.7644],
             'SX': [206.8826, -11.7644]
             }

# Run Prep Directory - copies raw images into working directory and creates dolphot object to be passed
a = dlc.prep_directory(orig_im_loc, orig_temp_loc, im_loc, ref_image,
                       dolphot_path, imroot, sn_ra_me, sn_dec_me, sexpath,
                       dolphot_params)

######################################################################

def test_prep_files_for_dolphot():
    """Test prep_files_for_dolphot, processes images through programs in Dolphot."""
    files = os.listdir(a.ORIG_IM_LOC)
    dlc.prep_files_for_dolphot('/dolphot_prepped', r_in=15, r_out=35, step=4, sig_low=2.25, sig_high=2.00, dlc_param=a)
    for entry in files:
        assert os.path.exists(f'{a.IMROOT}/Images/{entry}')

def test_dolphot_simultaneous():
    """Test dolphot_simultaneous, Creates Dolphot parameter file & runs on processed images.
        Ensure key parameters were correctly written using _mk_param function    """
    dlc.dolphot_simultaneous(a)

    #check keys in dolphot.param file properly updated
    for key in a.DOLPHOT_PARAMS:
        assert str(a.DOLPHOT_PARAMS[key]) == str(dolphot_params[key])

def test_dolphot_blot_back():
    """Test blot_back, blots coadded template image to distorted science images & creates difference image."""
    dlc.blot_back(r_in=15, 		        #Inner radius of sky annulus
                  r_out=35, 	        #Outer radius of sky annulus
                  step=4, 		        #How often is the sky value is sampled in pixels
                  sig_low=2.25, 	    #Low sigma under which samples will be rejected
                  sig_high=2.00, 	    #High sigma above which samples will be rejected
                  dlc_param=a) 	        #Parameter object from prep_directory function
    for im in a.IMAGES:  # all_images:
        diff_image = f'{im.name}_{a.SUFFIX}.fits'
        assert os.path.exists(f'{a.IMROOT}/diffs/{diff_image}')

def test_dolphot_force():
    """Test dolphot_force, runs Dolphot on difference images."""
    dlc.dolphot_force(apermag=False, force_same_mag=True, psfphot=1,
                      objCoords=objCoords, dlc_param=a)
    param_file = f'{a.IMROOT}/diffs/xytfile'
    f = open(param_file, 'r')
    for key in objCoords.keys():
        ra, dec = objCoords[key]
        hdu = fits.open(f'{ref_image}')[0]
        w = WCS(hdu.header)
        x, y = w.wcs_world2pix(ra, dec, 1)
        line = f.readline()
        assert str(line.strip()) == str(f'0 1 {x} {y} 2 10')

