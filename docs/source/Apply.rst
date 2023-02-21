************
Apply Your Own Images
************


Prepare to run your own images through the dolphot-lc pipeline
=============================

Jupyter Notebook Set Up
-----
Note that the following steps are also layed out in a Jupyter Notebook interface here: https://nbviewer.org/gist/whit5224/287af111f44bf83a23eaaf19a5121c75

In the command line be sure to activate your python=3.7, dolphot astroconda environment. ::

    conda activate /env name/

To use install this conda environment in Jupyter Notebook, run the following command. ::

    python -m ipykernel install --user --name=/env name/

Now navigate to Jupyter Notebook and change the kernel version of the notebook to your conda environment.


Python Test File Set Up
-----
::
    import dolphot_lc as dlc
    import time
    import os

    start = time.time()
    sexpath = 'source-extractor'

Use the example folder as the base directory. Within this directory create the folders "Images", "image_backup", and "ref". Within image_backup folder create the folders "ims", "ref", and "template". Within Images folder create the folders "ims" and "template".

Place fits image files in /example/Images/, /example/Images/ims/, and /example/image_backup/ims/.

Make sure your coadded image is named "coadd_CLASH_sci.fits". Place this file in /example/ref/, and /example/image_backup/ref/.

Make sure your coadded image is named "coadd_tweak_F110W_sci.fits". Place this file in /example/Images/, /example/Images/template, and /example/image_backup/template/.

(Note that the same file is used as the template and registration image.)

Now list these working directories and add your fits files to the correct locations. ::

    imroot = os.getcwd()                                # example folder / base directory
    orig_im_loc = f'{imroot}/image_backup/ims'          # image
    orig_temp_loc = f'{imroot}/image_backup/template'   # template
    im_loc = f'{imroot}/Images'                         # path to working directory
    ref_image = f'{imroot}/ref/coadd_CLASH_sci.fits'    # registration image
    dolphot_path = '/Users/rwhite/dolphot_lc/example/dolphot2.0/bin'

List the RA and DEC coordinates of the object in sexigesimal and degree units. ::

    sn_ra_me, sn_dec_me = '13:47:31.8180', '-11:45:51.914'
    objCoords = {'S1': [206.8833, -11.7644],
                 'S2': [206.8833, -11.7644],
                 'S3': [206.8833, -11.7644],
                 'S4': [206.8833, -11.7644],
                 'SX': [206.8833, -11.7644]
                 }

Input dolphot parameters.
(For help with this step, identify the HST instrument used to aquire your images and visit DOLPHOT user's guide for recommended settings. http://americano.dolphinsim.com/dolphot/dolphotWFC3.pdf) ::

    dolphot_params = {
        'UseWCS': 1,
        'raper': 5,
        'rchi': 1.5,
        'rsky0': 8,
        'rsky1': 3,
        'rpsf': 15,
        'WFC3IRpsfType': 1,      # wide field camera 3 IR
        }
        
Run Dolphot Process
-----
Nothing to change here. ::

    a = dlc.prep_directory(orig_im_loc, orig_temp_loc, im_loc, ref_image,
                           dolphot_path, imroot, sn_ra_me, sn_dec_me, sexpath,
                           dolphot_params)
    
    dlc.prep_files_for_dolphot('/dolphot_prepped',
                               r_in=15,
                               r_out=35,
                               step=4,
                               sig_low=2.25,
                               sig_high=2.00,
                               dlc_param=a)  #Processes images through masking/spliting/calcsky programs in Dolphot

    dlc.dolphot_simultaneous(a)         #Creates Dolphot parameter file & runs on processed images

    dlc.blot_back(r_in=15,              #Inner radius of sky annulus
                  r_out=35,             #Outer radius of sky annulus
                  step=4,               #How often is the sky value is sampled in pixels
                  sig_low=2.25,         #Low sigma under which samples will be rejected
                  sig_high=2.00,        #High sigma above which samples will be rejected
                  dlc_param=a)          #Parameter object from prep_directory function
            #Blots coadded template image to distorted science images & creates difference image
    
    dlc.dolphot_force(apermag=False, force_same_mag=True,  psfphot=1,
                      objCoords=objCoords, dlc_param=a)             #Runs Dolphot on difference images

    end = time.time()
    a = (end - start)/60
    m = 0
    while a > 1:
        a = a - 1
        m = m + 1
    print(f'Time: {m}:{str(int(60*a)).zfill(2)}')

Aftering running DOLPHOT forced photometry, the results are placed in the /example/diffs/ folder to be analyzed.
