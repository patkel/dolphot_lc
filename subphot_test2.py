import Reset_Dir
imroot = '/media/thomas/Drive_N1/HST_Icarus'
Reset_Dir.reset(imroot)

import subphot2 as sp

#import sys
#sys.stdout = open('output_log', 'w')

im_loc = '/media/thomas/Drive_N1/HST_Icarus/Images'
ref_image = '/media/thomas/Drive_N1/HST_Icarus/Images/IBF5B1DMQ_flt.fits'
dolphot_path = '/media/thomas/Drive_N1/dolphot2.0/bin'
sn_ra_me, sn_dec_me = '11:49:35.658', '+22:23:48.03'
imroot = '/media/thomas/Drive_N1/HST_Icarus'
sexpath = '/usr/local/bin/sex'

"""
sp.prep_dir(im_loc, ref_image, dolphot_path, imroot, sexpath, sn_ra_me, sn_dec_me)
sp.dolphot_file_prep(15, 35, 4, 2.25, 2.00)
sp.mk_diff(10, 'rscale', 3, 15)
#sp.dolphot_simultaneous()
sp.blot_back()
sp.make_difference_ims(90, objname='refsdal', big=True, single_sci=True)
"""

sp.prep_dir(im_loc, ref_image, dolphot_path, imroot, sexpath, sn_ra_me, sn_dec_me)

sp.mk_diff(threshold = 10,
    fitgeometry = 'rscale', 
    nclip = 3,
    minobj = 15)

sp.prep_files_for_dolphot('/dolphot_prepped',
    r_in = 15,
    r_out = 35,
    step = 4,
    sig_low = 2.25,
    sig_high = 2.00)
    
sp.dolphot_simultaneous()

sp.blot_back('/diffs',
    r_in = 15,
    r_out = 35,
    step = 4,
    sig_low = 2.25,
    sig_high = 2.00)

special = False

sp.dolphot_force(special=special, apermag=False, force_same_mag=False, psfphot=1)            
sp.dolphot_force(special=special, apermag=False, force_same_mag=True, psfphot=1)            
sp.dolphot_force(special=special, apermag=False, force_same_mag=False, psfphot=2)            
sp.dolphot_force(special=special, apermag=False, force_same_mag=True, psfphot=2)

print('Done Dude')

#sys.stdout.close()
