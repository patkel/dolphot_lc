import dolphot_lc as dlc
import time
import os

def download(items):
    for thing in items:
        os.system(f'wget {thing}')
        os.rename(thing.split('/Download/')[1].replace('/','%2F'), thing.split('/')[-1])

items1 = ['https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/jcdu55piq_flc.fits',
          'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/jcdu55pkq_flc.fits',
          'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/jcdu55poq_flc.fits',
          'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/jcdu55prq_flc.fits']

items2 = ['https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/j8qu08y9q_flc.fits',
          'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/j8qu08ykq_flc.fits',
          'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/j8qu08ysq_flc.fits',
          'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/j8qu08ywq_flc.fits',
          'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/j8qu08z0q_flc.fits',
          'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/j8qu08z4q_flc.fits']

download_files = True

if download_files:
    os.mkdir('image_backup')
    os.mkdir('image_backup/ims')
    os.mkdir('image_backup/template')
    os.mkdir('image_backup/ref')

    os.chdir('image_backup/ims')
    download(items1)

    os.chdir('../template')
    download(items2)

    os.chdir('../../')


start = time.time()

imroot = os.getcwd()

orig_im_loc = f'{imroot}/image_backup/ims'
orig_temp_loc = f'{imroot}/image_backup/template'
im_loc = f'{imroot}/Images'
ref_image = f'{imroot}/ref/registration.fits'
dolphot_path = f'{imroot}/dolphot2.0/bin'

# sn_ra_me, sn_dec_me = '11:49:35.658', '+22:24:48.03'
sn_ra_me, sn_dec_me = '11:49:35.5', '+22:23:44'

sexpath = 'source-extractor'

dolphot_params = {
    'UseWCS': 1,
    'raper': 3,
    'rchi': 2.0,
    'rsky0': 15,
    'rsky1': 35,
    'rpsf': 10,
    'WFC3UVISpsfType': 1,
    }

objCoords = {'S1': [177.398231, 22.395630],
             'S2': [177.397718, 22.395786],
             'S3': [177.397370, 22.395538],
             'S4': [177.397810, 22.395189],
             'SX': [177.400112, 22.396701]
             }


a = dlc.prep_directory(orig_im_loc, orig_temp_loc, im_loc, ref_image,
                       dolphot_path, imroot, sn_ra_me, sn_dec_me, sexpath,
                       dolphot_params)

dlc.coadd_ims(a)

dlc.align_sci(a)

dlc.prep_files_for_dolphot('/dolphot_prepped',
                           r_in=15,
                           r_out=35,
                           step=4,
                           sig_low=2.25,
                           sig_high=2.00,
                           dlc_param=a)

dlc.dolphot_simultaneous(a)


dlc.blot_back(r_in=15,
              r_out=35,
              step=4,
              sig_low=2.25,
              sig_high=2.00,
              dlc_param=a)

dlc.dolphot_force(apermag=False, force_same_mag=True,  psfphot=1,
                  objCoords=objCoords, dlc_param=a)

end = time.time()

a = (end - start)/60
m = 0
while a > 1:
    a = a - 1
    m = m + 1
print(f'Time: {m}:{str(int(60*a)).zfill(2)}')

