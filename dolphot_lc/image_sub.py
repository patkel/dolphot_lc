import os
import shutil
from .dolphot_scrips import prep_files_for_dolphot
from astropy.io import fits
from stsci.tools import teal
from drizzlepac.ablot import blot


def blot_back(r_in, r_out, step, sig_low, sig_high, dlc_param):
    '''
    Blots coadded template image to distorted science images and
    creates difference image

    Parameters
    ----------

    r_in : float
        Inner radius of sky annulus

    r_out : float
        Outer radius of sky annulus

    step : float
        How often is the sky value is sampled in pixels

    sig_low : float
        Low sigma under which samples will be rejected

    sig_high : float
        High sigma above which samples will be rejected

    dlc_param : obj
        Parameter object from prep_directory function
    '''

    _prep_imaging(dlc_param)
    im_drz_blot = f'{dlc_param.IMROOT}/Images/coadd_tweak_'\
                  f'{dlc_param.FILT}_sci.fits'

    diff_dir = f'{dlc_param.IMROOT}/diffs/'
    dolphot_prepped_dir = f'{dlc_param.IMROOT}/dolphot_prepped/'

    rescale_fac = 1

    try:
        os.mkdir(diff_dir)
    except FileExistsError:
        pass

    for im in dlc_param.IMAGES:  # all_images:
        im_to_blot = f'{dlc_param.IMROOT}/imaging/{im.name}_'\
                     f'{dlc_param.SUFFIX}.fits'
        im_to_blot_dolphot_prepped = f'{dlc_param.IMROOT}/dolphot_prepped/'\
                                     f'{im.name}_{dlc_param.SUFFIX}.fits'

        p = fits.open(im_to_blot)
        p_dol_prep = fits.open(im_to_blot_dolphot_prepped)

        im_diff = f'{diff_dir}{im.name}_{dlc_param.SUFFIX}.fits'

        for chip in dlc_param.CHIPS:
            outdata = f'{dlc_param.IMROOT}/imaging/{im.name}_'\
                      f'{dlc_param.SUFFIX}_bgblot_{chip:d}.fits'

            EXPTIME_DRZ = _gethead(im_drz_blot, 'EXPTIME')

            EXPTIME_BLT = _gethead(im_to_blot, 'EXPTIME')

            try:
                os.remove(outdata)
            except FileNotFoundError:
                pass

            blotobj = teal.load('ablot')

            if dlc_param.INST == 'ACS' or (dlc_param.INST == 'WFC3'
                                           and dlc_param.DETEC == 'UVIS'):
                blot(im_drz_blot, f'{im_to_blot}[sci,{chip:d}]', outdata,
                     addsky=False, in_units='cps', out_units='counts',
                     expout=1./EXPTIME_DRZ*EXPTIME_BLT*rescale_fac,
                     configObj=blotobj)

            elif dlc_param.INST == 'WFC3' and dlc_param.DETEC == 'IR':
                blot(im_drz_blot, f'{im_to_blot}[sci,{chip:d}]', outdata,
                     addsky=False, in_units='cps', out_units='counts',
                     expout=1./EXPTIME_DRZ*rescale_fac, configObj=blotobj)

            conv = p_dol_prep[1].data / p[1].data

            a = fits.open(outdata)

            chip_name = f'{im.name}_{dlc_param.SUFFIX.lower()}'\
                        f'.chip{chip:d}.fits'

            if chip == 1:
                p[1].data = 1. * (p[1].data - a[1].data)  # * conv

            elif chip == 2:
                p[4].data = 1. * (p[4].data - a[1].data)  # * conv

            shutil.copyfile(f'{dolphot_prepped_dir}'
                            f'{chip_name.replace(".fits", ".sky.fits")}',
                            f'{diff_dir}'
                            f'{chip_name.replace(".fits", ".sky.fits")}')

        p.writeto(im_diff, overwrite=True)

    prep_files_for_dolphot('/diffs', r_in, r_out,
                           step, sig_low, sig_high, dlc_param)


def _prep_imaging(dlc_param):
    try:
        os.mkdir(f'{dlc_param.IMROOT}/imaging')
    except FileExistsError:
        pass

    all_images = [a.loc for a in dlc_param.IMAGES]

    for im in all_images:
        shutil.copyfile(im, f'{dlc_param.IMROOT}/imaging/{im.split("/")[-1]}')


def _gethead(im, item):
    hdulist = fits.open(im)
    value = hdulist[0].header[item]
    hdulist.close()
    return value
