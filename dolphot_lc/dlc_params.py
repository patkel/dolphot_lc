import os
import shutil
from astropy.io import fits
from drizzlepac import skytopix


class _dlc_parameters:
    def __init__(self, ORIG_IM_LOC, ORIG_TEMP_LOC, IM_LOC, REF_IMAGE_FULL,
                 DOLPHOT_PATH, IMROOT, SN_RA_ME, SN_DEC_ME,
                 SEXPATH, DOLPHOT_PARAMS):
        self.ORIG_IM_LOC = ORIG_IM_LOC
        self.ORIG_TEMP_LOC = ORIG_TEMP_LOC
        self.IM_LOC = IM_LOC
        self.REF_IMAGE = os.path.split(REF_IMAGE_FULL)[1]
        self.REF_IMAGE_PATH = os.path.split(REF_IMAGE_FULL)[0]
        self.DOLPHOT_PATH = DOLPHOT_PATH
        self.IMROOT = IMROOT
        self.SN_RA_ME = SN_RA_ME
        self.SN_DEC_ME = SN_DEC_ME
        nodp = "_no_dolphot.fits"
        self.REF_IMAGE_NO_DOLPHOT = f'{self.REF_IMAGE_PATH}/'\
                                    f'{self.REF_IMAGE.replace(".fits", nodp)}'
        self.SEXPATH = SEXPATH
        self.IMAGES = _glob_image(self.ORIG_IM_LOC, f'{self.IM_LOC}/ims',
                                  self.REF_IMAGE_PATH, self.REF_IMAGE,
                                  self.IMROOT)
        self.TEMPLATE = _glob_image(self.ORIG_TEMP_LOC,
                                    f'{self.IM_LOC}/template',
                                    self.REF_IMAGE_PATH, self.REF_IMAGE,
                                    self.IMROOT)
        self.DOLPHOT_PARAMS = DOLPHOT_PARAMS
        self.INST = self.IMAGES[0].instrument
        self.DETEC = self.IMAGES[0].detector
        self.FILT = self.IMAGES[0].filter

        self.IMTYPE = 'fullarray'

        if self.INST == 'WFPC2':
            self.MASK = '/wfpc2mask'
            self.CHIPS = [1, 2, 3, 4]
            self.SUFFIX = 'c0m'

        if self.INST == 'WFC3':
            self.MASK = '/wfc3mask'
            if self.DETEC == 'UVIS':
                hdulist = fits.open(self.IMAGES[0].loc)
                if hdulist[1].header['NAXIS1'] < 3000. and \
                        hdulist[1].header['NAXIS2'] < 3000.:
                    self.CHIPS = [1]
                    self.IMTYPE = 'subarray'
                else:
                    self.CHIPS = [1, 2]
                self.SUFFIX = 'flc'

            if self.DETEC == 'IR':
                self.CHIPS = [1]
                self.SUFFIX = 'flt'

        if self.INST == 'ACS':
            self.MASK = '/acsmask'
            if self.DETEC == 'WFC':
                self.CHIPS = [1, 2]
                self.SUFFIX = 'flc'

        self.DR_SUFFIX = 'drz'
        if self.SUFFIX == 'flc':
            self.DR_SUFFIX = 'drc'

        self.ON_CHIP = _which_chip(self.IMAGES[0].loc, self.CHIPS,
                                   self.SN_RA_ME, self.SN_DEC_ME)


def _which_chip(im, chips, ra, dec):
    good_chip = None
    a = fits.open(im)
    for chip in chips:
        x, y = skytopix.rd2xy(f'{im}[sci,{chip}]', ra, dec)
        b = a[chip].header
        x_max = b['NAXIS1']
        y_max = b['NAXIS2']
        if x > 0 and x < x_max and y > 0 and y < y_max:
            good_chip = chip

    if good_chip is None:
        raise ValueError('Transient not in image')

    return good_chip


# Image object with useful properties
class _Image:
    def __init__(self, loc, instrument, detector, filter, typ, imroot):
        self.loc = loc
        self.instrument = instrument
        self.detector = detector
        self.filter = filter
        self.typ = typ
        self.prep_loc = f'{imroot}/dolphot_prepped/{loc.split("/")[-1]}'
        self.name = loc.split('/')[-1].split('.fits')[0].split('_')[0]


# Build IMAGES list with image objects
def _glob_image(orig_im_loc, im_loc, ref_image_path, ref_image, imroot):
    image_list = os.listdir(orig_im_loc)
    N = len(image_list)
    for i in range(0, N):
        image_list[i] = im_loc+'/'+image_list[i]

    image_details = [0]*N

    for i in range(0, N):
        im = fits.open(image_list[i])[0].header
        inst = im['INSTRUME']
        detec = im['DETECTOR']

        if inst == 'WFC3':
            filt = im['FILTER']
        if inst == 'ACS':
            filt = im['FILTER2']

        if image_list[i] == f'{ref_image_path}/{ref_image}':
            typ = 'ref'
        else:
            typ = 'sci'

        image_details[i] = _Image(image_list[i], inst, detec,
                                  filt, typ, imroot)

    return(image_details)


def _Images_Setup(im_loc, orig_im_loc, orig_temp_loc):

    try:
        os.mkdir(im_loc)
    except FileExistsError:
        pass

    try:
        os.mkdir(f'{im_loc}/ims')
        a = os.listdir(orig_im_loc)

        b = []

        for x in range(0, len(a)):
            if a[x][-4:] == 'fits':
                b.append(a[x])

        for x in range(0, len(b)):
            shutil.copyfile(f'{orig_im_loc}/{b[x]}', f'{im_loc}/ims/{b[x]}')

    except FileExistsError:
        pass

    try:
        os.mkdir(f'{im_loc}/template')
        a = os.listdir(orig_temp_loc)

        b = []

        for x in range(0, len(a)):
            if a[x][-4:] == 'fits':
                b.append(a[x])

        for x in range(0, len(b)):
            shutil.copyfile(f'{orig_temp_loc}/{b[x]}',
                            f'{im_loc}/template/{b[x]}')

    except FileExistsError:
        pass


def prep_directory(ORIG_IM_LOC, ORIG_TEMP_LOC, IM_LOC, REF_IMAGE_FULL,
                   DOLPHOT_PATH, IMROOT, SN_RA_ME, SN_DEC_ME,
                   SEXPATH, DOLPHOT_PARAMS):
    '''
    Copies raw images into a working directory and creates an object to be
    passed to other functions

    Parameters
    ----------
    ORIG_IM_LOC : str
        Path to original science images

    ORIG_TEMP_LOC : str
        Path to original template images

    IM_LOC : str
        Path of working directory

    REF_IMAGE_FULL : str
        Path to reference image

    DOLPHOT_PATH : str
        Path to Dolphot bin directory

    IMROOT : str
        Base directory where program will run

    SN_RA_ME : str
        Right ascension of transient object (Ex: 12:23:56.7)

    SN_DEC_ME : str
        Declination of transient object (Ex: -12:23:56.7)

    SEXPATH : str
        Path to Source-Extractor executable file

    DOLPHOT_PARAMS : dict
        Dictionary of Dolphot parameters and their values

    Returns
    -------
    param_obj : obj
        Parameter object to be passed to other functions
    '''

    _Images_Setup(IM_LOC, ORIG_IM_LOC, ORIG_TEMP_LOC)

    param_obj = _dlc_parameters(ORIG_IM_LOC, ORIG_TEMP_LOC, IM_LOC,
                                REF_IMAGE_FULL, DOLPHOT_PATH, IMROOT, SN_RA_ME,
                                SN_DEC_ME, SEXPATH, DOLPHOT_PARAMS)

    return param_obj
