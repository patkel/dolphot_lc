import os
import subprocess
import shutil
from glob import glob
from astropy.io import fits
from astropy import units as u
import astropy.coordinates as coord
from astropy.wcs import WCS


def dolphot_simultaneous(dlc_param):
    '''
    Creates Dolphot parameter file and runs Dolphot on processed images

    Parameters
    ----------

    dlc_param : obj
        Parameter object from prep_directory function
    '''

    imdir_simultaneous = f'{dlc_param.IMROOT}/dolphot_prepped'
    try:
        os.mkdir(imdir_simultaneous)
    except FileExistsError:
        pass

    ref_image_use = f'{dlc_param.REF_IMAGE}'

    ''' now set up dolphot parameter files '''
    imgNum = 0
    extra_params = {}
    info_params = {}

    files = [a.loc for a in dlc_param.IMAGES]

    imdir = f'{dlc_param.IMROOT}/dolphot/'
    imdir_dolphot_prepped = f'{dlc_param.IMROOT}/dolphot_prepped/'

    try:
        os.mkdir(imdir)
    except FileExistsError:
        pass

    for file in files:
        for chip in dlc_param.CHIPS:

            if imgNum >= 99:
                raise Exception

            name_root = file.split("/")[-1].split("_")[0]

            fname_prepped = f'{imdir_dolphot_prepped}{name_root}_'\
                            f'{dlc_param.SUFFIX}.chip{str(chip)}'

            if chip == dlc_param.ON_CHIP:
                imgNum += 1

                fname_simultaneous = f'{imdir_simultaneous}/{name_root}_'\
                                     f'{dlc_param.SUFFIX}.chip{chip:d}'

                if not glob(f'{fname_simultaneous}.fits'):
                    shutil.copy(f'{fname_prepped}.fits',
                                f'{fname_simultaneous}.fits')

                extra_params[f'img{imgNum}_file'] = f'{name_root}_'\
                                                    f'{dlc_param.SUFFIX}'\
                                                    f'.chip{str(chip)}'
                extra_params[f'img{imgNum}_shift'] = '0 0'
                extra_params[f'img{imgNum}_xform'] = '1 0 0'

                if dlc_param.INST != 'WFPC2':
                    info_params[f'img{imgNum}_instrument'] = dlc_param.INST
                    info_params[f'img{imgNum}_detector'] = dlc_param.DETEC
                    info_params[f'img{imgNum}_filt'] = dlc_param.FILT

                    orig = f'{dlc_param.IMROOT}/imaging/{file}_'\
                           f'{dlc_param.SUFFIX}.fits'
                    masked = f'{dlc_param.IMROOT}/dolphot/{file}_'\
                             f'{dlc_param.SUFFIX}.fits'

                    info_params[f'img{imgNum}_orig'] = orig
                    info_params[f'img{imgNum}_masked'] = masked

                if dlc_param.INST == 'WFPC2':
                    command = f'gethead {fname_simultaneous}.fits EXPNAME'
                    namef = subprocess.getoutput(command)

                    orig_crclean = f'{dlc_param.IMROOT}/imaging/'\
                                   f'{namef}_C0M_crclean.fits'
                    orig = f'{dlc_param.IMROOT}/imaging/{namef}_C0M.fits'
                    orig_dq = f'{dlc_param.IMROOT}/imaging/{namef}_C1M.fits'
                    fn_dp_masked = f'{dlc_param.IMROOT}/dolphot/'\
                                   f'{namef}_C0M.fits'

                    info_params[f'img{imgNum}_orig_crclean'] = orig_crclean
                    info_params[f'img{imgNum}_orig'] = orig
                    info_params[f'img{imgNum}_orig_dq'] = orig_dq
                    info_params[f'img{imgNum}_dolphot_masked'] = fn_dp_masked

                    info_params[f'img{imgNum}_instrument'] = dlc_param.INST
                    info_params[f'img{imgNum}_detector'] = dlc_param.DETEC
                    info_params[f'img{imgNum}_filt'] = dlc_param.FILT

    extra_params['Nimg'] = imgNum

    if dlc_param.IMTYPE == 'subarray':
        output_dir = f'{dlc_param.IMROOT}/coadd/'
        ref_image_use = f'{output_dir}{dlc_param.FILT}glass_drz.fits'

    if dlc_param.IMTYPE == 'subarray':
        if not glob(f'{imdir_simultaneous}{ref_image_use.split("/")[-1]}'):

            shutil.copy(ref_image_use, imdir_simultaneous)

            os.chdir(imdir_simultaneous)

            cmd = f'{dlc_param.DOLPHOT_PATH}wfc3mask '\
                  f'{ref_image_use.split("/")[-1]}'
            os.system(cmd)
            cmd = f'{dlc_param.DOLPHOT_PATH}splitgroups '\
                  f'{ref_image_use.split("/")[-1]}'
            os.system(cmd)
            for chip in [1]:
                cmd = f'{dlc_param.DOLPHOT_PATH}calcsky '\
                      f'{ref_image_use.split("/")[-1].replace(".fits","")}'\
                      f'.chip{chip} 15 35 4 2.25 2.00'
                os.system(cmd)

    else:
        refim_c1 = ref_image_use.replace(".fits", ".chip1.fits").split("/")[-1]
        if not glob(f'{imdir_simultaneous}{refim_c1}'):
            pass

    if True:
        extra_params['img0_file'] = ref_image_use.replace('.fits', '.chip1')
        extra_params['img0_RAper'] = '4'
        extra_params['img0_RChi'] = '2.0'
        extra_params['img0_RSky'] = '15 35'
        extra_params['img0_RPSF'] = '15'

    shutil.copyfile(f'{dlc_param.REF_IMAGE_PATH}/'
                    f'{ref_image_use.replace(".fits",".chip1.fits")}',
                    f'{imdir_simultaneous}/'
                    f'{ref_image_use.replace(".fits",".chip1.fits")}')

    param_file = f'{imdir_simultaneous}/dolphot.params'

    os.chdir(imdir_simultaneous)

    '''
    here using the recommended settings for a WFC3 UVIS registration image
    '''

    _mk_param(dlc_param.INST, dlc_param.DETEC, param_file, extra_params,
              dlc_param.IMTYPE, dlc_param)

    cmd = f'{dlc_param.DOLPHOT_PATH}/dolphot output -pdolphot.params'
    os.system(cmd)


def _mk_param(instrument, detector, param_file,
              extra_params, imtype, dlc_param):
    f = open(param_file, 'w')
    from copy import copy

    if imtype == 'subarray':
        dlc_param.DOLPHOT_PARAMS['UseWCS'] = 2

    for key in dlc_param.DOLPHOT_PARAMS:
        f.write(f'{key} = {str(dlc_param.DOLPHOT_PARAMS[key])}\n')
    for key in extra_params:
        f.write(f'{key} = {str(extra_params[key])}\n')
    f.close()


def prep_files_for_dolphot(image_directory, r_in, r_out,
                           step, sig_low, sig_high, dlc_param):
    '''
    Processes images in a directory through masking, spliting, and calcsky
    programs in Dolphot

    Parameters
    ----------
    image_directory : str
        Path to directory with images

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

    prepped_dir = f'{dlc_param.IMROOT}{image_directory}'

    # FLAG
    if image_directory == '/dolphot_prepped':
        try:
            os.mkdir(prepped_dir)
        except FileExistsError:
            shutil.rmtree(prepped_dir)
            os.mkdir(prepped_dir)

        for im in dlc_param.IMAGES:
            shutil.copyfile(im.loc, im.prep_loc)

    os.chdir(prepped_dir)

    for im in dlc_param.IMAGES:
        os.system(f'{dlc_param.DOLPHOT_PATH}{dlc_param.MASK} '
                  f'{im.name}_{dlc_param.SUFFIX}.fits')
        os.system(f'{dlc_param.DOLPHOT_PATH}/splitgroups '
                  f'{im.name}_{dlc_param.SUFFIX}.fits')
        for chip in dlc_param.CHIPS:
            os.system(f'{dlc_param.DOLPHOT_PATH}/calcsky '
                      f'{im.name}_{dlc_param.SUFFIX}.chip{chip} {r_in} '
                      f'{r_out} {step} {sig_low} {sig_high}')

    os.chdir(dlc_param.REF_IMAGE_PATH)

    os.system(f'{dlc_param.DOLPHOT_PATH}/splitgroups {dlc_param.REF_IMAGE}')

    os.chdir(dlc_param.IMROOT)

    # Refresh the files as AstroDrizzle changed the header info
    files = os.listdir(dlc_param.ORIG_IM_LOC)
    for file in files:
        shutil.copyfile(f'{dlc_param.ORIG_IM_LOC}/{file}',
                        f'{dlc_param.IMROOT}/Images/{file}')


def dolphot_force(objCoords, dlc_param, apermag=False,
                  force_same_mag=True, psfphot=1):
    '''
    Runs Dolphot on difference images

    Parameters
    ----------

    objCoords : dict
        Dictionary of transient object coordinates
        (Ex: {'S1': [177.398231, 22.395630], ...})

    dlc_param : obj
        Parameter object from prep_directory function

    apermag : bool
        Run aperture photometry

    force_same_mag : bool
        Run photometry assuming same count rate

    psfphot : float
        Type of photometry to run
        0 = aperture
        1 = standard PSF-fit
        2 = PSF-fit weighted for central pixels

    '''

    imdir_prepped = f'{dlc_param.IMROOT}/dolphot_prepped'
    imdir_simultaneous = f'{dlc_param.IMROOT}/diffs'

    files = os.listdir(imdir_prepped)
    N = len(files)
    i = N - 1
    while i > 0:
        if 'output' not in files[i]:
            del files[i]
        i = i - 1

    files.append('registration.chip1.fits')

    for fname in files:
        shutil.copy(f'{imdir_prepped}/{fname}',
                    f'{imdir_simultaneous}/{fname}')

    imgNum = 0
    extra_params = {}
    info_params = {}

    ref_image_use = dlc_param.REF_IMAGE

    files = [a.loc for a in dlc_param.IMAGES]

    imdir = f'{dlc_param.IMROOT}/dolphot/'
    imdir_dolphot_prepped = f'{dlc_param.IMROOT}/diffs/'

    for file in files:
        for chip in dlc_param.CHIPS:

            if imgNum >= 99:
                raise Exception

            fname_prepped = f'{imdir_dolphot_prepped}'\
                            f'{file.split("/")[-1].split("_")[0]}_'\
                            f'{dlc_param.SUFFIX}'\
                            f'.chip{str(chip)}'

            if chip == dlc_param.ON_CHIP:
                imgNum += 1

                fname_simultaneous = f'{imdir_simultaneous}{file}_'\
                                     f'{dlc_param.SUFFIX}.chip{chip}'

                cmd = f'cp {fname_prepped}.sky.fits '\
                      f'{fname_simultaneous}.sky.fits'

                extra_params[f'img{imgNum}_file'] = \
                    f'{file.split("/")[-1].split("_")[0]}_'\
                    f'{dlc_param.SUFFIX}.chip{str(chip)}'
                extra_params[f'img{imgNum}_shift'] = '0 0'
                extra_params[f'img{imgNum}_xform'] = '1 0 0'

                if dlc_param.INST != 'WFPC2':
                    info_params[f'img{imgNum}_instrument'] = dlc_param.INST
                    info_params[f'img{imgNum}_detector'] = dlc_param.DETEC
                    info_params[f'img{imgNum}_filt'] = dlc_param.FILT

                    orig = f'{dlc_param.IMROOT}/imaging/{file}_'\
                           f'{dlc_param.SUFFIX}.fits'
                    masked = f'{dlc_param.IMROOT}/dolphot/{file}_'\
                             f'{dlc_param.SUFFIX}.fits'

                    info_params[f'img{imgNum}_orig'] = orig
                    info_params[f'img{imgNum}_masked'] = masked

                if dlc_param.INST == 'WFPC2':
                    command = f'gethead {fname_simultaneous}.fits EXPNAME'

                    namef = subprocess.getoutput(command)
                    orig_crclean = f'{dlc_param.IMROOT}/imaging/'\
                                   f'{namef}_C0M_crclean.fits'
                    orig = f'{dlc_param.IMROOT}/imaging/{namef}_C0M.fits'
                    orig_dq = f'{dlc_param.IMROOT}/imaging/{namef}_C1M.fits'
                    fn_dp_masked = f'{dlc_param.IMROOT}/dolphot/'\
                                   f'{namef}_C0M.fits'

                    info_params[f'img{imgNum}_orig_crclean'] = orig_crclean
                    info_params[f'img{imgNum}_orig'] = orig
                    info_params[f'img{imgNum}_orig_dq'] = orig_dq
                    info_params[f'img{imgNum}_dolphot_masked'] = fn_dp_masked

                    info_params[f'img{imgNum}_instrument'] = dlc_param.INST
                    info_params[f'img{imgNum}_detector'] = dlc_param.DETEC
                    info_params[f'img{imgNum}_filt'] = dlc_param.FILT

    extra_params['Nimg'] = imgNum

    if True:
        extra_params['img0_file'] = ref_image_use.replace('.fits', '.chip1')
        extra_params['img0_RAper'] = '4'
        extra_params['img0_RChi'] = '2.0'
        extra_params['img0_RSky'] = '15 35'
        extra_params['img0_RPSF'] = '15'

    param_file = f'{imdir_simultaneous}/dolphot.params'

    os.chdir(imdir_simultaneous)

    '''
    here using the recommended settings for a WFC3 UVIS registration image
    '''
    _mk_param(dlc_param.INST, dlc_param.DETEC, param_file, extra_params,
              dlc_param.IMTYPE, dlc_param)
    # End

    xytfile = open('xytfile', 'w')

    objs = []

    for key in objCoords.keys():

        small_ra, small_dec = objCoords[key]

        from astropy.wcs import WCS

        from astropy.io import fits

        try:
            w = WCS(fits.open(f'{dlc_param.REF_IMAGE_PATH}/'
                              f'{dlc_param.REF_IMAGE}')['SCI'])
        except KeyError:
            w = WCS(fits.open(f'{dlc_param.REF_IMAGE_PATH}/'
                              f'{dlc_param.REF_IMAGE}'))

        import astropy.units as u

        import astropy.coordinates as coord
        ra = coord.Angle(small_ra, unit=u.hour)  # pylint: disable = no-member
        ra_deg = ra.degree

        dec = coord.Angle(small_dec, unit=u.degree)\
            # pylint: disable = no-member
        dec_deg = dec.degree

        ''' need to translate '''
        big_x, big_y = w.wcs_world2pix(small_ra, small_dec, 1,
                                       ra_dec_order=True)

        objs.append([key, big_x, big_y])

        xytfile.write(f'0 1 {big_x} {big_y} 2 10\n')
    xytfile.close()

    if apermag:
        cmd = f'{dlc_param.DOLPHOT_PATH}/dolphot singlestar -pdolphot.params '\
              f'xytfile=xytfile usephot=output PSFPhot=0 Force1=1 SigFind=-99'\
              f' Force1=1 SigFindMult=1.0 SigFinal=-99'
    else:
        cmd = f'{dlc_param.DOLPHOT_PATH}/dolphot singlestar -pdolphot.params '\
              f'xytfile=xytfile usephot=output PSFPhot={psfphot} Force1=1 '\
              f'FitSky=1 SigFind=-99 SigFindMult=1.0 SigFinal=-99'

        if force_same_mag:
            cmd += ' ForceSameMag=1'
        else:
            cmd += ' ForceSameMag=0'

    os.system(cmd)

    statinfo = os.stat('singlestar')
