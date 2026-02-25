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

        output = subprocess.run( [f'gethead',f'{file}', 'FILTER'], capture_output=True, text=True )
        filt = output.stdout[:-1]

        output = subprocess.run( [f'gethead',f'{file}', 'PROPOSID'], capture_output=True, text=True )
        proposid = output.stdout[:-1]



        if True: # filt == 'F350LP' and proposid == '15936':

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

        cmd = f'{dlc_param.DOLPHOT_PATH}wfc3mask '\
                  f'{ref_image_use.split("/")[-1]}'
        #os.system('pwd')
        print(cmd)
        os.system(cmd)



        refim_c1 = ref_image_use.replace(".fits", ".chip1.fits").split("/")[-1]
        if not glob(f'{imdir_simultaneous}{refim_c1}'):
            pass

    if True:
        extra_params['img0_file'] = ref_image_use.replace('.fits', '.chip1')
        extra_params['img0_RAper'] = '4'
        extra_params['img0_RChi'] = '2.0'
        extra_params['img0_RSky'] = '15 35'
        extra_params['img0_RPSF'] = '15'

    shutil.copyfile(f'{dlc_param.REF_IMAGE_PATH}/'      # /example/Aligned_Images
                    f'{ref_image_use.replace(".fits",".chip1.fits")}',
                    f'{imdir_simultaneous}/'            # /example/dolphot_prepped
                    f'{ref_image_use.replace(".fits",".chip1.fits")}') #/example/dolphot_prepped

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


def run_mask_split(inarr):    

    dolphot_path, mask, im, suffix = inarr

    cmd = f'{dolphot_path}{mask} {im}_{suffix}.fits'
              

    print(cmd)
    os.system(cmd)


    cmd = f'{dolphot_path}/splitgroups {im}_{suffix}.fits'

    print(cmd)
    os.system(cmd)



def run_calcsky(inarr):    

    dolphot_path, im, suffix, chip, r_in, r_out, step, sig_low, sig_high = inarr

    cmd = f'{dolphot_path}/calcsky {im}_{suffix}.chip{chip} {r_in} {r_out} {step} {sig_low} {sig_high}'

    print(cmd)
    os.system(cmd)





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

    from concurrent.futures import ProcessPoolExecutor

    
    
    

    
    

    

    

    if True:

        from copy import deepcopy                                                                                                                         
        inarr = [ [deepcopy(dlc_param.DOLPHOT_PATH), deepcopy(dlc_param.MASK), deepcopy(im.name), deepcopy(dlc_param.SUFFIX)] for im in dlc_param.IMAGES]
                                                                                                                                                          

                                                                                                                                                          
                                                                                                                                                          
        with ProcessPoolExecutor() as executor:
            for output in executor.map(run_mask_split, inarr):
                print(output)
                                                                                                                                                          
        inarr = []
        for im in dlc_param.IMAGES:
            for chip in dlc_param.CHIPS:
                inarr.append([dlc_param.DOLPHOT_PATH, im.name, dlc_param.SUFFIX, chip, r_in, r_out, step, sig_low, sig_high])
                                                                                                                                                          
        with ProcessPoolExecutor() as executor:
            for output in executor.map(run_calcsky, inarr):
                print(output)

    print(dlc_param)

    print(os.getcwd())    

    os.chdir(dlc_param.IMROOT)

    os.chdir(dlc_param.REF_IMAGE_PATH)

    print(dlc_param.REF_IMAGE_PATH)

    ''' set GAIN to EXPTIME in header since coadd in units of e- / s '''
    import subprocess
    output = subprocess.run( [f'gethead',f'{dlc_param.REF_IMAGE}', 'EXPTIME'], capture_output=True, text=True )

    exptime = output.stdout[:-1]
                                                                                                           
    print('exptime', exptime)
    print(f'dlc_param.REF_IMAGE')

    cmd = f'sethead {dlc_param.REF_IMAGE} GAIN={exptime}'
    print(cmd)                                                                                                           
    os.system( cmd )
                                                                                                           


    os.system(f'{dlc_param.DOLPHOT_PATH}/splitgroups {dlc_param.REF_IMAGE}')

    os.chdir(dlc_param.IMROOT)

    # Refresh the files as AstroDrizzle changed the header info
    files = os.listdir(dlc_param.ORIG_IM_LOC)
    for file in files:
        shutil.copyfile(f'{dlc_param.ORIG_IM_LOC}/{file}',
                        f'{dlc_param.IMROOT}/Images/{file}')


def dolphot_force(objCoords, dlc_param, apermag=False,
                  force_same_mag=True, psfphot=1, runFake=False):
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

    os.system(f'cp {imdir_prepped}/output* {imdir_simultaneous}')

    #imdir_simultaneous = f'{dlc_param.IMROOT}/dolphot_prepped'

    files = os.listdir(imdir_prepped)
    N = len(files)
    i = N - 1
    while i > 0:
        if 'output' not in files[i]:
            del files[i]
        i = i - 1

    ### files.append('registration.chip1.fits')
    files.append(dlc_param.REF_IMAGE.replace('.fits','.chip1.fits'))


    #for fname in files:
    #    if imdir_prepped != imdir_simultaneous:
    #        shutil.copy(f'{imdir_prepped}/{fname}',
    #                f'{imdir_simultaneous}/{fname}') # copies 1st file to replace 2nd


    imgNum = 0
    extra_params = {}
    info_params = {}

    ref_image_use = dlc_param.REF_IMAGE

    files = [a.loc for a in dlc_param.IMAGES]

    imdir = f'{dlc_param.IMROOT}/dolphot/'
    imdir_dolphot_prepped = f'{dlc_param.IMROOT}/dolphot_prepped/'

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

        print('cwd', os.getcwd() )

        try:
            w = WCS(fits.open(f'../{dlc_param.REF_IMAGE_PATH}/'
                              f'{dlc_param.REF_IMAGE}')['SCI']) ###
            print("w is chosen here 1")
        except KeyError:
            try:
                w = WCS(fits.open(f'../{dlc_param.REF_IMAGE_PATH}/'   ###
                              f'{dlc_param.REF_IMAGE}')[0])
                print("w is chosen here 2")

            except:
                w = WCS(fits.open(f'../{dlc_param.REF_IMAGE_PATH}/'   ###
                              f'{dlc_param.REF_IMAGE}'))
                print("w is chosen here 3")


        ra = coord.Angle(small_ra, unit=u.hour)  # pylint: disable = no-member
        ra_deg = ra.degree

        dec = coord.Angle(small_dec, unit=u.degree)\
            # pylint: disable = no-member
        dec_deg = dec.degree

        ''' need to translate '''
        #big_x, big_y = 2170.59, 2713.410 # w.wcs_world2pix(small_ra, small_dec, 1,
                                       #ra_dec_order=True)


        big_x, big_y = w.wcs_world2pix(small_ra, small_dec, 1,
                                       ra_dec_order=True)

        print("HERE IN D_FORCE big_x, big_y = ", big_x, big_y)

        objs.append([key, big_x, big_y])

        xytfile.write(f'0 1 {big_x} {big_y} 2 10\n')


        side = 4
        photsec = '0 1 %d %d %d %d' % (big_x - side, big_y - side, big_x + side, big_y + side)



    
    xytfile.write(f'0 1 2164.6681 2709.3879 2 10\n')
    xytfile.write(f'0 1 2166.2951 2714.7855 2 10\n')
    xytfile.close()

    if apermag:
        cmd = f'{dlc_param.DOLPHOT_PATH}/dolphot singlestar -pdolphot.params '\
              f'xytfile=xytfile usephot=output PSFPhot=0 Force1=1 SigFind=-99'\
              f' Force1=1 SigFindMult=1.0 SigFinal=-99'
    else:
        cmd = f'{dlc_param.DOLPHOT_PATH}/dolphot singlestar -pdolphot.params '\
              f'xytfile=xytfile usephot=output PSFPhot={psfphot} Force1=1 '\
              f'FitSky=2 SigFind=-99 SigFindMult=1.0 SigFinal=-99'

        #cmd = f'{dlc_param.DOLPHOT_PATH}/dolphot singlestar -pdolphot.params '\
        #f'usephot=output PSFPhot={psfphot} '\
        #      f'FitSky=2 SigFind=5 SigFindMult=1.0 SigFinal=-99 photsec="' + photsec + '"'




        if force_same_mag:
            cmd += ' ForceSameMag=1'
        else:
            cmd += ' ForceSameMag=0'

    print('Running DOLPHOT forced photometry')
    print(cmd)
    os.system(cmd)
    print('Finished running DOLPHOT forced photometry')

    #os.system("pwd") # prints current working directory, ensure 'singlestar' present

    statinfo = os.stat('singlestar') # path



    #points = get_random_points_in_regions('my_regions.reg', 10, shape=(2000, 2000))



def _sample_points_in_region(region_file, n, wcs):
    """
    Sample n random pixel positions from within the shapes defined in a DS9
    region file.

    Parameters
    ----------
    region_file : str
        Path to DS9 region file (.reg)
    n : int
        Number of random points to sample
    wcs : astropy.wcs.WCS
        WCS of the reference image, used to convert sky regions to pixel coords

    Returns
    -------
    list of (x, y) tuples
        Pixel coordinates (1-indexed, FITS convention)
    """
    from regions import Regions, PixCoord
    import random as _random

    regs = Regions.read(region_file, format='ds9')
    if not regs:
        raise ValueError(f"No regions found in {region_file}")

    # Convert any sky regions to pixel regions
    pixel_regs = []
    for reg in regs:
        if hasattr(reg, 'to_pixel'):
            pixel_regs.append(reg.to_pixel(wcs))
        else:
            pixel_regs.append(reg)

    points = []
    max_attempts = n * 500

    for _ in range(max_attempts):
        if len(points) >= n:
            break
        reg = _random.choice(pixel_regs)
        bb = reg.bounding_box
        x = _random.uniform(bb.ixmin, bb.ixmax)
        y = _random.uniform(bb.iymin, bb.iymax)
        if reg.contains(PixCoord(x, y)):
            # regions uses 0-indexed pixels; DOLPHOT/FITS expects 1-indexed
            points.append((x + 1.0, y + 1.0))

    if len(points) < n:
        raise ValueError(
            f"Could only sample {len(points)}/{n} points from region file "
            f"'{region_file}' after {max_attempts} attempts. "
            f"Region may be too small."
        )
    return points


def fake_phot(dlc_param, objCoords, psfphot, filt, key, nstar=50,
              region_file=None):

    from astropy.io import fits


    imdir_simultaneous = f'{dlc_param.IMROOT}/diffs'

    #imdir_simultaneous = f'{dlc_param.IMROOT}/dolphot_prepped'
    os.chdir(imdir_simultaneous)

    ''' generate fakelist '''

    import subprocess
    #cmd = f'{dlc_param.DOLPHOT_PATH}/fakelist output {filt} WFC3_F350LP 28.0 28.1 0 0.1 -NSTAR={nstar} &> Fake.list'
    #print(cmd)

    #output = subprocess.run([f'{dlc_param.DOLPHOT_PATH}/fakelist', 'output', f'{filt}', 'WFC3_F350LP', '28.0', '28.1', '0', '0.1', f'-NSTAR={nstar}'], capture_output=True)
    #print(output)


    #f = open('Fake.list','wb')
    #f.write(output.stdout)
    #f.close()

    ''' create Fake.list with a bunch of fakes '''
    #x, y = 2176, 2711 

    small_ra, small_dec = objCoords[key]
                                                                                
    print('cwd', os.getcwd() )
                                                                                
    try:
        w = WCS(fits.open(f'../{dlc_param.REF_IMAGE_PATH}/'
                          f'{dlc_param.REF_IMAGE}')['SCI']) ###
        print("w is chosen here 1")
    except KeyError:
        try:
            w = WCS(fits.open(f'../{dlc_param.REF_IMAGE_PATH}/'   ###
                          f'{dlc_param.REF_IMAGE}')[0])
            print("w is chosen here 2")
                                                                                
        except:
            w = WCS(fits.open(f'../{dlc_param.REF_IMAGE_PATH}/'   ###
                          f'{dlc_param.REF_IMAGE}'))
            print("w is chosen here 3")
                                                                                
                                                                                
    ra_star = coord.Angle(small_ra, unit=u.hour)  # pylint: disable = no-member
    ra_deg = ra_star.degree
                                                                                
    dec_star = coord.Angle(small_dec, unit=u.degree)\
        # pylint: disable = no-member
    dec_deg = dec_star.degree
                                                                                
    ''' need to translate '''
    #big_x, big_y = 2170.59, 2713.410 # w.wcs_world2pix(small_ra, small_dec, 1,
                                   #ra_dec_order=True)
                                                                                
                                                                                
    x_transient, y_transient = w.wcs_world2pix(small_ra, small_dec, 1,
                                   ra_dec_order=True)
                                                                                
    print("HERE IN D_FORCE big_x, big_y = ", x_transient, y_transient)



    mag = 28.5

    import random

    f = open('Fake.list','w')
    f2 = open('Fake.xy','w')

    xps = []
    yps = []

    if region_file is not None:
        print(f"Sampling {nstar} fake-star positions from region file: {region_file}")
        sampled = _sample_points_in_region(region_file, nstar, w)
        for xp, yp in sampled:
            ra, dec = w.wcs_pix2world(xp, yp, 1)
            xps.append(ra)
            yps.append(dec)
            f.write(f'0 1 {xp:.2f} {yp:.2f} {mag:.2f}\n')
            f2.write(f'0 1 {xp:.2f} {yp:.2f} 1 10\n')
    else:
        for i in range(nstar):
            xp_shift = random.choice(range(-200,200))
            yp_shift = random.choice(range(-200,200))

            xp = x_transient + xp_shift
            yp = y_transient + yp_shift

            ra, dec = w.wcs_pix2world( xp, yp, 1 )

            xps.append(ra)
            yps.append(dec)

            f.write(f'0 1 {xp:.2f} {yp:.2f} {mag:.2f}\n')
            f2.write(f'0 1 {xp:.2f} {yp:.2f} 1 10\n')

    f.close()
    f2.close()


    ''' addstars '''
    cmd = f'{dlc_param.DOLPHOT_PATH}/addstars output -pdolphot.params '\
          f'usephot=output FakeStars=Fake.list '  # check FakeStarPSF
    print(cmd)
    os.system(cmd)

    ''' create dolphotFake.params and reset filenames '''
    f = open('dolphot.params','r').readlines()
    o = open('dolphotFake.params','w')
                                                                                                 
    for l in f:
        import re
        res = re.split('\s+', l)
        if res[0][:3] == 'img' and l.find('_file') != -1 and l.find('img0') == -1:

            imInd = int(float(res[0].split('_')[0][3:])) 

            ifile = res[2] + '.fits' 
            ofile = 'output.%d.mod.fits' % imInd

            print('cwd', os.getcwd() )

            print(ifile)



            if ifile.find('chip1') != -1:
                w = WCS(fits.open('../Images/' + ifile.replace('.chip1',''))['SCI',1].header) ###
            else:
                print(ifile.replace('.chip2',''))


                fl = '../Images/' + ifile.replace('.chip2','')

                print(fl)

                ''' use DOLPHOT routines instead !!! '''
                w = WCS(fits.open(fl)[2].header) ###

            ''' need to translate '''
            #big_x, big_y = 2170.59, 2713.410 # w.wcs_world2pix(small_ra, small_dec, 1,
                                           #ra_dec_order=True)
                                                                                        
                                                                                        
            
            
                                                                                        
            print("HERE IN D_FORCE big_x, big_y = ", x_transient, y_transient)


            from astropy.io import fits

            ''' now copy over negative pixels, since added sources can demask by making negative pixels positive '''
            ih = fits.open(ifile)
            oh = fits.open(ofile)

            oh[0].data[ ih[0].data < 0 ] = ih[0].data[ ih[0].data < 0 ]

            ''' also copy over negative pixels '''


            
            x_transient, y_transient= w.wcs_world2pix(small_ra, small_dec, 1,
                               ra_dec_order=True)

            print(x_transient, y_transient) 




            width = 10 
            m_neg = oh[0].data[ int(y_transient) - width : int(y_transient) + width , int(x_transient) - width : int(x_transient) + width  ] < 0

            print('m_neg', m_neg)


            for ra, dec in zip(xps, yps):


                x, y= w.wcs_world2pix(ra, dec, 1,
                               ra_dec_order=True)


                print(x, y)




                #oh[0].data[ int(y_transient) + ysh - 4 : int(y_transient) + ysh + 4, int(x_transient) + xsh - 4 : int(x_transient) + xsh + 4][m_neg] = 99999

                import numpy

                #print('here',99999. * numpy.array( m_neg, dtype=int))

                print(m_neg.shape)

                oh[0].data[ int(y) - width : int(y) + width, int(x) - width : int(x) + width] += -99999. * numpy.array( m_neg, dtype=int)


                print(x_transient, y_transient)

            oh.writeto(ofile, overwrite=True)

            os.system('pwd')
            print(ofile)

            o.write(res[0] + ' = output.%d.mod\n' % imInd)
        else:
            o.write(l)
    o.close()

    ''' run forced photometry on .mod files '''
    cmd = f'{dlc_param.DOLPHOT_PATH}/dolphot fakes -pdolphotFake.params '\
              f'xytfile=Fake.xy usephot=output FakeStarPSF=1 Force1=1 SigFind=-99 SigFindMult=0.85 SigFinal=-99' # PSFPhot={psfphot} ' #Force1=1 '\
#              f'FitSky=2 SigFind=-99 SigFindMult=0.5 SigFinal=-99'
    print(cmd)
    os.system(cmd)

    os.system('pwd')


    ''' feed into loadtxt and compute errors '''
    import numpy



    lf = numpy.loadtxt('fakes')                                                                                                                                  
    print(lf[:,27])
    print(lf[:,24])
    print(lf[:,32])

    std1 = (numpy.percentile(lf[:,27], 84) - numpy.percentile(lf[:,27], 16) ) / 2.
    print(std1)

    print('median flux', numpy.median(lf[:,27]))
    print(numpy.median(lf[:,28]))
    print('----')
    print(lf[:,40])
    print(lf[:,37])
    print(lf[:,45])
    std2 = (numpy.percentile(lf[:,40], 84) - numpy.percentile(lf[:,40], 16) ) / 2.
    print(std2)

    print('median flux', numpy.median(lf[:,40]))
    print(numpy.median(lf[:,41]))

    print('----')
    std3 = (numpy.percentile(lf[:,53], 84) - numpy.percentile(lf[:,53], 16) ) / 2.
    print(std3)

    print('median flux', numpy.median(lf[:,53]))
    print(numpy.median(lf[:,54]))

    print('----')
    std4 = (numpy.percentile(lf[:,66], 84) - numpy.percentile(lf[:,66], 16) ) / 2.
    print(std4)

    print('median flux', numpy.median(lf[:,66]))
    print(numpy.median(lf[:,67]))


    a = open('stats.txt','w')
    a.write(f'{std1}\n{std2}\n{std3}\n{std4}\n')
    a.close()

    for i in range(1,5):
        os.system(f'gethead output.{i}.mod.fits EXPTIME')

    #print(numpy.std(lf[:,40]))
    #print(numpy.std(lf[:,53]))
    #print(numpy.std(lf[:,66]))






