import os
import shutil
from drizzlepac import astrodrizzle
from drizzlepac import tweakreg
from drizzlepac import tweakback
import scipy


def coadd_ims(dlc_param):
    '''
    Creates coadded template image from raw template images

    Parameters
    ----------
    dlc_param : obj
        Parameter object from prep_directory function
    '''

    image_list = os.listdir('image_backup')

    shutil.copyfile('default.conv', './Images/default.conv')
    shutil.copyfile('default.param', './Images/default.param')

    try:
        os.mkdir('coadd')
    except FileExistsError:
        pass

    flist = os.listdir(f'{dlc_param.IM_LOC}/template')
    files_found_backup = [f'{dlc_param.IM_LOC}/template/{x}'
                          for x in flist]

    i = len(files_found_backup) - 1
    while i > -1:
        if '.fits' not in files_found_backup[i]:
            del files_found_backup[i]
        i = i - 1

    refim = f'{dlc_param.REF_IMAGE_PATH}/{dlc_param.REF_IMAGE}'

    crCleanFirst = 1
    crCleanFirstUseSaved = 0
    runtweakreg = 1
    runtweakregUseSaved = 0
    coaddastrom = 1
    coaddastromUseSaved = 0
    adjustwithcoadd = 1
    recoaddwithadjust = 1
    meassky = 1
    recoaddskyfix = 1
    recoaddglobalmin = 0
    recoaddnorthup = 0

    files_backup = []

    '''
    threshold = 10
    refim_align ='./ref/registration.fits'
    _make_template_cat(refim_align, 1, filt, threshold)
    sys.exit()
    #'''

    for fname in files_found_backup:
        from astropy.io import fits
        p = fits.open(fname)
        try:
            filt_exp = p[0].header['FILTER']
        except KeyError:
            filt_exp = p[0].header['FILTER2']

        if dlc_param.FILT == filt_exp:
            files_backup.append(fname)

    files = [x.replace('_backup', '') for x in files_backup][:]

    ims = files

    if crCleanFirst:
        if not crCleanFirstUseSaved:

            if len(files) <= 6:
                combine_type = 'minmed'
            else:
                combine_type = 'imedian'

            os.chdir(dlc_param.IM_LOC)

            ''' localmin since we just want this for
                source detection and astrometry'''
            for cfiles in [ims]:
                astrodrizzle.AstroDrizzle(cfiles, driz_cr_corr=True,
                                          output=str('output'),
                                          num_cores=8,
                                          combine_type=combine_type,
                                          final_wcs=True, skysub=True,
                                          skymethod='localmin',
                                          preserve=False, )
            try:
                os.mkdir('TEMP_nocr')
            except FileExistsError:
                pass

            cfiles_nocr = [x.replace('_flc.fits', '_crclean.fits') for
                           x in (ims)]

            for fname in cfiles_nocr:
                shutil.copy(fname, f'./TEMP_nocr/')

        else:
            shutil.copytree('./TEMP_nocr/', './IMS/')
            shutil.copytree('./Images/', './IMS/')

    fitgeometry = 'rscale'
    nclip = 3
    minobj = 10

    if runtweakreg:
        try:
            os.mkdir('TEMP_aligned')
        except FileExistsError:
            pass

        if not runtweakregUseSaved:
            for cfiles in [ims]:

                cfiles_nocr = [x.replace('_flc.fits', '_crclean.fits') for
                               x in cfiles]

                catName = 'astdriz_catfile.list'

                f = open(catName, 'w')

                fnames = []

                for im in cfiles_nocr:
                    ''' manually compute sigma '''

                    im = os.path.abspath(im)

                    im_short = im.replace('//', '/')

                    threshold = 10

                    _make_sextractor_cat1(im, 1, dlc_param.FILT, threshold)
                    _make_sextractor_cat1(im, 4, dlc_param.FILT, threshold)

                    a = im_short.replace('.fits', '')
                    b = im_short.replace('.fits', '')
                    f.write(f'{im_short} {a}_ref_1.cat {b}_ref_4.cat\n')

                    fnames.append(im_short)

                f.close()

                tweakreg.TweakReg(fnames[:], catfile=catName, xcol=1,
                                  ycol=2, updatehdr=True,  nclip=5,
                                  peakmax=50000, sigma=2.5, searchrad=10.0,
                                  tolerance=5.0, writecat=True,
                                  headerlet=True, attach=False,
                                  clobber=True, minobj=-1,
                                  fitgeometry='general',  wcsname="TWEAK1",
                                  shiftfile=True,
                                  outshifts='shift_file.txt',
                                  residplot='No plot', see2dplot=False,
                                  interactive=False)

                from astropy.table import Table
                shift_tab = Table.read('shift_file.txt',
                                       format='ascii.no_header',
                                       names=['file', 'dx', 'dy', 'rot',
                                              'scale', 'xrms', 'yrms'])

                formats = ['.2f', '.2f', '.3f', '.5f', '.2f', '.2f']
                for i, col in enumerate(shift_tab.colnames[1:]):
                    shift_tab[col].format = formats[i]

                cfiles_nocr = [x.replace('_flc.fits', '_crclean.fits') for
                               x in cfiles]

                from stwcs.wcsutil import headerlet
                for fname, fname_orig in zip(cfiles_nocr, cfiles):

                    '''  some of these settings are necessary so that
                            DOLPHOT doesn't choke later'''

                    headerlet.apply_headerlet_as_primary(
                        fname_orig, fname.replace('.fits', '_hlet.fits'),
                        attach=False, archive=False)

            for fname in files:
                shutil.copy(fname, f'./TEMP_aligned/')

        else:
            shutil.copytree('./TEMP_aligned/', './IMS/')

    if coaddastrom:

        from stwcs.wcsutil import headerlet

        if not coaddastromUseSaved:

            for file in ims:
                default_wcsname = fits.getval(file, 'wcsname', ext=1)

            for cfiles in [[ims]]:
                if len(cfiles) <= 6:
                    combine_type = 'minmed'
                else:
                    combine_type = 'imedian'

                ''' both together '''

                ofile = f'coadd_{dlc_param.FILT}.fits'

                astrodrizzle.AstroDrizzle(
                    [f for f in cfiles], output=ofile, final_wcs=True,
                    driz_cr_corr=True, num_cores=48,
                    combine_type=combine_type, skysub=True,
                    skymethod='localmin', skystat='mode', build=True,
                    wcskey='TWEAK1', preserve=False)

                os.mkdir('./coadded_ims')
                shutil.copy(ofile, f'./coadded_ims/')

        else:
            shutil.copytree('./coadded_ims/', './')

    threshold = 10

    im_tweak = f'coadd_{dlc_param.FILT}.fits'
    refim_align = f'{dlc_param.IMROOT}/ref/registration.fits'

    if adjustwithcoadd:

        catName = 'astdriz_catfile.list'
        refcat = 'astdriz_catfile_ref.list'

        f = open(catName, 'w')
        _make_sextractor_cat1(im_tweak, 1, dlc_param.FILT,
                              threshold, maxobjs=3000)
        a = im_tweak.replace('.fits', '')
        f.write(f'{im_tweak} {a}_ref_1.cat\n')
        f.close()

        _make_template_cat(refim_align, 1, dlc_param.FILT,
                           threshold, maxobjs=3000)
        a = refim_align.replace('.fits', '')
        refcat = f'{a}_ref_1.cat'

        a = im_tweak.replace('.fits', '_orig.fits')
        shutil.copy(im_tweak, a)

        tweakreg.TweakReg([im_tweak], catfile=catName, xcol=1, ycol=2,
                          interactive=False, refimage=refim_align,
                          refcat=refcat, refxcol=1, refycol=2,
                          refxyunits='pixels', updatehdr=True, nclip=3,
                          peakmax=50000, sigma=2.5, searchrad=1,
                          wcsname='TWEAKCOADD', writecat=False,
                          headerlet=False, attach=False,  clobber=True,
                          fitgeometry='general', tolerance=5.0,
                          shiftfile=True, outshifts='shift_file_coadd.txt',
                          residplot='No plot', see2dplot=False, minobj=-1)

        for file in [im_tweak]:
            from stwcs.wcsutil import headerlet

            default_wcsname = fits.getval(file, 'wcsname', ext=1)

        ''' reversing the WCS names '''
        tweakback.tweakback(im_tweak, input=ims, verbose=True, force=True,
                            wcsname='TWEAK1', newname='TWEAK2')
    coadd_nohff_tweak = f'coadd_nohff_tweak_{dlc_param.FILT}.fits'

    if recoaddwithadjust:

        for cfiles, year in [[files, '2021']]:

            if len(cfiles) <= 6:
                combine_type = 'minmed'
            else:
                combine_type = 'imedian'

            for file in cfiles:

                from stwcs.wcsutil import headerlet

                default_wcsname = fits.getval(file, 'wcsname', ext=1)

            ''' both together '''
            astrodrizzle.AstroDrizzle(
                [f for f in cfiles[:]], output=coadd_nohff_tweak,
                final_wcs=True, final_refimage=refim,  driz_cr_corr=True,
                num_cores=48, combine_type=combine_type, skysub=True,
                skymethod='localmin', skystat='mode', build=True,
                preserve=False)

    if meassky:

        from astropy.io import fits

        from scipy import stats

        for cfiles in [ims]:

            skyfile = open('skyfile.txt', 'w')

            skyfile_mode = open('skyfile_mode.txt', 'w')

            cfiles_nocr = [x.replace('_flc.fits', '_crclean.fits') for
                           x in cfiles]

            for im_nocr, im in zip(cfiles_nocr, cfiles):

                im_nocr = os.path.abspath(im_nocr)

                threshold = 5
                fitscat, fname = _make_sextractor_cat1(
                                    im_nocr, 1, dlc_param.FILT, threshold,
                                    bgChipEstimate=True)

                bg_1 = fits.open(fitscat)[2].data['BACKGROUND'][0]

                f = fits.open(im, mode='update')
                f[1].header['PKSKY'] = bg_1

                mask = scipy.array(f[3].data).flatten() == 0
                mode_1 = scipy.stats.mode(
                            scipy.array(f[1].data).flatten()[mask])[0][0]

                f[1].header['PKSKYMODE'] = mode_1

                fitscat, fname = _make_sextractor_cat1(
                                    im_nocr, 4, dlc_param.FILT, threshold,
                                    bgChipEstimate=True)

                bg_4 = fits.open(fitscat)[2].data['BACKGROUND'][0]

                f[4].header['PKSKY'] = bg_4

                mask = scipy.array(f[6].data).flatten() == 0
                mode_4 = scipy.stats.mode(
                            scipy.array(f[4].data).flatten()[mask])[0][0]

                f[4].header['PKSKYMODE'] = mode_4

                f.flush()

                skyfile.write(f'{os.path.abspath(im)} {bg_1} {bg_4}\n')

                skyfile_mode.write(
                    f'{os.path.abspath(im)} {mode_1} {mode_4}\n')

            skyfile.close()

            skyfile_mode.close()

    if recoaddskyfix:

        for cfiles in [ims]:

            if len(cfiles) <= 6:
                combine_type = 'minmed'
            else:
                combine_type = 'imedian'

            astrodrizzle.AstroDrizzle(
                [f for f in cfiles[:]],
                output=f'coadd_tweak_{dlc_param.FILT}.fits',
                final_wcs=True, driz_cr_corr=True, num_cores=64,
                combine_type=combine_type, skysub=True, build=False,
                skyfile='skyfile_mode.txt', skyuser='', preserve=False,
                final_refimage=refim)

    if recoaddglobalmin:

        for cfiles, year in [[ims, '2020']]:

            if len(cfiles) <= 6:
                combine_type = 'minmed'
            else:
                combine_type = 'imedian'

            astrodrizzle.AstroDrizzle(
                [f for f in cfiles[:]],
                output=f'coadd_glotweak_{dlc_param.FILT}.fits', final_wcs=True,
                final_refimage=refim,  driz_cr_corr=True, num_cores=8,
                combine_type=combine_type, skysub=True, build=False,
                skyfile='', skymethod='globalmin+match', skyuser='')

    if recoaddnorthup:

        for cfiles in [[ims, '2020']]:

            if len(cfiles) <= 6:
                combine_type = 'minmed'
            else:
                combine_type = 'imedian'

            from astropy.io import fits

            f = fits.open(refim)
            pix_scale = (f[1].header['CD1_1']**2. +
                         f[1].header['CD1_2']**2.)**0.5 * 3600.
            naxis1 = f[1].header['NAXIS1']
            naxis2 = f[1].header['NAXIS2']

            g = fits.open(refim_align)
            final_ra = g[0].header['RA_TARG']
            final_dec = g[0].header['DEC_TARG']

            astrodrizzle.AstroDrizzle(
                [f for f in cfiles[:]],
                output=f'coadd_upnorth_tweak_{dlc_param.FILT}.fits',
                final_rot=0, final_scale=pix_scale, final_outnx=naxis1,
                final_outny=naxis2, final_ra=final_ra, final_dec=final_dec,
                driz_cr_corr=True, num_cores=64, combine_type=combine_type,
                skysub=True, build=False, skyfile='skyfile.txt',
                skyuser='')


def _make_template_cat(image, extension, filt, threshold, maxobjs=100):
    '''run sextractor and generate numbered set of detections, and reg file'''

    fitscat = image.replace('.fits', '_sex.cat')

    ''' PHOT_APERTURES is a DIAMETER !! '''

    command = f'source-extractor {image}[{extension}] -PHOT_APERTURES 6 '\
              f'-FLAG_IMAGE "" -CATALOG_TYPE FITS_LDAC -DETECT_THRESH '\
              f'{threshold} -DEBLEND_MINCONT 0.001 -PARAMETERS_NAME '\
              f'default.param -FILTER_NAME default.conv '\
              f'-CATALOG_NAME {fitscat}'

    os.system(command)

    from astropy.io import fits as pyfits
    p = pyfits.open(fitscat)
    array = p[2].data

    fname = image.replace('.fits', '_ref_' + str(extension) + '.cat')

    reg = open(fname, 'w')

    data = []
    for i in range(len(array)):
        if array['FLAGS'][i] == 0:
            data.append([array['FLUX_AUTO'][i], [array['X_IMAGE'][i],
                        array['Y_IMAGE'][i], array['FLUX_AUTO'][i]]])

    data.sort(reverse=True)

    added_objects = 0
    for i in range(len(data)):
        if added_objects > maxobjs:
            break
        a = round(data[i][1][0], 6)
        b = round(data[i][1][1], 6)
        c = round(data[i][1][2], 6)
        reg.write(f'{a} {b} {c}\n')
        added_objects += 1
    reg.close()

    if True:
        reg = open(image + '_' + str(extension) + '_IMAGE.reg', 'w')
        reg.write('global color=green dashlist=8 3 width=1 font="helvetica 10 '
                  'normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 '
                  'delete=1 include=1 source=1\nphysical\n')
        for i in range(len(data)):
            a = round(data[i][1][0], 2)
            b = round(data[i][1][1], 2)
            reg.write(f'circle({a},{b},20) # font="times 19" color="green"\n')
        reg.close()


def _make_sextractor_cat1(image, extension, filt, threshold,
                          maxobjs=100, bgChipEstimate=False):

    '''run sextractor and generate numbered set of detections, and reg file'''

    if bgChipEstimate:
        fitscat = image.replace('.fits', '_sex_bg.cat')
    else:
        fitscat = image.replace('.fits', '_sex.cat')

    from astropy.io import fits
    f = fits.open(image)
    back_size = 2. * f[1].data.shape[0]

    ''' PHOT_APERTURES is a DIAMETER !! '''
    import os

    if bgChipEstimate:
        bg_suffix = f' -BACK_SIZE {back_size}'
    else:
        bg_suffix = ' '

    command = f'source-extractor {image}[{extension}] -PHOT_APERTURES 6 '\
              f'-FLAG_IMAGE "" -CATALOG_TYPE FITS_LDAC '\
              f'-DETECT_THRESH {threshold} -DEBLEND_MINCONT 0.001 '\
              f'-PARAMETERS_NAME default.param -FILTER_NAME default.conv '\
              f'-CATALOG_NAME {fitscat}{bg_suffix}'
    import os
    os.system(command)

    from astropy.io import fits as pyfits
    p = pyfits.open(fitscat)
    array = p[2].data

    if bgChipEstimate:
        fname = image.replace('.fits', '_ref_' + str(extension) + '_bg.cat')
    else:
        fname = image.replace('.fits', '_ref_' + str(extension) + '.cat')

    reg = open(fname, 'w')

    data = []
    for i in range(len(array)):
        if array['FLAGS'][i] == 0:
            data.append([array['FLUX_AUTO'][i], [array['X_IMAGE'][i],
                         array['Y_IMAGE'][i], array['FLUX_MAX'][i]]])

    data.sort(reverse=True)

    added_objects = 0
    for i in range(len(data)):
        if added_objects > maxobjs:
            break
        a = round(data[i][1][0], 6)
        b = round(data[i][1][1], 6)
        c = round(data[i][1][2], 6)
        reg.write(f'{a} {b} {c}\n')
        added_objects += 1
    reg.close()

    if True:
        reg = open(image + '_' + str(extension) + '_IMAGE.reg', 'w')
        reg.write('global color=green dashlist=8 3 width=1 '
                  'font="helvetica 10 normal" select=1 highlite=1 dash=0 '
                  'fixed=0 edit=1 move=1 delete=1 include=1 '
                  'source=1\nphysical\n')
        for i in range(len(data)):
            a = round(data[i][1][0], 2)
            b = round(data[i][1][1], 2)
            reg.write(f'circle({a},{b},20) # font="times 19" color="green"\n')
        reg.close()

    return fitscat, fname


def align_sci(dlc_param):
    '''
    Aligns all science images to coadded template image

    Parameters
    ----------
    dlc_param : obj
        Parameter object from prep_directory function
    '''

    ims = []
    for im in dlc_param.IMAGES:
        ims.append(im.loc)

    refim = f'{dlc_param.IMROOT}/Images/coadd_tweak_{dlc_param.FILT}_sci.fits'

    if len(ims) <= 6:
        combine_type = 'minmed'
    else:
        combine_type = 'imedian'

    os.chdir(dlc_param.IM_LOC)

    ''' localmin since we just want this for
        source detection and astrometry'''
    for cfiles in [ims]:
        astrodrizzle.AstroDrizzle(cfiles, driz_cr_corr=True,
                                  output=str('output'),
                                  num_cores=8,
                                  combine_type=combine_type,
                                  final_wcs=True, skysub=True,
                                  skymethod='localmin',
                                  preserve=False, )

    try:
        os.mkdir('IMS_aligned')
    except FileExistsError:
        pass

    if True:
        for cfiles in [ims]:

            cfiles_nocr = [x.replace('_flc.fits', '_crclean.fits') for
                           x in cfiles]

            catName = 'astdriz_catfile.list'

            f = open(catName, 'w')

            fnames = []

            for im in cfiles_nocr:
                ''' manually compute sigma '''

                im = os.path.abspath(im)

                im_short = im.replace('//', '/')

                threshold = 10

                _make_sextractor_cat1(im, 1, dlc_param.FILT, threshold,
                                      maxobjs=1000)
                _make_sextractor_cat1(im, 4, dlc_param.FILT, threshold,
                                      maxobjs=1000)

                a = im_short.replace('.fits', '')
                b = im_short.replace('.fits', '')
                f.write(f'{im_short} {a}_ref_1.cat {b}_ref_4.cat\n')

                fnames.append(im_short)

            f.close()

            tweakreg.TweakReg(fnames[:], catfile=catName, xcol=1,
                              ycol=2, updatehdr=True,  nclip=5,
                              peakmax=50000, sigma=2.5, searchrad=10.0,
                              tolerance=5.0, writecat=True,
                              headerlet=True, attach=False,
                              clobber=True, minobj=-1,
                              fitgeometry='general',  wcsname="TWEAK1",
                              shiftfile=True, refimage=refim,
                              outshifts='shift_file.txt',
                              residplot='no plot', see2dplot=False,
                              interactive=False)

            from astropy.table import Table
            shift_tab = Table.read('shift_file.txt',
                                   format='ascii.no_header',
                                   names=['file', 'dx', 'dy', 'rot',
                                          'scale', 'xrms', 'yrms'])

            formats = ['.2f', '.2f', '.3f', '.5f', '.2f', '.2f']
            for i, col in enumerate(shift_tab.colnames[1:]):
                shift_tab[col].format = formats[i]

            cfiles_nocr = [x.replace('_flc.fits', '_crclean.fits') for
                           x in cfiles]

            from stwcs.wcsutil import headerlet
            for fname, fname_orig in zip(cfiles_nocr, cfiles):

                '''  some of these settings are necessary so that
                        DOLPHOT doesn't choke later'''

                headerlet.apply_headerlet_as_primary(
                    fname_orig, fname.replace('.fits', '_hlet.fits'),
                    attach=False, archive=False)

        for fname in ims:
            shutil.copy(fname, f'./IMS_aligned/')
