import os
from drizzlepac import astrodrizzle
from drizzlepac import tweakreg
from drizzlepac import tweakback
import shutil
from astropy.io import fits
import astropy.coordinates as coord
from astropy import units as u
from astropy.wcs import WCS
import sys
from stwcs.wcsutil import headerlet
from drizzlepac.ablot import blot
from sndrizpipe import badpix
from numpy import where, isfinite
from stsci.tools import teal
import subprocess
from glob import glob


# Image object with useful properties
class Image:
    def __init__(self, loc, instrument, detector, filter, typ):
        self.loc = loc
        self.instrument = instrument
        self.detector = detector
        self.filter = filter
        self.typ = typ
        self.prep_loc = f'{IMROOT}/dolphot_prepped/{loc.split("/")[-1]}'
        self.name = loc.split('/')[-1].split('.fits')[0].split('_')[0]


# Declare useful global variables and check for consistancy among images
def prep_dir(Orig_im_loc, Im_loc, Ref_image, Dolphot_path, Imroot, Sexpath,
             Sn_ra_me, Sn_dec_me, dolphot_params):
    global ORIG_IM_LOC
    ORIG_IM_LOC = Orig_im_loc
    global IM_LOC
    IM_LOC = Im_loc

    Images_Setup()

    global REF_IMAGE
    REF_IMAGE = os.path.split(Ref_image)[1]
    global REF_IMAGE_PATH
    REF_IMAGE_PATH = os.path.split(Ref_image)[0]
    global DOLPHOT_PATH
    DOLPHOT_PATH = Dolphot_path
    global IMROOT
    IMROOT = Imroot
    global SN_RA_ME
    SN_RA_ME = Sn_ra_me
    global SN_DEC_ME
    SN_DEC_ME = Sn_dec_me
    global REF_IMAGE_NO_DOLPHOT
    REF_IMAGE_NO_DOLPHOT = f'{REF_IMAGE_PATH}/'\
                           f'{REF_IMAGE.replace(".fits", "_no_dolphot.fits")}'
    global SEXPATH
    SEXPATH = Sexpath
    global IMAGES
    IMAGES = glob_image()

    global DOLPHOT_PARAMS
    DOLPHOT_PARAMS = dolphot_params

    check_for_consistancy()

    global INST
    INST = IMAGES[0].instrument
    global DETEC
    DETEC = IMAGES[0].detector
    global FILT
    FILT = IMAGES[0].filter

    shutil.copyfile(f'{REF_IMAGE_PATH}/{REF_IMAGE}', REF_IMAGE_NO_DOLPHOT)

    global MASK
    global CHIPS
    global SUFFIX
    global IMTYPE

    IMTYPE = 'fullarray'

    if INST == 'WFPC2':
        MASK = '/wfpc2mask'
        CHIPS = [1, 2, 3, 4]
        SUFFIX = 'c0m'

    if INST == 'WFC3':
        MASK = '/wfc3mask'
        if DETEC == 'UVIS':
            hdulist = fits.open(IMAGES[0].loc)
            if hdulist[1].header['NAXIS1'] < 3000. and \
                    hdulist[1].header['NAXIS2'] < 3000.:
                CHIPS = [1]
                IMTYPE = 'subarray'
            else:
                CHIPS = [1, 2]
            SUFFIX = 'flc'

        if DETEC == 'IR':
            CHIPS = [1]
            SUFFIX = 'flt'

    if INST == 'ACS':
        MASK = '/acsmask'
        if DETEC == 'WFC':
            CHIPS = [1, 2]
            SUFFIX = 'flc'

    global DR_SUFFIX
    DR_SUFFIX = 'drz'
    if SUFFIX == 'flc':
        DR_SUFFIX = 'drc'


def Images_Setup():
    os.mkdir(IM_LOC)

    a = os.listdir(ORIG_IM_LOC)

    b = []

    for x in range(0, len(a)):
        if a[x][-4:] == 'fits':
            b.append(a[x])

    for x in range(0, len(b)):
        shutil.copyfile(f'{ORIG_IM_LOC}/{b[x]}', f'{IM_LOC}/{b[x]}')


# Build IMAGES list with image objects
def glob_image():
    print('Building Image Objects')
    image_list = os.listdir(IM_LOC)
    N = len(image_list)
    for i in range(0, N):
        image_list[i] = IM_LOC+'/'+image_list[i]

    image_details = [0]*N

    for i in range(0, N):
        im = fits.open(image_list[i])[0].header
        inst = im['INSTRUME']
        detec = im['DETECTOR']

        if inst == 'WFC3':
            filt = im['FILTER']
        if inst == 'ACS':
            filt = im['FILTER2']

        if image_list[i] == f'{REF_IMAGE_PATH}/{REF_IMAGE}':
            typ = 'ref'
        else:
            typ = 'sci'

        image_details[i] = Image(image_list[i], inst, detec, filt, typ)

    return(image_details)


# Run SExtractor
def make_sextractor_cat(image, extension, threshold, maxobjs=50):
    '''run sextractor and generate numbered set of detections, and reg file'''

    fitscat = image.replace('.fits', '_sex.cat')

    # FLAG sex vs source-extractor
    ''' PHOT_APERTURES is a DIAMETER !! '''
    command = f'{SEXPATH} {image}[{str(extension)}] -PHOT_APERTURES 6 '\
              f'-FLAG_IMAGE "" -CATALOG_TYPE FITS_LDAC -DETECT_THRESH '\
              f'{str(threshold)} -DEBLEND_MINCONT 0.001 -PARAMETERS_NAME '\
              f'{IMROOT}/default.param -FILTER_NAME {IMROOT}/default.conv '\
              f'-CATALOG_NAME {fitscat}'

    print(command, '<-------------------------------098')
    os.system(command)

    p = fits.open(fitscat)
    array = p[2].data

    fname = image.replace('.fits', f'_ref_{str(extension)}.cat')

    print(fname)
    reg = open(fname, 'w')
    print(fname)

    data = []
    for i in range(len(array)):
        if array['FLAGS'][i] == 0:
            data.append([array['FLUX_AUTO'][i], [array['X_IMAGE'][i],
                         array['Y_IMAGE'][i], array['FLUX_MAX'][i]]])

    data.sort(reverse=True)

    print(len(data))

    added_objects = 0
    for i in range(len(data)):
        if added_objects > maxobjs:
            break
        print(data[i][1])
        reg.write(f'{data[i][1][0]:.6f} {data[i][1][1]:.6f} '
                  f'{data[i][1][2]:.6f}\n')
        added_objects += 1
    reg.close()

    reg = open(f'{image}_{str(extension)}_IMAGE.reg', 'w')
    print(f'{image}_{str(extension)}_IMAGE.reg')
    reg.write('global color=green dashlist=8 3 width=1 font="helvetica 10 '
              'normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 '
              'delete=1 include=1 source=1\nphysical\n')
    for i in range(len(data)):
        reg.write(f'circle({data[i][1][0]:.2f},{data[i][1][1]:.2f},20) # '
                  f'font="times 19" color="green"\n')
    reg.close()
    print('done')


# Check that all images have the same instument, detector, and filter
def check_for_consistancy():
    N = len(IMAGES)
    insts = [0]*N
    detecs = [0]*N
    filts = [0]*N

    for i in range(0, N):
        insts[i] = IMAGES[i].instrument
        detecs[i] = IMAGES[i].detector
        filts[i] = IMAGES[i].filter

    for i in range(1, N):
        print(insts[0], detecs[0], filts[0])
        if insts[0] != insts[i]:
            print('Instuments do not match!')
            sys.exit()
        if detecs[0] != detecs[i]:
            print('Detectors do not match!')
            sys.exit()
        if filts[0] != filts[i]:
            print('Filters do not match!')
            sys.exit()

    print('All instument, detector, and filter header entries match!')


# Shortcut to print and run a shell command
def shell_command(cmd):
    print('SHELL COMMAND:', cmd)
    os.system(cmd)


# Copy files to dolphot_prepped and run through Dolphot masking, splitgroups,
# and calsky
def prep_files_for_dolphot(image_directory, r_in, r_out,
                           step, sig_low, sig_high):

    prepped_dir = f'{IMROOT}{image_directory}'

    # FLAG
    if image_directory == '/dolphot_prepped':
        try:
            os.mkdir(prepped_dir)
        except FileExistsError:
            shutil.rmtree(prepped_dir)
            os.mkdir(prepped_dir)

        for im in IMAGES:
            shutil.copyfile(im.loc, im.prep_loc)

    os.chdir(prepped_dir)

    for im in IMAGES:
        shell_command(f'{DOLPHOT_PATH}{MASK} {im.name}_{SUFFIX}.fits')
        shell_command(f'{DOLPHOT_PATH}/splitgroups {im.name}_{SUFFIX}.fits')
        for chip in CHIPS:
            shell_command(f'{DOLPHOT_PATH}/calcsky '
                          f'{im.name}_{SUFFIX}.chip{chip} {r_in} '
                          f'{r_out} {step} {sig_low} {sig_high}')

    os.chdir(REF_IMAGE_PATH)

    shell_command(f'{DOLPHOT_PATH}/splitgroups {REF_IMAGE}')

    os.chdir(IMROOT)

    # Refresh the files as AstroDrizzle changed the header info
    files = os.listdir(ORIG_IM_LOC)
    for file in files:
        shutil.copyfile(f'{ORIG_IM_LOC}/{file}', f'{IMROOT}/Images/{file}')


# Create SExtactor libraries and coadded image
def mk_diff(threshold, fitgeometry, nclip, minobj, final_scale_coadd,
            ra_coadd, dec_coadd, nx_coadd, ny_coadd):
    N = len(IMAGES)
    group = [0]*N
    for i in range(0, N):
        group[i] = IMAGES[i].loc

    # crCleanFirst
    camera = f'{INST}-{DETEC}'

    # FLAG
    # This is just a default setting
    # Non IR images dont get assigned one
    docombine = True

    combine_type = 'imedian'
    if DETEC != 'IR' and len(group) < 7:
        combine_type = 'iminmed'

    combine_nhigh = 0
    if combine_type in ['median', 'imedian']:
        nflt = len(group)

        docombine = True
        if nflt == 1:
            docombine = False

        print(nflt)

        if nflt <= 3:
            combine_nhigh = 0
        elif camera == 'WFC3-IR':
            # For WFC3-IR set combine_nhigh to 1 or 2 for CRs that slip through
            # the up-the-ramp sampling
            combine_nhigh = (nflt > 5) + (nflt > 9)
        else:
            # For ACS and UVIS set combine_nhigh to 1, 2, 3, or 4, keeping
            # an odd number of pixels for the median each time
            combine_nhigh = (1 + nflt % 2)*(1 + 2*(nflt > 7)*(1 - nflt % 2) +
                                            (nflt > 11)*(nflt % 2))

    # This may be supposed to loop
    # FLAG
    os.chdir(IMROOT+'/Images')
    astrodrizzle.AstroDrizzle(group, driz_cr_corr=True, driz_combine='no',
                              output='OUTPUT', combine_type=combine_type,
                              combine_nhigh=combine_nhigh, median=docombine,
                              blot=docombine)
    os.chdir(IMROOT)

    # runtweakreg
    group_nocr = [x.replace('.fits', '_crclean.fits') for x in group]

    catName = f'{IM_LOC}/astdriz_catfile.list'

    f = open(catName, 'w')
    ''' repurposing for custom IMAGEFIND '''
    for im in group_nocr[:]:
        ''' manually compute sigma '''

        im = im.replace('_flt', '')
        im = im.replace('_flc', '')

        im_short = im.replace('//', '/')

        print(im)

        if INST == 'WFC3' and DETEC == 'IR':
            make_sextractor_cat(im, 1, threshold=threshold)
            f.write(f'{im_short.replace("_flt","")} '
                    f'{im_short.replace("_flt","").replace(".fits","")}'
                    f'_ref_1.cat\n')
        else:
            make_sextractor_cat(im, 1, threshold=threshold)
            make_sextractor_cat(im, 4, threshold=threshold)
            f.write(f'{im_short} {im_short.replace(".fits","")}_ref_1.cat '
                    f'{im_short.replace(".fits","")}_ref_4.cat\n')

    f.close()

    for x in range(0, len(group_nocr)):
        group_nocr[x] = group_nocr[x].replace('_flt', '')
        group_nocr[x] = group_nocr[x].replace('_flc', '')

    ''' for reasons unknown, the astrometry became awful '''
    os.chdir(IM_LOC)
    tweakreg.TweakReg(group_nocr[:], catfile=catName, xcol=1, ycol=2,
                      fluxcol=3, updatehdr=False,  nclip=5, peakmax=50000,
                      sigma=2.5, searchrad=1.0, writecat=True, headerlet=True,
                      attach=False,  shiftfile=True,  clobber=True, minobj=-1,
                      fitgeometry='rscale', residplot='No plot',
                      see2dplot=False, wcsname="TWEAK1")
    os.chdir(IMROOT)

    # tweakback
    for fname, fname_orig in zip(group_nocr, group):
        print(fname_orig, fname.replace('.fits', '_hlet.fits'))

        '''some of these settings are necessary so that DOLPHOT doesn't
           choke later '''
        headerlet.apply_headerlet_as_primary(fname_orig,
                                             fname.replace('.fits',
                                                           '_hlet.fits'),
                                             attach=False, archive=False)

    output_dir = f'{IMROOT}/coadd'
    os.mkdir(output_dir)

    # coaddastrom
    if len(group) <= 6:
        combine_type = 'minmed'
    else:
        combine_type = 'imedian'

    output_filename = f'/{FILT}glassastrom'
    nflt = len(group)

    docombine = True
    if nflt == 1:
        docombine = False

    if nflt <= 3:
        combine_nhigh = 0
    elif camera == 'WFC3-IR':
        # For WFC3-IR set combine_nhigh to 1 or 2 for CRs that slip through
        # the up-the-ramp sampling
        combine_nhigh = (nflt > 5) + (nflt > 9)
    else:
        # For ACS and UVIS set combine_nhigh to 1, 2, 3, or 4, keeping
        # an odd number of pixels for the median each time
        combine_nhigh = (1 + nflt % 2)*(1 + 2*(nflt > 7)*(1 - nflt % 2) +
                                        (nflt > 11)*(nflt % 2))

    astrodrizzle.AstroDrizzle(group[:],
                              output=f'{output_dir}{output_filename}',
                              skysub=False, driz_cr_corr=True, final_wcs=True,
                              final_refimage=REF_IMAGE_NO_DOLPHOT,
                              combine_type=combine_type, final_pixfrac=0.75,
                              combine_nhigh=combine_nhigh, median=docombine,
                              blot=docombine, build=True, static=False)

    # adjustwithcoadd
    coadded_image = f'{output_dir}{output_filename}_{DR_SUFFIX}.fits'

    # FLAG
    minobj = 15
    threshold = 20

    catName = 'astdriz_catfile.list'

    f = open(catName, 'w')

    for im in [coadded_image]:
        im_short = im.replace('//', '/')
        if INST == 'WFC3' and DETEC == 'IR':
            make_sextractor_cat(im, 1, threshold=threshold)
            f.write(f'{im_short} {im_short.replace(".fits","")}_ref_1.cat\n')
        else:
            make_sextractor_cat(im, 1, threshold=threshold)
            make_sextractor_cat(im, 4, threshold=threshold)
            f.write(f'{im_short} {im_short.replace(".fits","")}_ref_1.cat '
                    f'{im_short.replace(".fits","")}_ref_4.cat\n')
    f.close()

    tweakback.tweakback(coadded_image, input=group, wcsname='TWEAK',
                        newname='TWEAK75', verbose=True, force=True)

    combine_type = 'imedian'

    if DETEC != 'IR' and len(group) < 7:
        combine_type = 'iminmed'

    drizpipe = f'{IMROOT}/imaging_drzpipe'
    os.mkdir(drizpipe)

    # FLAG
    # Look at if statements
    driz_cr = 1

    for final_rot in [0, 90]:
        if final_rot == 0:
            if driz_cr == 1:
                output_filename = f'{FILT}glass'
            else:
                output_filename = f'{FILT}glass_nocr'
        else:
            if driz_cr == 1:
                output_filename = f'{FILT}glass_{str(final_rot)}'
            else:
                output_filename = f'{FILT}glass_nocr_{str(final_rot)}'

        if final_rot == 0:
            singlesci = True
        else:
            singlesci = False

        secondDrizzle(fltlist=group, outroot=output_dir,
                      output_filename=output_filename, driz_cr=driz_cr,
                      driz_cr_snr='5 4.5', ra=ra_coadd, dec=dec_coadd,
                      rot=final_rot, naxis12=f'{nx_coadd:d},{ny_coadd:d}',
                      pixfrac=1.00, pixscale=final_scale_coadd,
                      combine_type=combine_type, refimage=REF_IMAGE_NO_DOLPHOT,
                      singlesci=singlesci, clobber=True, build=False)


# FLAG
# Just copied over with only slight mondification
def secondDrizzle(fltlist='*fl?.fits', outroot='final',
                  output_filename='coadd.fits', refimage='',
                  ra=None, dec=None, rot=0, imsize_arcsec=None,
                  naxis12=None, driz_cr=False,  driz_cr_snr='5.0 4.5',
                  singlesci=False, pixscale=None, pixfrac=None,
                  wht_type='IVM', combine_type='imedian',
                  clean=True, clobber=True, verbose=True, build=False,
                  debug=False):

    """
    Run astrodrizzle on a pile of flt images.

    If the user does not specify pixscale, pixfrac, or imsize_arcsec
    then these are set to reasonable defaults for the camera.

    Returns the names of the output sci and wht.fits images.
    """

    hdulist = fits.open(fltlist[0])
    hdr = hdulist[0].header
    hdulist.close()

    # For image sets with fewer than 5 images :
    # if the exposure time in the image set varies by more than a factor of 10
    # then disable CR rejection and wipe out existing CR flags, because the
    # drizzlepac driz_cr step will flag sky noise as CRs.
    if len(fltlist) < 5:
        etimelist = [fits.getval(flt, 'EXPTIME') for flt in fltlist]
        if max(etimelist) / min(etimelist) > 10:
            driz_cr = -1

    # define the default astrodrizzle parameters for this camera
    # Note that we fake the number of exposures to be 2, so that we get
    # consistent default pixel scales across all epochs, regardless of the
    # varying number of exposures per epoch.  This can of course be
    # over-ridden by the user specifying pixscale and pixfrac.
    instrument = hdr['INSTRUME']
    detector = hdr['DETECTOR']
    camera = f'{instrument}-{detector}'
    drizpar = getdrizpar(instrument, detector, nexposures=2)

    if not pixscale:
        pixscale = drizpar['pixscale']
    if not pixfrac:
        pixfrac = drizpar['pixfrac']

    # the ra and the dec are the desired ra and dec for the center of the frame
    if ra is None and dec is None and refimage == '':
        if verbose:
            print("WARNING: No ra,dec or refimage provided. Using target"
                  "coordinates of first image.")
        # grab the target ra and dec from the header of the first file
        ra, dec = hdr['RA_TARG'], hdr['DEC_TARG']

    # If we only have one image, skip the median,blot,and driz_cr steps
    docombine = True
    if len(fltlist) == 1:
        docombine = False
    combine_nhigh = 0
    if combine_type in ['median', 'imedian']:
        nflt = len(fltlist)
        if nflt <= 3:
            combine_nhigh = 0
        elif camera == 'WFC3-IR':
            # For WFC3-IR set combine_nhigh to 1 or 2 for CRs that slip through
            # the up-the-ramp sampling
            combine_nhigh = (nflt > 5) + (nflt > 9)
        else:
            # For ACS and UVIS set combine_nhigh to 1, 2, 3, or 4, keeping
            # an odd number of pixels for the median each time
            combine_nhigh = (1 + nflt % 2)*(1 + 2*(nflt > 7)*(1 - nflt % 2) +
                                            (nflt > 11)*(nflt % 2))

    if imsize_arcsec is None and naxis12 is None:
        imsize_arcsec = drizpar['imsize_arcsec']
    if naxis12 is not None:
        naxis1 = int(naxis12.split(',')[0])
        naxis2 = int(naxis12.split(',')[1])
    else:
        naxis1 = imsize_arcsec/pixscale
        naxis2 = imsize_arcsec/pixscale
    if driz_cr:
        resetbits = 4096
    else:
        resetbits = 0

    a_dir = os.getcwd()

    # FLAG
    # Fix Later
    os.chdir(f'{IMROOT}/coadd/')
    outroot = './'

    if True:
        astrodrizzle.AstroDrizzle(
            fltlist, output=f'{outroot}{output_filename}',
            runfile=f'{outroot}_astrodriz.log', updatewcs=False,
            resetbits=resetbits, restore=False, preserve=False,
            overwrite=False, clean=clean, median=docombine, blot=docombine,
            driz_cr=(driz_cr > 0 and docombine), driz_cr_snr=driz_cr_snr,
            build=build, combine_type=combine_type,
            combine_nhigh=combine_nhigh, driz_sep_wcs=True,
            driz_sep_pixfrac=1.0, driz_sep_scale=pixscale, driz_sep_ra=ra,
            driz_sep_dec=dec, driz_sep_rot=rot,
            driz_sep_bits=drizpar['drizbits'], driz_sep_outnx=naxis1,
            driz_sep_outny=naxis2, final_wcs=True, final_pixfrac=pixfrac,
            final_scale=pixscale, final_bits=drizpar['drizbits'], final_ra=ra,
            final_dec=dec, final_rot=rot, final_outnx=naxis1,
            final_outny=naxis2, final_wht_type=wht_type)

    if not build:

        if fltlist[0].find('_flc.fits') > 0:
            drzsfx = '_drc'
        elif fltlist[0].find('_flm.fits') > 0:
            drzsfx = '_drc'
        else:
            drzsfx = '_drz'
        outscifile = f'{outroot}{output_filename}{drzsfx}_sci.fits'
        outwhtfile = f'{outroot}{output_filename}{drzsfx}_wht.fits'

        if(not os.path.isfile(outscifile)) or (not os.path.isfile(outwhtfile)):
            if os.path.isfile(outscifile.replace('drc', 'drz')):
                os.rename(outscifile.replace('drc', 'drz'), outscifile)
            if os.path.isfile(outwhtfile.replace('drc', 'drz')):
                os.rename(outwhtfile.replace('drc', 'drz'), outwhtfile)

        print(outscifile, outwhtfile)

        scrubnans(outscifile)
        scrubnans(outwhtfile)

        scilist = [outscifile]
        whtlist = [outwhtfile]
        if singlesci:
            if True:
                astrodrizzle.AstroDrizzle(
                    fltlist, output=outroot, updatewcs=False, resetbits=0,
                    restore=False, preserve=False, overwrite=False,
                    clean=False, driz_separate=True, median=False, blot=False,
                    driz_cr=False, driz_combine=False, driz_sep_wcs=True,
                    driz_sep_pixfrac=1.0, driz_sep_scale=pixscale,
                    driz_sep_ra=ra, driz_sep_dec=dec, driz_sep_rot=rot,
                    driz_sep_bits=drizpar['drizbits'], driz_sep_outnx=naxis1,
                    driz_sep_outny=naxis2)

            # FLAG
            # give the output single_sci.fits files some more helpful names
            for fltfile in fltlist:
                print('fltfile', fltfile)

                if fltfile.endswith('_flc.fits'):
                    fltsfx = '_flc.fits'
                elif fltfile.endswith('_flm.fits'):
                    fltsfx = '_flm.fits'
                else:
                    fltsfx = '_flt.fits'
                scifile0 = fltfile.replace(fltsfx, '_single_sci.fits')
                scifile1 = scifile0.replace('_single_sci.fits',
                                            '_keep_single_sci.fits')

                print(scifile0, scifile1)

                os.rename(scifile0, scifile1)
                whtfile0 = scifile0.replace('_sci.fits', '_wht.fits')
                whtfile1 = scifile1.replace('_sci.fits', '_wht.fits')
                os.rename(whtfile0, whtfile1)
                scilist.append(scifile1)
                whtlist.append(whtfile1)
                if clean:
                    maskfile1 = scifile0.replace('_single_sci.fits',
                                                 '_sci1_single_mask.fits')
                    if os.path.isfile(maskfile1):
                        os.remove(maskfile1)
                    maskfile2 = scifile0.replace('_single_sci.fits',
                                                 '_sci2_single_mask.fits')
                    if os.path.isfile(maskfile2):
                        os.remove(maskfile2)

        bpxlist = []
        for whtfile in whtlist:
            bpxfile = whtfile.replace('_wht', '_bpx')
            bpxfile = badpix.zerowht2badpix(
                whtfile, bpxfile, verbose=verbose, clobber=clobber)
            bpxlist.append(bpxfile)

        if clean:
            for scifile in scilist:
                ctxfile = scifile.replace('_sci.fits', '_ctx.fits')
                if os.path.isfile(ctxfile):
                    os.remove(ctxfile)

        os.chdir(a_dir)

        return(scilist, whtlist, bpxlist)


def getdrizpar(instrument, detector, nexposures=None):
    """
    return a dict with defaults for a few key astrodrizzle parameters,
    based on the instrument and detector
    """
    if nexposures is None:
        nexposures = 2  # set a middle-of-the-road pixscale as the default

    if instrument == 'WFC3':
        if detector.startswith('IR'):
            if nexposures == 1:
                pixscale = 0.13
            elif nexposures == 2:
                pixscale = 0.09
            elif nexposures >= 3:
                pixscale = 0.06
            return({'pixscale': pixscale, 'pixfrac': 1.0, 'imsize_arcsec': 30,
                    'drizbits': '8192,512'})
        elif detector.startswith('UV'):
            if nexposures == 1:
                pixscale = 0.04
            elif nexposures >= 2:
                pixscale = 0.03
            return({'pixscale': pixscale, 'pixfrac': 1.0, 'imsize_arcsec': 30,
                    'drizbits': '32'})
    elif instrument == 'ACS':
        if detector.startswith('WFC'):
            if nexposures == 1:
                pixscale = 0.05
            elif nexposures == 2:
                pixscale = 0.04
            elif nexposures >= 3:
                pixscale = 0.03
            return({'pixscale': pixscale, 'pixfrac': 1.0, 'imsize_arcsec': 30,
                    'drizbits': '32'})
    else:
        raise RuntimeError(f'Unknown instrument+detector:  {instrument} '
                           f'{detector}')


# FLAG
# Just copied over with little modification
def scrubnans(filename, fillval=0):
    """Locate any pixels in the given fits file that have values of NaN,
    indef, or inf. Replace them all with the given fillval.
    """

    hdulist = fits.open(filename, mode='update')
    imdata = hdulist[0].data
    ybad, xbad = where(1 - isfinite(imdata))
    imdata[ybad, xbad] = fillval
    hdulist.flush()
    hdulist.close()
    return


def prep_imaging():
    try:
        os.mkdir(f'{IMROOT}/imaging')
    except FileExistsError:
        print(f'{IMROOT}/imaging already exists')

    all_images = [a.loc for a in IMAGES]

    for im in all_images:
        shutil.copyfile(im, f'{IMROOT}/imaging/{im.split("/")[-1]}')


def gethead(im, item):
    hdulist = fits.open(im)
    value = hdulist[0].header[item]
    hdulist.close()
    return value


def blot_back(r_in, r_out, step, sig_low, sig_high):
    prep_imaging()
    im_drz_blot = f'{IMROOT}/coadd/{FILT}glass_{DR_SUFFIX}_sci.fits'

    diff_dir = f'{IMROOT}/diffs/'
    dolphot_prepped_dir = f'{IMROOT}/dolphot_prepped/'

    rescale_fac = 1

    try:
        os.mkdir(diff_dir)
    except FileExistsError:
        print(f'{diff_dir} already exists')

    for im in IMAGES:  # all_images:
        im_to_blot = f'{IMROOT}/imaging/{im.name}_{SUFFIX}.fits'
        im_to_blot_dolphot_prepped = f'{IMROOT}/dolphot_prepped/{im.name}_'\
                                     f'{SUFFIX}.fits'

        p = fits.open(im_to_blot)
        p_dol_prep = fits.open(im_to_blot_dolphot_prepped)

        im_diff = f'{diff_dir}{im.name}_{SUFFIX}.fits'

        for chip in CHIPS:
            outdata = f'{IMROOT}/imaging/{im.name}_{SUFFIX}_bgblot_'\
                      f'{chip:d}.fits'

            EXPTIME_DRZ = gethead(im_drz_blot, 'EXPTIME')
            print(EXPTIME_DRZ)

            EXPTIME_BLT = gethead(im_to_blot, 'EXPTIME')
            print(EXPTIME_BLT)

            try:
                os.remove(outdata)
            except FileNotFoundError:
                print(outdata, 'does not exist')

            blotobj = teal.load('ablot')

            if INST == 'ACS' or (INST == 'WFC3' and DETEC == 'UVIS'):
                blot(im_drz_blot, f'{im_to_blot}[sci,{chip:d}]', outdata,
                     addsky=False, in_units='cps', out_units='counts',
                     expout=1./EXPTIME_DRZ*EXPTIME_BLT*rescale_fac,
                     configObj=blotobj)

            elif INST == 'WFC3' and DETEC == 'IR':
                blot(im_drz_blot, f'im_to_blot[sci,{chip:d}]', outdata,
                     addsky=False, in_units='cps', out_units='counts',
                     expout=1./EXPTIME_DRZ*rescale_fac, configObj=blotobj)

            conv = p_dol_prep[1].data / p[1].data

            a = fits.open(outdata)
            print(im_to_blot_dolphot_prepped, outdata)
            print(p_dol_prep[1].data, a[1].data, conv)
            print([a.name for a in p])

            chip_name = f'{im.name}_{SUFFIX.lower()}.chip{chip:d}.fits'

            if chip == 1:
                p[1].data = 1. * (p[1].data - a[1].data)  # * conv

            elif chip == 2:
                p[4].data = 1. * (p[4].data - a[1].data)  # * conv

            shutil.copyfile(f'{dolphot_prepped_dir}'
                            f'{chip_name.replace(".fits", ".sky.fits")}',
                            f'{diff_dir}'
                            f'{chip_name.replace(".fits", ".sky.fits")}')

        print('im diff', im_diff)
        p.writeto(im_diff)

        print(im_diff, im_drz_blot, f'{im_to_blot}[sci,1]')

    prep_files_for_dolphot('/diffs', r_in, r_out,
                           step, sig_low, sig_high)


def dolphot_simultaneous():
    imdir_simultaneous = f'{IMROOT}/dolphot_prepped'
    try:
        os.mkdir(imdir_simultaneous)
    except FileExistsError:
        print(f'{imdir_simultaneous} already exists')

    ref_image_use = f'{REF_IMAGE}'

    ''' now set up dolphot parameter files '''
    imgNum = 0
    extra_params = {}
    info_params = {}

    files = [a.loc for a in IMAGES]

    imdir = f'{IMROOT}/dolphot/'
    imdir_dolphot_prepped = f'{IMROOT}/dolphot_prepped/'
    os.mkdir(imdir)

    for file in files:
        for chip in CHIPS:

            if imgNum >= 99:
                print('cannot have more than 99 images')
                raise Exception

            name_root = file.split("/")[-1].split("_")[0]

            fname_prepped = f'{imdir_dolphot_prepped}{name_root}_'\
                            f'{SUFFIX}.chip{str(chip)}'

            print(f'{fname_prepped}.fits')
            x, y, x_size, y_size = sky2xy(f'{fname_prepped}.fits')
            print(x, y, x_size, y_size)
            if 0 < x < x_size and 0 < y < y_size:
                imgNum += 1

                fname_simultaneous = f'{imdir_simultaneous}/{name_root}_'\
                                     f'{SUFFIX}.chip{chip:d}'

                cmd = f'cp {fname_prepped}.fits {fname_simultaneous}.fits'
                print(cmd)
                if not glob(f'{fname_simultaneous}.fits'):
                    os.system(cmd)

                cmd = f'cp {fname_prepped}.sky.fits '\
                      f'{fname_simultaneous}.sky.fits'
                print(cmd)

                extra_params[f'img{imgNum}_file'] = f'{name_root}_{SUFFIX}'\
                                                    f'.chip{str(chip)}'
                extra_params[f'img{imgNum}_shift'] = '0 0'
                extra_params[f'img{imgNum}_xform'] = '1 0 0'

                if INST != 'WFPC2':
                    info_params[f'img{imgNum}_instrument'] = INST
                    info_params[f'img{imgNum}_detector'] = DETEC
                    info_params[f'img{imgNum}_filt'] = FILT

                    orig = f'{IMROOT}/imaging/{file}_{SUFFIX}.fits'
                    masked = f'{IMROOT}/dolphot/{file}_{SUFFIX}.fits'

                    info_params[f'img{imgNum}_orig'] = orig
                    info_params[f'img{imgNum}_masked'] = masked

                if INST == 'WFPC2':
                    command = f'gethead {fname_simultaneous}.fits EXPNAME'
                    namef = subprocess.getoutput(command)

                    orig_crclean = f'{IMROOT}/imaging/{namef}_C0M_crclean.fits'
                    orig = f'{IMROOT}/imaging/{namef}_C0M.fits'
                    orig_dq = f'{IMROOT}/imaging/{namef}_C1M.fits'
                    fn_dp_masked = f'{IMROOT}/dolphot/{namef}_C0M.fits'

                    info_params[f'img{imgNum}_orig_crclean'] = orig_crclean
                    info_params[f'img{imgNum}_orig'] = orig
                    info_params[f'img{imgNum}_orig_dq'] = orig_dq
                    info_params[f'img{imgNum}_dolphot_masked'] = fn_dp_masked

                    info_params[f'img{imgNum}_instrument'] = INST
                    info_params[f'img{imgNum}_detector'] = DETEC
                    info_params[f'img{imgNum}_filt'] = FILT

    extra_params['Nimg'] = imgNum

    if IMTYPE == 'subarray':
        output_dir = f'{IMROOT}/coadd/'
        ref_image_use = f'{output_dir}{FILT}glass_drz.fits'

    if IMTYPE == 'subarray':
        if not glob(f'{imdir_simultaneous}{ref_image_use.split("/")[-1]}'):

            os.system(f'cp {ref_image_use} {imdir_simultaneous}')

            os.chdir(imdir_simultaneous)

            cmd = f'{DOLPHOT_PATH}wfc3mask {ref_image_use.split("/")[-1]}'
            print(cmd)
            os.system(cmd)
            cmd = f'{DOLPHOT_PATH}splitgroups {ref_image_use.split("/")[-1]}'
            print(cmd)
            os.system(cmd)
            for chip in [1]:
                cmd = f'{DOLPHOT_PATH}calcsky '\
                      f'{ref_image_use.split("/")[-1].replace(".fits","")}'\
                      f'.chip{chip} 15 35 4 2.25 2.00'
                print(cmd)
                os.system(cmd)

    else:
        refim_c1 = ref_image_use.replace(".fits", ".chip1.fits").split("/")[-1]
        if not glob(f'{imdir_simultaneous}{refim_c1}'):
            os.system(f'cp {ref_image_use.replace(".fits",".chip1.fits")} '
                      f'{imdir_simultaneous}')
            os.system(f'cp {ref_image_use.replace(".fits",".chip1.sky.fits")} '
                      f'{imdir_simultaneous}')

    if True:
        extra_params['img0_file'] = ref_image_use.replace('.fits', '.chip1')
        extra_params['img0_RAper'] = '4'
        extra_params['img0_RChi'] = '2.0'
        extra_params['img0_RSky'] = '15 35'
        extra_params['img0_RPSF'] = '15'

    shutil.copyfile(f'{REF_IMAGE_PATH}/'
                    f'{ref_image_use.replace(".fits",".chip1.fits")}',
                    f'{imdir_simultaneous}/'
                    f'{ref_image_use.replace(".fits",".chip1.fits")}')

    param_file = f'{imdir_simultaneous}/dolphot.params'

    os.chdir(imdir_simultaneous)

    '''
    here using the recommended settings for a WFC3 UVIS registration image
    '''

    mk_param(INST, DETEC, param_file, extra_params, IMTYPE)

    print(param_file)
    os.system('pwd')

    print(ref_image_use)

    print(imdir_simultaneous)
    os.system('pwd')

    cmd = f'{DOLPHOT_PATH}/dolphot output -pdolphot.params'
    print(cmd)
    os.system(cmd)

    print('END')


def sky2xy(img):
    ra = coord.Angle(SN_RA_ME, unit=u.hour)  # pylint: disable = no-member
    ra_deg = ra.degree

    dec = coord.Angle(SN_DEC_ME, unit=u.degree)  # pylint: disable = no-member
    dec_deg = dec.degree

    print(ra_deg, dec_deg)

    header = fits.open(img)[0].header

    naxis1 = header['NAXIS1']
    print(naxis1)
    naxis2 = header['NAXIS2']
    print(naxis2)

    try:
        w = WCS(header, fix=True)
        print('reference image', REF_IMAGE)
        print(SN_RA_ME, SN_DEC_ME)
        print(ra_deg, dec_deg)

        coords = coord.SkyCoord(ra_deg*u.degree, dec_deg*u.degree,
                                equinox='J2000')  # pylint: disable = no-member
        print(coords.ra, coords.dec)

        print(w.wcs_world2pix(0., 0., 1))
        snx, sny = w.wcs_world2pix(ra_deg, dec_deg, 1, ra_dec_order=True)
        print('pix coords', snx, sny)
    except:
        cmd = f'sky2xy {img} {ra_deg:f} {dec_deg:f}'
        output = subprocess.getoutput(cmd)
        print(output)

        # FLAG
        import re
        res = re.split('\\s+', output.replace(' (off image)',
                       ''))  # pylint: disable = anomalous-backslash-in-string
        snx = float(res[-2])
        sny = float(res[-1])

    return snx, sny, naxis1, naxis2


def mk_param(instrument, detector, param_file, extra_params, imtype):
    f = open(param_file, 'w')
    from copy import copy

    if imtype == 'subarray':
        DOLPHOT_PARAMS['UseWCS'] = 2

    for key in DOLPHOT_PARAMS:
        f.write(f'{key} = {str(DOLPHOT_PARAMS[key])}\n')
    for key in extra_params:
        f.write(f'{key} = {str(extra_params[key])}\n')
    f.close()


def dolphot_force(objCoords, apermag=False,
                  force_same_mag=True, psfphot=1):
    imdir_prepped = f'{IMROOT}/dolphot_prepped'
    imdir_simultaneous = f'{IMROOT}/diffs'

    ref_image_subarray = f'{FILT}glass_drz'

    for fname in ['output.columns', 'output.warnings', 'output.info',
                  'output.apcor', 'output.psfs', 'output', 'output.data',
                  'registration.chip1.fits', 'registration.chip1.sky.fits',
                  'output.*.psf.fits', f'{ref_image_subarray}.chip1.fits',
                  f'{ref_image_subarray}.chip1.sky.fits',
                  'F125Wglass_100004_drz.chip1.fits',
                  'F160Wglass_100006_drz.chip1.fits']:
        cmd = f'cp {imdir_prepped}/{fname} {imdir_simultaneous}'
        print('command', cmd)
        os.system(cmd)

    # Start
    imgNum = 0
    extra_params = {}
    info_params = {}

    ref_image_use = REF_IMAGE

    files = [a.loc for a in IMAGES]

    imdir = f'{IMROOT}/dolphot/'
    imdir_dolphot_prepped = f'{IMROOT}/diffs/'

    for file in files:
        for chip in CHIPS:

            if imgNum >= 99:
                print('cannot have more than 99 images')
                raise Exception

            fname_prepped = f'{imdir_dolphot_prepped}'\
                            f'{file.split("/")[-1].split("_")[0]}_{SUFFIX}'\
                            f'.chip{str(chip)}'

            print(f'{fname_prepped}.fits')
            x, y, x_size, y_size = sky2xy(f'{fname_prepped}.fits')
            print(x, y, x_size, y_size)
            if 0 < x < x_size and 0 < y < y_size:
                imgNum += 1

                fname_simultaneous = f'{imdir_simultaneous}{file}_'\
                                     f'{SUFFIX}.chip{chip}'
                cmd = f'cp {fname_prepped}.fits {fname_simultaneous}.fits'
                print(cmd)
                if not glob(f'{fname_simultaneous}.fits'):
                    os.system(cmd)

                cmd = f'cp {fname_prepped}.sky.fits '\
                      f'{fname_simultaneous}.sky.fits'
                print(cmd)

                extra_params[f'img{imgNum}_file'] = \
                    f'{file.split("/")[-1].split("_")[0]}_'\
                    f'{SUFFIX}.chip{str(chip)}'
                extra_params[f'img{imgNum}_shift'] = '0 0'
                extra_params[f'img{imgNum}_xform'] = '1 0 0'

                if INST != 'WFPC2':
                    info_params[f'img{imgNum}_instrument'] = INST
                    info_params[f'img{imgNum}_detector'] = DETEC
                    info_params[f'img{imgNum}_filt'] = FILT

                    orig = f'{IMROOT}/imaging/{file}_{SUFFIX}.fits'
                    masked = f'{IMROOT}/dolphot/{file}_{SUFFIX}.fits'

                    info_params[f'img{imgNum}_orig'] = orig
                    info_params[f'img{imgNum}_masked'] = masked

                if INST == 'WFPC2':
                    command = f'gethead {fname_simultaneous}.fits EXPNAME'
                    print(command)
                    namef = subprocess.getoutput(command)
                    orig_crclean = f'{IMROOT}/imaging/{namef}_C0M_crclean.fits'
                    orig = f'{IMROOT}/imaging/{namef}_C0M.fits'
                    orig_dq = f'{IMROOT}/imaging/{namef}_C1M.fits'
                    fn_dp_masked = f'{IMROOT}/dolphot/{namef}_C0M.fits'

                    info_params[f'img{imgNum}_orig_crclean'] = orig_crclean
                    info_params[f'img{imgNum}_orig'] = orig
                    info_params[f'img{imgNum}_orig_dq'] = orig_dq
                    info_params[f'img{imgNum}_dolphot_masked'] = fn_dp_masked

                    info_params[f'img{imgNum}_instrument'] = INST
                    info_params[f'img{imgNum}_detector'] = DETEC
                    info_params[f'img{imgNum}_filt'] = FILT

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
    mk_param(INST, DETEC, param_file, extra_params, IMTYPE)
    # End

    print(f'directory: {imdir_simultaneous}')

    xytfile = open('xytfile', 'w')

    objs = []

    for key in objCoords.keys():

        _, _, small_ra, small_dec = objCoords[key]

        from astropy.wcs import WCS

        from astropy.io import fits

        w = WCS(fits.open(f'{REF_IMAGE_PATH}/{REF_IMAGE}')['SCI'])

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
        print(big_x, big_y, key)

        objs.append([key, big_x, big_y])

        xytfile.write(f'0 1 {big_x} {big_y} 2 10\n')
    xytfile.close()

    if apermag:
        cmd = f'{DOLPHOT_PATH}/dolphot singlestar -pdolphot.params '\
              f'xytfile=xytfile usephot=output PSFPhot=0 Force1=1 SigFind=99 '\
              f'Force1=1 SigFindMult=1.0 SigFinal=99'
    else:
        cmd = f'{DOLPHOT_PATH}/dolphot singlestar -pdolphot.params '\
              f'xytfile=xytfile usephot=output PSFPhot={psfphot} Force1=1 '\
              f'FitSky=1 SigFind=99 SigFindMult=1.0 SigFinal=99'

        if force_same_mag:
            cmd += ' ForceSameMag=1'
        else:
            cmd += ' ForceSameMag=0'

    print(os.getcwd())
    print(cmd)
    os.system(cmd)

    statinfo = os.stat('singlestar')
