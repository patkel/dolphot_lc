import os
import drizzlepac
from drizzlepac import astrodrizzle
from drizzlepac import tweakreg
from drizzlepac import tweakback
import shutil
import sewpy
from astropy.io import ascii
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

final_scale_coadd = 0.03
ra_coadd, dec_coadd = 177.397192792, 22.3937313215
nx_coadd, ny_coadd = 5500, 5500

dolphot_params_WFPC2 = {
    #'xform' : '"1 0 0"',
    'UseWCS': 1,
    'raper' : 2,
    'rchi' : 1.5,
    'rsky0' : 8,
    'rsky1' : 20,
    'rpsf' : 10,
    'WFC3IRpsfType': 1,
    #'ref2img' : "20 Value Array",
    }

dolphot_params_WFC3_UVIS = {
    'UseWCS': 1,
    'raper' : 3,
    'rchi' : 2.0,
    'rsky0' : 15,
    'rsky1' : 35,
    'rpsf' : 10,
    'WFC3UVISpsfType': 1,
    #'ref2img' : 20 Value Array,
    #'ref2img' : 20 Value Array,
    }

dolphot_params_ACS_WFC = {
    'raper' : 3,
    'rchi' : 2.0,
    'rsky0' : 15,
    'rsky1' : 35,
    'rpsf' : 10,
    'ACSpsfType': 1,
    #'ref2img' : '20 Value Array',
    #'ref2img' : '20 Value Array',    
    }

dolphot_params_WFC3_IR = {
    #'xform' : '"1 0 0"',
    'UseWCS': 1,
    'raper' : 2,
    'rchi' : 1.5,
    'rsky0' : 8,
    'rsky1' : 20,
    'rpsf' : 10,
    'WFC3IRpsfType': 1,
    #'ref2img' : "20 Value Array",
    }

objCoords = {'Icarus': [495.70, 1004.34, 177.398584, 22.396684], # measured from ACS_WFC_F814W/F814Wglass_small100001_drz_sci.fits using hstphot 
    'Iapyx': [485.99, 1007.85, 177.398671, 22.396713], # measured from DIFFERENCE IMAGE WFC3_IR_F125W/F814Wglass_diff_drz.fits using hstphot 
    'Perdix': [488.82, 1009.72, 177.398646, 22.396729], # measured from DIFFERENCE IMAGE WFC3_IR_F125W/F814Wglass_diff_drz.fits using hstphot 
    'S1': [2635.27, 2978.30, 177.398231, 22.395630], # for big image; measured from F125W 133 
    'S2': [2692.26, 2997.02, 177.397718, 22.395786], # for big image; measured from 134 
    'S3': [2730.87, 2967.26, 177.397370, 22.395538], # for big image; measured from 134 
    'S4': [2682.07, 2925.37, 177.397810, 22.395189], # for big image; measured from 134 
    'SX': [2426.58, 3106.85, 177.400112, 22.396701], # for big image; measured from 194
    'STAR1': [1124.680000,2007.160000,177.422940,22.389870],
    'STAR2': [2829.910000,1996.640000,177.402448,22.389754],
    'STAR3': [3082.560000,4425.110000,177.399412,22.416737],
    'STAR4': [3838.270000,2101.550000,177.390331,22.390919],
    'STAR5': [3976.490000,4263.650000,177.388667,22.414942],
    'STAR6': [4150.350000,2744.150000,177.386580,22.398059],
    'STAR101': [-99,-99,177.420382,22.372373],
    'STAR102': [-99,-99,177.406276,22.384212],
    'STAR103': [-99,-99,177.412897,22.391200],
    'STAR104': [-99,-99,177.415413,22.409081],
    'STAR105': [-99,-99,177.379507,22.411876],
    'STAR106': [-99,-99,177.397187,22.372679],
    'STAR107': [-99,-99,177.387891,22.384909],
    'STAR108': [-99,-99,177.406167,22.387714],
    'STAR109': [-99,-99,177.375904,22.401404],
    'STAR110': [-99,-99,177.385290,22.415757],
    'STAR201': [-99,-99,177.412897,22.391200],
    'STAR202': [-99,-99,177.379507,22.411877],
    'STAR203': [-99,-99,177.387891,22.384909],
    'STAR204': [-99,-99,177.406276,22.384211],
    'STAR205': [-99,-99,177.406168,22.387714],
    'STAR206': [-99,-99,177.415414,22.409081],
    'STAR207': [-99,-99,177.395361,22.410062],
    'STAR208': [-99,-99,177.406947,22.408331],
    'STAR209': [-99,-99,177.402394,22.407392],
    'STAR210': [-99,-99,177.418077,22.417575],
    'STAR211': [-99,-99,177.388402,22.393247],
    'STAR212': [-99,-99,177.388777,22.411933],
    'STAR213': [-99,-99,177.407241,22.410524],
    'STAR214': [-99,-99,177.422670,22.390170],
    'STAR215': [-99,-99,177.408522,22.392026],
    'STAR216': [-99,-99,177.386934,22.395965],
    'STAR217': [-99,-99,177.409056,22.393866],
    'STAR218': [-99,-99,177.407859,22.395983],
    'STAR219': [-99,-99,177.4146934,22.4028456],
    'STAR220': [-99,-99,177.388666,22.414936],
    'STAR221': [-99,-99,177.416606,22.398750],
    'STAR222': [-99,-99,177.411282,22.403750],
    'STAR223': [-99,-99,177.397255,22.387050],
    'STAR224': [-99,-99,177.411654,22.415529],
    'STAR225': [-99,-99,177.385809,22.396707],
    'STAR226': [-99,-99,177.388105,22.390316],
    }

#Image object with useful properties
class Image:
    def __init__(self,loc,instrument,detector,filtergrismprism,typ):
        self.loc = loc
        self.instrument = instrument
        self.detector = detector
        self.filtergrismprism = filtergrismprism
        self.typ = typ
        self.prep_loc = imroot + '/dolphot_prepped/' + loc.split('/')[-1]
        self.name = loc.split('/')[-1].split('.fits')[0].split('_')[0]

#Declair useful global variables and check for consistancy among images
def prep_dir(Im_loc, Ref_image, Dolphot_path, Imroot, Sexpath, Sn_ra_me, Sn_dec_me):
    global im_loc
    im_loc = Im_loc
    global ref_image
    ref_image = Ref_image
    global dolphot_path
    dolphot_path = Dolphot_path
    global imroot
    imroot = Imroot
    global sn_ra_me
    sn_ra_me = Sn_ra_me
    global sn_dec_me
    sn_dec_me = Sn_dec_me
    global ref_image_no_dolphot
    ref_image_no_dolphot = ref_image.replace('.fits','_no_dolphot.fits')
    global ref_cat
    ref_cat = ref_image_no_dolphot.replace('.fits','_ref_1.cat')
    global sexpath
    sexpath = Sexpath
    global IMAGES
    IMAGES = glob_image()

    check_for_consistancy()

    global INST
    INST = IMAGES[0].instrument
    global DETEC
    DETEC = IMAGES[0].detector
    global FILT
    FILT = IMAGES[0].filtergrismprism

    shutil.copyfile(ref_image,ref_image_no_dolphot)

    global MASK
    global CHIPS
    global SUFFIX

    if INST == 'WFPC2':
        MASK = '/wfpc2mask'
        CHIPS = [1,2,3,4]
        SUFFIX = 'c0m'

    if INST == 'WFC3':
        MASK = '/wfc3mask'
        if DETEC == 'UVIS':
            hdulist = fits.open(IMAGES[0].loc)
            if hdulist[1].header['NAXIS1'] < 3000. and  hdulist[1].header['NAXIS2'] < 3000.:
                CHIPS = [1]
            else:
                CHIPS = [1,2]
            SUFFIX = 'flc'

        if DETEC == 'IR':
            CHIPS = [1]
            SUFFIX = 'flt'

    if INST == 'ACS':
        MASK = '/acsmask'
        if DETEC == 'WFC':
            CHIPS = [1,2]
            SUFFIX = 'flc'

#Build IMAGES list with image objects
def glob_image():
    print('Building Image Objects')
    image_list = os.listdir(im_loc)
    N = len(image_list)    
    for i in range(0,N):
        image_list[i] = im_loc+'/'+image_list[i]

    image_details = [0]*N

    for i in range(0,N):
        im = fits.open(image_list[i])[0].header
        inst = im['INSTRUME']
        detec = im['DETECTOR']
        filt = im['FILTER']

        if image_list[i] == ref_image:
            typ = 'ref'
        else:
            typ = 'sci'

        image_details[i] = Image(image_list[i], inst, detec, filt, typ)

    return(image_details)

'''
#Run SExtractor
def make_sextractor_cat(im, extension, threshold):

    fitscat = im.replace('.fits','_sex.cat')

    config = {"DETECT_THRESH":threshold, 
        "PHOT_APERTURES":'6', 
        "FLAG_IMAGE":'""', 
        "CATALOG_TYPE":'FITS_LDAC',
        "DEBLEND_MINCONT":'0.001',
        "PARAMETERS_NAME":imroot +'/default.param',
        "FILTER_NAME":imroot + '/default.conv',
        "CATALOG_NAME":fitscat}

    sew = sewpy.SEW(params=['FLUX_AUTO' ,'X_IMAGE', 'Y_IMAGE', 'FLUX_MAX'],config=config,sexpath=sexpath)
    out = sew(im)
    print(type(out))
    fname = im.replace('.fits','_ref_' + str(extension) + '.cat')
    ascii.write(out['table'],fname,overwrite=True)
    
    with open(fname, 'r') as fil:
        data = fil.read().splitlines(True)
    with open(fname, 'w') as fal:
        fal.writelines(data[1:])
'''

#Run SExtractor
#FLAG
def make_sextractor_cat(image, extension, threshold, maxobjs=50):
    ''' run sextractor and generate numbered set of detections, and reg file ''' 

    fitscat = image.replace('.fits','_sex.cat')

    ''' PHOT_APERTURES is a DIAMETER !! '''   
    command = 'sex ' + image + '[' + str(extension) + '] -PHOT_APERTURES 6 -FLAG_IMAGE "" -CATALOG_TYPE FITS_LDAC -DETECT_THRESH ' + str(threshold) + ' -DEBLEND_MINCONT 0.001 -PARAMETERS_NAME '+ imroot +'/default.param -FILTER_NAME ' + imroot + '/default.conv -CATALOG_NAME ' + fitscat  # -MAG_ZEROPOINT %f' % ZPS[filt]           
    print(command,'<-------------------------------098')
    os.system(command)

    p = fits.open(fitscat)
    #print p[2].columns
    array = p[2].data


    fname = image.replace('.fits','_ref_' + str(extension) + '.cat')

    print(fname)
    reg = open(fname,'w')
    print(fname)

    data = [] 
    for i in range(len(array)):
        if array['FLAGS'][i] == 0:
            data.append( [array['FLUX_AUTO'][i], [array['X_IMAGE'][i], array['Y_IMAGE'][i], array['FLUX_MAX'][i]]] )

    data.sort(reverse=True)

    #print data
    print(len(data))
    #raw_input()

    added_objects = 0
    for i in range(len(data)):
        if added_objects > maxobjs: break
        print(data[i][1])
        reg.write('%.6f %.6f %.6f\n' % (data[i][1][0], data[i][1][1], data[i][1][2]) )
        added_objects += 1
    reg.close()

    reg = open(image + '_' + str(extension) + '_IMAGE.reg','w')    
    print(image + '_' + str(extension) + '_IMAGE.reg')
    reg.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nphysical\n')
    for i in range(len(data)):
        reg.write('circle(%.2f' % data[i][1][0] + ',%.2f' % data[i][1][1] + ',20) # font="times 19" color="green"\n')
    reg.close()
    print('done')

#Check that all images have the same instument, detector, and filter
def check_for_consistancy():
    N = len(IMAGES)
    insts = [0]*N
    detecs = [0]*N
    filts = [0]*N

    for i in range(0,N):
        insts[i] = IMAGES[i].instrument
        detecs[i] = IMAGES[i].detector
        filts[i] = IMAGES[i].filtergrismprism

    for i in range(1,N):
        print(insts[0],detecs[0],filts[0])
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

#Shortcut to print and run a shell command
def shell_command(cmd):
    print('SHELL COMMAND:',cmd)
    os.system(cmd)

#Copy files to dolphot_prepped and run through Dolphot masking, splitgroups, and calsky
def prep_files_for_dolphot(image_directory, r_in, r_out, step, sig_low, sig_high):
    prepped_dir = imroot + image_directory

    #FLAG
    if image_directory == '/dolphot_prepped':
        try:
            os.mkdir(prepped_dir)
        except:
            shutil.rmtree(prepped_dir)
            os.mkdir(prepped_dir)

        for im in IMAGES:
            shutil.copyfile(im.loc,im.prep_loc)

    for im in IMAGES:
        shell_command(dolphot_path + MASK + ' ' + prepped_dir + '/' +'%s_%s.fits' % (im.name, SUFFIX))
        shell_command(dolphot_path + '/splitgroups' + ' ' + prepped_dir + '/' +'%s_%s.fits' % (im.name, SUFFIX))
        for chip in CHIPS:
            shell_command(dolphot_path + '/calcsky' + ' ' + prepped_dir + '/' +'%s_%s.chip%s %s %s %s %s %s' % (im.name, SUFFIX, chip, r_in, r_out, step, sig_low, sig_high))

#Create SExtactor libraries and coadded image
def mk_diff(threshold, fitgeometry, nclip, minobj):
    N = len(IMAGES)
    group = [0]*N
    for i in range(0,N):
        group[i] = IMAGES[i].loc

    #crCleanFirst
    camera = INST + '-' + DETEC

    combine_type = 'imedian' 
    if DETEC != 'IR' and len(group) < 7:
        combine_type = 'iminmed'       

    combine_nhigh = 0
    if combine_type in ['median', 'imedian'] :
        nflt = len(group)
        
        docombine=True
        if nflt==1:
            docombine=False

        print(nflt)

        if nflt <= 3 :
            combine_nhigh=0
        elif camera=='WFC3-IR' :
            # For WFC3-IR set combine_nhigh to 1 or 2 for CRs that slip through
            # the up-the-ramp sampling
            combine_nhigh= (nflt>5) + (nflt>9)
        else :
            # For ACS and UVIS set combine_nhigh to 1, 2, 3, or 4, keeping
            # an odd number of pixels for the median each time
            combine_nhigh =  (1 + nflt%2)*( 1 + 2*(nflt>7)*(1-nflt%2) + (nflt>11)*(nflt%2))

    #This may be supposed to loop
    #FLAG
    os.chdir(imroot+'/Images')
    astrodrizzle.AstroDrizzle(group, driz_cr_corr=True, driz_combine='no', output='OUTPUT', combine_type=combine_type, combine_nhigh=combine_nhigh, median=docombine, blot=docombine)
    os.chdir(imroot)

    #runtweakreg
    group_nocr = [x.replace('.fits','_crclean.fits') for x in group]

    catName = im_loc + '/astdriz_catfile.list'

    f = open(catName, 'w')                                                                                 
    ''' repurposing for custom IMAGEFIND '''
    for im in group_nocr[:]:
        ''' manually compute sigma '''                  

        im_short = im.replace('//','/')

        im = im.replace('_flt', '')

        print(im)

        if INST == 'WFC3' and DETEC == 'IR':
            make_sextractor_cat(im, 1, threshold=threshold)
            f.write('%s %s_ref_1.cat\n' % (im_short.replace('_flt',''), im_short.replace('_flt','').replace('.fits','')) )
        else:
            make_sextractor_cat(im, 1, threshold=threshold)
            make_sextractor_cat(im, 4, threshold=threshold)
            f.write('%s %s_ref_1.cat %s_ref_4.cat\n' % (im_short, im_short.replace('.fits',''), im_short.replace('.fits','')) )

    f.close()

    for x in range(0, len(group_nocr)):
        group_nocr[x] = group_nocr[x].replace('_flt','')

    ''' for reasons unknown, the astrometry became awful '''
    os.chdir(im_loc)
    tweakreg.TweakReg( group_nocr[:],    catfile=catName, xcol=1, ycol=2, fluxcol=3, updatehdr=False,  nclip=5, peakmax=50000, sigma=2.5, searchrad=1.0, writecat=True, headerlet=True, attach=False,  shiftfile=True,  clobber=True, minobj=-1, fitgeometry='rscale', residplot='No plot', see2dplot=False, wcsname="TWEAK1")
    os.chdir(imroot)

    #tweakback
    for fname, fname_orig in zip(group_nocr, group):
        print(fname_orig, fname.replace('.fits','_hlet.fits'))

        '''  some of these settings are necessary so that DOLPHOT doesn't choke later '''
        headerlet.apply_headerlet_as_primary(fname_orig,fname.replace('.fits','_hlet.fits'), attach=False, archive=False)
    
    output_dir = imroot + '/coadd'
    os.mkdir(output_dir)

    #coaddastrom
    if len(group) <=6: 
        combine_type = 'minmed'
    else: 
        combine_type =  'imedian'

    output_filename = '/' + FILT + 'glassastrom'
    nflt = len(group)

    docombine=True
    if nflt==1:
        docombine=False

    if nflt <= 3 :
        combine_nhigh=0
    elif camera=='WFC3-IR' :
        # For WFC3-IR set combine_nhigh to 1 or 2 for CRs that slip through
        # the up-the-ramp sampling
        combine_nhigh= (nflt>5) + (nflt>9)
    else :
        # For ACS and UVIS set combine_nhigh to 1, 2, 3, or 4, keeping
        # an odd number of pixels for the median each time
        combine_nhigh =  (1 + nflt%2)*( 1 + 2*(nflt>7)*(1-nflt%2) + (nflt>11)*(nflt%2))

    astrodrizzle.AstroDrizzle(group[:], output=output_dir + output_filename, skysub=False, driz_cr_corr=True, final_wcs=True, final_refimage=ref_image_no_dolphot, combine_type=combine_type, final_pixfrac=0.75,  combine_nhigh=combine_nhigh, median=docombine, blot=docombine, build=True, static=False)

    #adjustwithcoadd
    coadded_image = output_dir + output_filename + '_drz.fits'

    #FLAG
    minobj = 15                        
    threshold = 20

    catName = 'astdriz_catfile.list'

    f = open(catName, 'w')

    for im in [coadded_image]:
        im_short = im.replace('//','/')
        if INST == 'WFC3' and DETEC == 'IR':
            make_sextractor_cat(im, 1, threshold=threshold)
            f.write('%s %s_ref_1.cat\n' % (im_short, im_short.replace('.fits','')) )
        else:
            make_sextractor_cat(im, 1, threshold=threshold)
            make_sextractor_cat(im, 4, threshold=threshold)
            f.write('%s %s_ref_1.cat %s_ref_4.cat\n' % (im_short, im_short.replace('.fits',''), im_short.replace('.fits','')) )
    f.close()

    tweakback.tweakback( coadded_image, input=group, wcsname='TWEAK', newname='TWEAK75', verbose=True, force=True )

    #coadd_drizpipe

    combine_type = 'imedian'

    if DETEC != 'IR' and len(group) < 7:
        combine_type = 'iminmed'

    drizpipe = imroot + '/imaging_drzpipe'
    os.mkdir(drizpipe)

    #FLAG
    #Look at if statements
    driz_cr = 1

    for final_rot in [0, 90]:
        if final_rot == 0:
            if driz_cr == 1:
                output_filename = FILT + 'glass'
            else:
                output_filename = FILT + 'glass' + '_nocr'
        else:
            if driz_cr == 1:
                output_filename = FILT + 'glass_' + str(final_rot)
            else:
                output_filename = FILT + 'glass_' + 'nocr_' + str(final_rot)

        if final_rot == 0:
            singlesci = True
        else:
            singlesci = False 

        secondDrizzle(fltlist=group, outroot=output_dir, output_filename=output_filename, driz_cr=driz_cr,   driz_cr_snr='5 4.5', ra=ra_coadd, dec=dec_coadd, rot=final_rot, naxis12='%d,%d' % (nx_coadd, ny_coadd), pixfrac=1.00, pixscale=final_scale_coadd, combine_type=combine_type, refimage=ref_image_no_dolphot, singlesci=singlesci, clobber=True, build=False)

#FLAG
#Just copied over with only slight mondification
def secondDrizzle( fltlist='*fl?.fits', outroot='final', output_filename = 'coadd.fits', refimage='',
                   ra=None, dec=None, rot=0, imsize_arcsec=None,
                   naxis12=None, driz_cr=False,  driz_cr_snr='5.0 4.5',
                   singlesci=False, pixscale=None, pixfrac=None,
                   wht_type='IVM', combine_type='imedian',
                   clean=True, clobber=True, verbose=True, build=False, debug=False ) :
    """ 
    Run astrodrizzle on a pile of flt images.
    
    If the user does not specify pixscale, pixfrac, or imsize_arcsec
    then these are set to reasonable defaults for the camera.

    Returns the names of the output sci and wht.fits images.
    """

    hdulist=fits.open(fltlist[0])
    hdr=hdulist[0].header
    hdulist.close()

    # For image sets with fewer than 5 images :
    # if the exposure time in the image set varies by more than a factor of 10
    # then disable CR rejection and wipe out existing CR flags, because the
    # drizzlepac driz_cr step will flag sky noise as CRs.
    if len(fltlist) < 5 :
        etimelist = [ fits.getval(flt,'EXPTIME') for flt in fltlist ]
        if max( etimelist ) / min( etimelist ) > 10 :
            driz_cr = -1

    # define the default astrodrizzle parameters for this camera
    # Note that we fake the number of exposures to be 2, so that we get
    # consistent default pixel scales across all epochs, regardless of the
    # varying number of exposures per epoch.  This can of course be
    # over-ridden by the user specifying pixscale and pixfrac.
    instrument = hdr['INSTRUME']
    detector = hdr['DETECTOR']
    camera = instrument + '-' + detector
    drizpar = getdrizpar( instrument, detector, nexposures=2 )

    if not pixscale : pixscale = drizpar['pixscale']
    if not pixfrac : pixfrac = drizpar['pixfrac']

    # the ra and the dec are the desired ra and dec for the center of the frame
    if ra==None and dec==None and refimage=='' :
        if verbose : print("WARNING: No ra,dec or refimage provided. Using target coordinates of first image.")
        #grab the target ra and dec from the header of the first file
        ra,dec = hdr['RA_TARG'],hdr['DEC_TARG']

    # If we only have one image, skip the median,blot,and driz_cr steps
    docombine=True
    if len(fltlist)==1:
        docombine=False
    combine_nhigh = 0
    if combine_type in ['median', 'imedian'] :
        nflt = len(fltlist)
        if nflt <= 3 :
            combine_nhigh=0
        elif camera=='WFC3-IR' :
            # For WFC3-IR set combine_nhigh to 1 or 2 for CRs that slip through
            # the up-the-ramp sampling
            combine_nhigh= (nflt>5) + (nflt>9)
        else :
            # For ACS and UVIS set combine_nhigh to 1, 2, 3, or 4, keeping
            # an odd number of pixels for the median each time
            combine_nhigh =  (1 + nflt%2)*( 1 + 2*(nflt>7)*(1-nflt%2) + (nflt>11)*(nflt%2))

    if imsize_arcsec is None and naxis12 is None :
        imsize_arcsec = drizpar['imsize_arcsec']
    if naxis12 is not None :
        naxis1 = int(naxis12.split(',')[0])
        naxis2 = int(naxis12.split(',')[1])
    else :
        naxis1 = imsize_arcsec/pixscale
        naxis2 = imsize_arcsec/pixscale
    if driz_cr:
        resetbits = 4096
    else:
        resetbits = 0

    a_dir = os.getcwd()

    #Fix Later
    os.chdir(imroot + '/coadd/')
    outroot = './'

    if True:
        astrodrizzle.AstroDrizzle(
        fltlist, output=outroot + output_filename, runfile=outroot+'_astrodriz.log',
        updatewcs=False, resetbits=resetbits,
        restore=False, preserve=False, overwrite=False, clean=clean,
        median=docombine, blot=docombine,
        driz_cr=(driz_cr>0 and docombine),  
        driz_cr_snr=driz_cr_snr,
        build=build,
        combine_type=combine_type, combine_nhigh=combine_nhigh,
        driz_sep_wcs=True, driz_sep_pixfrac=1.0, driz_sep_scale=pixscale,
        driz_sep_ra=ra, driz_sep_dec=dec, driz_sep_rot=rot,
        driz_sep_bits=drizpar['drizbits'],
        driz_sep_outnx=naxis1, driz_sep_outny=naxis2,
        final_wcs=True, final_pixfrac=pixfrac, final_scale=pixscale,
        final_bits=drizpar['drizbits'],
        final_ra=ra, final_dec=dec, final_rot=rot,
        final_outnx=naxis1, final_outny=naxis2,
        final_wht_type=wht_type)

    if not build:
                    
        if fltlist[0].find('_flc.fits') > 0:                                                                                                                
            drzsfx = '_drc'
        elif fltlist[0].find('_flm.fits') > 0:
            drzsfx = '_drc'
        else:
            drzsfx = '_drz'
        outscifile = outroot+output_filename + drzsfx+'_sci.fits'
        outwhtfile = outroot+output_filename + drzsfx+'_wht.fits'
                                                                                                                                                            
        if (not os.path.isfile(outscifile)) or (not os.path.isfile(outwhtfile)):
            if os.path.isfile(outscifile.replace('drc','drz')):
                os.rename(outscifile.replace('drc','drz'), outscifile)
            if os.path.isfile(outwhtfile.replace('drc','drz')):
                os.rename(outwhtfile.replace('drc','drz'), outwhtfile)
                                                                                                                                                            
        print(outscifile, outwhtfile)

        '''        
        import exceptions
        if (not os.path.isfile(outscifile)) or (not os.path.isfile(outwhtfile)):
            raise exceptions.RuntimeError( "astrodriz.py says : Astrodrizzle did not produce the expected output files %s and %s"%(outscifile,outwhtfile) )
        '''

                                                                                                                                                  
        scrubnans( outscifile ) 
        scrubnans( outwhtfile )
                                                                                                                                                            
        scilist = [ outscifile ]
        whtlist = [ outwhtfile ]
        if singlesci :
            if True: 
                astrodrizzle.AstroDrizzle(
                fltlist, output=outroot, updatewcs=False, resetbits=0,
                restore=False, preserve=False, overwrite=False, clean=False,
                driz_separate=True,
                median=False, blot=False, driz_cr=False, driz_combine=False,
                driz_sep_wcs=True, driz_sep_pixfrac=1.0, driz_sep_scale=pixscale,
                driz_sep_ra=ra, driz_sep_dec=dec, driz_sep_rot=rot,
                driz_sep_bits=drizpar['drizbits'],
                driz_sep_outnx=naxis1, driz_sep_outny=naxis2 )
                                                                                                                                                            
                                                                                                                                                            
            #raw_input('check')            
                                                                                                                                                            
            # give the output single_sci.fits files some more helpful names
            for fltfile in fltlist :
                                                                                                                                                            
                #fltfile = fltfile.lower()
                                                                                                                                                            
                print('fltfile', fltfile)
                                                                                                                                                            
                if fltfile.endswith('_flc.fits') :
                    fltsfx = '_flc.fits'
                elif fltfile.endswith('_flm.fits') :
                    fltsfx = '_flm.fits'
                else:
                    fltsfx = '_flt.fits'
                scifile0 = fltfile.replace(fltsfx,'_single_sci.fits')
                scifile1 = scifile0.replace('_single_sci.fits', '_keep_single_sci.fits') #outroot + '_' + scifile0
                                                                                                                                                            
                print(scifile0, scifile1)
                                                                                                                                                            
                #import exceptions
                #if not os.path.isfile( scifile0 ) :
                    #raise exceptions.RuntimeError('Missing _single_sci.fits file %s'%scifile0)
                os.rename( scifile0, scifile1 )
                whtfile0 = scifile0.replace( '_sci.fits','_wht.fits')
                whtfile1 = scifile1.replace( '_sci.fits','_wht.fits')
                os.rename( whtfile0, whtfile1 )
                scilist.append( scifile1 )
                whtlist.append( whtfile1 )
                if clean :
                    maskfile1 = scifile0.replace( '_single_sci.fits','_sci1_single_mask.fits')
                    if os.path.isfile( maskfile1 ) : os.remove( maskfile1 )
                    maskfile2 = scifile0.replace( '_single_sci.fits','_sci2_single_mask.fits')
                    if os.path.isfile( maskfile2 ) : os.remove( maskfile2 )
                                                                                                                                                            
        bpxlist = []
        for whtfile in whtlist :
            bpxfile = whtfile.replace('_wht','_bpx')
            bpxfile = badpix.zerowht2badpix(
                whtfile, bpxfile, verbose=verbose, clobber=clobber )
            bpxlist.append( bpxfile )
                                                                                                                                                            
        if clean :
            for scifile in scilist :
                ctxfile = scifile.replace( '_sci.fits','_ctx.fits')
                if os.path.isfile( ctxfile ) :
                    os.remove( ctxfile )

        os.chdir(a_dir)
                                                                                                                                                            
        return( scilist, whtlist, bpxlist )

#FLAG
#Just copied over with little modification
#May want to add user entered values
def getdrizpar( instrument, detector, nexposures=None ) :
    """
    return a dict with defaults for a few key astrodrizzle parameters,
    based on the instrument and detector 
    """
    if nexposures is None :
        nexposures = 2 # set a middle-of-the-road pixscale as the default

    if instrument=='WFC3':
        if detector.startswith('IR'): 
            if nexposures == 1 :
                pixscale=0.13
            elif nexposures == 2 :
                pixscale=0.09
            elif nexposures >= 3 :
                pixscale=0.06
            # drizbits DQ flags allowed as OK :
            # 8192 = up-the-ramp CR; 512 = bad flat (blobs); 64 = warm pixel
            return( {'pixscale':pixscale, 'pixfrac':1.0, 'imsize_arcsec':30,
                     'drizbits':'8192,512'} )
        elif detector.startswith('UV'):
            if nexposures == 1 :
                pixscale=0.04
            elif nexposures >= 2 :
                pixscale=0.03
            # drizbits DQ flags allowed as OK :
            # 64 = warm pixel; 32 = hot pix CTE tail
            return( {'pixscale':pixscale, 'pixfrac':1.0, 'imsize_arcsec':30,
                     'drizbits':'32' } )
    elif instrument=='ACS': 
        if detector.startswith('WFC'):
            if nexposures == 1 :
                pixscale=0.05
            elif nexposures == 2 :
                pixscale=0.04
            elif nexposures >= 3 :
                pixscale=0.03
            # drizbits DQ flags allowed as OK :
            # 64 = warm pixel; 32 = hot pix CTE tail
            return( {'pixscale':pixscale, 'pixfrac':1.0, 'imsize_arcsec':30,
                     'drizbits':'32'} )
    else :
        raise RuntimeError('Unknown instrument+detector:  %s %s'%(instrument, detector ) )

#FLAG
#Just copied over with little modification
def scrubnans( filename, fillval=0 ):
    """Locate any pixels in the given fits file that have values of NaN,
    indef, or inf. Replace them all with the given fillval.
    """

    hdulist = fits.open( filename, mode='update' )
    imdata = hdulist[0].data
    ybad, xbad  = where( 1-isfinite( imdata ) )
    imdata[ybad, xbad] = fillval
    hdulist.flush()
    hdulist.close()
    return


def prep_imaging():
    try:
        os.mkdir(imroot + '/imaging')
    except:
        print(imroot+'/imaging already exists')

    all_images = [a.loc for a in IMAGES]

    for im in all_images:
        shutil.copyfile(im, imroot + '/imaging/' + im.split('/')[-1])

def gethead(im,item):
    hdulist = fits.open(im)
    value = hdulist[0].header[item]
    hdulist.close()
    return value


def blot_back(image_directory, r_in, r_out, step, sig_low, sig_high):
    prep_imaging()

    #all_images = [a.loc for a in IMAGES]

    im_drz_blot = imroot + '/coadd/%sglass_drz.fits[sci,1]' % (FILT)
    im_drz_blot = imroot + '/coadd/%sglass_drz_sci.fits' % (FILT)
    diff_dir = imroot + '/diffs/'
    dolphot_prepped_dir = imroot + '/dolphot_prepped/'

    rescale_fac = 1

    try:
        os.mkdir(diff_dir)
    except:
        print(diff_dir + ' already exists')

    for im in IMAGES: #all_images:
        im_to_blot = imroot + '/imaging/%s_%s.fits' % (im.name, SUFFIX)
        im_to_blot_dolphot_prepped = imroot + '/dolphot_prepped/%s_%s.fits' % (im.name, SUFFIX)

        p = fits.open(im_to_blot)
        p_dol_prep = fits.open(im_to_blot_dolphot_prepped)

        im_diff = diff_dir + '%s_%s.fits' % (im.name, SUFFIX)

        for chip in CHIPS:
            outdata = imroot + '/imaging/%s_%s_bgblot_%d.fits' % (im.name, SUFFIX, chip)

            EXPTIME_DRZ = gethead(im_drz_blot,'EXPTIME')
            print(EXPTIME_DRZ)

            EXPTIME_BLT = gethead(im_to_blot,'EXPTIME')
            print(EXPTIME_BLT)

            try:
                os.remove(outdata)
            except:
                print(outdata,'does not exist')

            blotobj = teal.load('ablot')

            print(1./EXPTIME_DRZ*EXPTIME_BLT*rescale_fac)
            print(EXPTIME_DRZ, EXPTIME_BLT, rescale_fac)
            print(im_drz_blot, im_to_blot + '[sci,%d]' % chip, outdata, 1./EXPTIME_DRZ*rescale_fac,blotobj)

            if INST == 'ACS' or (INST == 'WFC3' and DETEC == 'UVIS'):
                blot(im_drz_blot, im_to_blot + '[sci,%d]' % chip, outdata, addsky=False, in_units='cps', out_units='counts', expout=1./EXPTIME_DRZ*EXPTIME_BLT*rescale_fac,configObj=blotobj)
            elif INST == 'WFC3' and DETEC == 'IR':
                blot(im_drz_blot, im_to_blot + '[sci,%d]' % chip, outdata, addsky=False, in_units='cps', out_units='counts', expout=1./EXPTIME_DRZ*rescale_fac,configObj=blotobj)

            conv =  p_dol_prep[1].data / p[1].data

            a = fits.open(outdata)
            print(im_to_blot_dolphot_prepped, outdata)
            print(p_dol_prep[1].data , a[1].data , conv)
            print([a.name for a in p])

            chip_name = '%s_%s.chip%d.fits' % (im.name, SUFFIX.lower(), chip)

            if chip == 1:
                p[1].data = 1. * (p[1].data - a[1].data) #* conv

            elif chip == 2:
                p[4].data = 1. * (p[4].data - a[1].data) #* conv

            shutil.copyfile(dolphot_prepped_dir + chip_name.replace('.fits','.sky.fits'), diff_dir + chip_name.replace('.fits','.sky.fits'))

        print('im diff', im_diff)
        #os.remove(im_diff)
        p.writeto(im_diff)

        print(im_diff, im_drz_blot, im_to_blot + '[sci,1]')

    prep_files_for_dolphot(image_directory, r_in, r_out, step, sig_low, sig_high)

#Not needed right now
def make_difference_ims(groupnum=None, instrument='WFC3', detector='IR', filt='F125W', objname='refsdal', big=True, driz_cr=1, single_sci=False):
    #Subtracts coadded images from templates to make difference images
    #os.chdir( os.environ['icarus'] + '/processing/')

    #FLAG
    #Why 90?
    #changed /coadd/ to /dolphot_prepped/
    im_drz_template = imroot + '/coadd/%sglass_%d_drz_sci.fits' % (FILT, 90)

    fac_template = 1

    #im_drz_template = os.environ[name] + '/coadd_backup/%s_%s_%s/%sglass_small%d_drz_sci.fits' % (instrument, detector, filt, filt, group_template)
    print(im_drz_template, 'im_drz_template')

    if single_sci:

        im_glass = [a.loc for a in IMAGES]

        drizpipe = imroot + '/imaging_drzpipe/'

        im_drizpipe = [drizpipe + a.split('/')[-1].replace('_FLT.fits','_keep_single_sci.fits').replace('_FLC.fits','_keep_single_sci.fits').replace('_FLM.fits','_keep_single_sci.fits') for a in im_glass]

        print(im_drizpipe)

        for im_drz in im_drizpipe:

            im_drz_diff = im_drz.replace('single_sci','diff')
            sub_images(im_drz, im_drz_template, im_drz_diff, fac_template, big)
        
    #FLAG
    #changed /coadd/ to /dolphot_prepped/
    im_drz = imroot + '/coadd/%sglass_%d_drz_sci.fits' % (filt, groupnum)
    im_drz_diff = imroot + '/coadd/%sglass_%d_%s_bigdiff_drz.fits' % (filt, groupnum, objname)

    sub_images(im_drz, im_drz_template, im_drz_diff, fac_template, big)          

def sub_images(im_drz, im_drz_template, im_drz_diff, fac_template, big):
    print(im_drz, 'im_drz', im_drz)

    if im_drz:                                                                                                         
        a = fits.open(im_drz)
        template = fits.open(im_drz_template)

        print('-----------------------------------------------')
        print(im_drz)
        print(im_drz_template)
        print('-----------------------------------------------')

        a.writeto('Tommy1.fits')
        template.writeto('Tommy2.fits')

        a[0].data = a[0].data - fac_template*template[0].data

        try:                                                                                                                        
            os.remove(im_drz_diff)
        except:
            print(im_drz_diff, 'does not exist')
        print('saving %s' % im_drz_diff)
        a.writeto(im_drz_diff)
        print('finished saving')

def dolphot_simultaneous():
    imdir_simultaneous = imroot + '/dolphot_prepped'
    try:
        os.mkdir(imdir_simultaneous)
    except:
        print(imdir_simultaneous + ' already exists')

    ''' GET WAY MORE ASTROMETRIC MATCHES IF USE SAME FILTER !! -- HAS A SUBSTANTIAL EFFECT ON PHOTOMETRY '''

    ref_image_use = ref_image
            
    ''' now set up dolphot parameter files '''
    imgNum = 0
    extra_params = {}
    info_params = {}

    imtype = 'fullarray'
    
    #FLAG
    #There is a cleaner way to do this
    files_imaging_dir = [a.loc for a in IMAGES]

    print(files_imaging_dir)

    imdir = imroot + '/dolphot/'
    imdir_dolphot_prepped = imroot + '/dolphot_prepped/'
    os.mkdir(imdir)

    ''' remove SUFFIX '''        
    files = [a for a in files_imaging_dir]

    if INST == 'WFPC2':                                                                                                 
        CHIPS = [1,2,3,4]
    
        for file in files:
            for chip in CHIPS:

                if imgNum >= 99: 
                    print('cannot have more than 99 images')
                    raise Exception 

                info_dict = {'file': file, 'chip': chip, 'SUFFIX': SUFFIX}

                fname_prepped = imdir_dolphot_prepped + '%(file)s_%(SUFFIX)s.chip%(chip)d' % info_dict
                
                ''' check to see if SN position falls within or near the chip boundary '''

                x,y, x_size, y_size = sky2xy(fname_prepped + '.fits')
                print(x, y, x_size, y_size)
                if 0 < x < x_size and 0 < y < y_size:
                    imgNum += 1
                    info_dict['imgNum'] = imgNum

                    fname_simultaneous = imdir_simultaneous + '%(file)s_%(SUFFIX)s.chip%(chip)d' % info_dict         
                                                                                                                        
                    cmd = 'cp %s.fits %s.fits' % (fname_prepped, fname_simultaneous)
                    print(cmd)
                    if not glob(fname_simultaneous + '.fits'): os.system(cmd)
                                                                                                                        
                    cmd = 'cp %s.sky.fits %s.sky.fits' % (fname_prepped, fname_simultaneous)
                    print(cmd)
                    if not glob(fname_simultaneous + '.sky.fits'): os.system(cmd)
                                                                                                                        
                    if not glob(fname_simultaneous + '.fits') or not glob(fname_simultaneous + '.fits'):
                        raise Exception
                                                                                                                        
                    extra_params['img%(imgNum)d_file' % info_dict ] = '%(file)s_%(SUFFIX)s.chip%(chip)d' % info_dict
                    extra_params['img%(imgNum)d_shift' % info_dict ] = '0 0'
                    extra_params['img%(imgNum)d_xform' % info_dict ] = '1 0 0'


                    ''' record information about original file for use later '''
                    ''' get name of original file ''' 
                    command = 'gethead %s.fits EXPNAME' % fname_simultaneous
                    print(command)
                    namef = subprocess.getoutput(command)
                    orig_crclean = imroot + '/imaging/%s_C0M_crclean.fits' % namef
                    orig = imroot + '/imaging/%s_C0M.fits' % namef
                    orig_dq = imroot + '/imaging/%s_C1M.fits' % namef
                    fname_dolphot_masked = imroot + '/dolphot/%s_C0M.fits' % namef

                    print(orig_crclean, orig, orig_dq)

                    info_params['img%(imgNum)d_orig_crclean' % info_dict] = orig_crclean
                    info_params['img%(imgNum)d_orig' % info_dict] = orig
                    info_params['img%(imgNum)d_orig_dq' % info_dict] = orig_dq
                    info_params['img%(imgNum)d_dolphot_masked' % info_dict] = fname_dolphot_masked 

                    info_params['img%(imgNum)d_instrument' % info_dict] =  INST                            
                    info_params['img%(imgNum)d_detector' % info_dict] = DETEC 
                    info_params['img%(imgNum)d_filt' % info_dict] = FILT 

    elif INST == 'WFC3':

        for file in files:
            if DETEC == 'UVIS': #CHIPS = [1,2]                                                 
                                                                                                    
                hdulist = fits.open(imdir_dolphot_prepped + file + '_FLC.fits')

                ''' this is for subarrays '''                                                    
                if hdulist[1].header['NAXIS1'] < 3000. and  hdulist[1].header['NAXIS2'] < 3000.:
                    CHIPS = [1]
                    imtype = 'subarray'
                else:
                    CHIPS = [1,2]  
                                                                                                    
                                                                                                    
                                                                                                    
            elif DETEC == 'IR': CHIPS = [1]  
                                                                                                                                
            for chip in CHIPS:

                if imgNum >= 99: 
                    print('cannot have more than 99 images')
                    raise Exception 

                info_dict = {'file': file, 'chip': chip, 'SUFFIX': SUFFIX}

                fname_prepped = imdir_dolphot_prepped + info_dict['file'].split('/')[-1].split('_')[0] + '_' + SUFFIX +'.chip' + str(chip)
                print('----------------------------------------------------')
                print('FNP:',fname_prepped)
                print('IDDP:',imdir_dolphot_prepped)
                print('Other:', '%(file)s_%(SUFFIX)s.chip%(chip)d' % info_dict)
            
                print(fname_prepped + '.fits')
                x,y, x_size, y_size = sky2xy(fname_prepped + '.fits')                        
                print(x, y, x_size, y_size)
                if 0 < x < x_size and 0 < y < y_size:
                    imgNum += 1
                    info_dict['imgNum'] = imgNum

                    fname_simultaneous = imdir_simultaneous + '%(file)s_%(SUFFIX)s.chip%(chip)d' % info_dict         
                    cmd = 'cp %s.fits %s.fits' % (fname_prepped, fname_simultaneous)
                    print(cmd)
                    if not glob(fname_simultaneous + '.fits'): os.system(cmd)

                    cmd = 'cp %s.sky.fits %s.sky.fits' % (fname_prepped, fname_simultaneous)
                    print(cmd)
                    #if not glob(fname_simultaneous + '.sky.fits'): os.system(cmd)
                                                                                                                
                    #if not glob(fname_simultaneous + '.fits') or not glob(fname_simultaneous + '.sky.fits'):
                    #    raise Exception
                                                                                                                        
                    #extra_params['img%(imgNum)d_file' % info_dict ] = '%(file)s_%(SUFFIX)s.chip%(chip)d' % info_dict
                    extra_params['img%(imgNum)d_file' % info_dict ] = imdir_dolphot_prepped + info_dict['file'].split('/')[-1].split('_')[0] + '_' + SUFFIX +'.chip' + str(chip)
                    extra_params['img%(imgNum)d_shift' % info_dict ] = '0 0'
                    extra_params['img%(imgNum)d_xform' % info_dict ] = '1 0 0'

                    info_params['img%(imgNum)d_instrument' % info_dict] =  INST                            
                    info_params['img%(imgNum)d_detector' % info_dict] = DETEC 
                    info_params['img%(imgNum)d_filt' % info_dict] = FILT 

                    orig = imroot + '/imaging/%s_%s.fits' % (info_dict['file'], info_dict['SUFFIX'])
                    masked = imroot + '/dolphot/%s_%s.fits' % (info_dict['file'], info_dict['SUFFIX'] )


                    info_params['img%(imgNum)d_orig' % info_dict] =  orig 
                    info_params['img%(imgNum)d_masked' % info_dict] =  masked 


    elif INST == 'ACS': 
        if DETEC == 'WFC': CHIPS = [1,2] 
        print(files)
        for file in files:
            for chip in CHIPS:

                if imgNum >= 99: 
                    print('cannot have more than 99 images')
                    raise Exception 

                info_dict = {'file': file, 'chip': chip, 'SUFFIX': SUFFIX}

                fname_prepped = imdir_dolphot_prepped + '%(file)s_%(SUFFIX)s.chip%(chip)d' % info_dict

                print(fname_prepped + '.fits')
                x,y, x_size, y_size = sky2xy(fname_prepped + '.fits')                        
                print(x, y, x_size, y_size)
                if 0 < x < x_size and 0 < y < y_size:
                    imgNum += 1
                    info_dict['imgNum'] = imgNum

                    fname_simultaneous = imdir_simultaneous + '%(file)s_%(SUFFIX)s.chip%(chip)d' % info_dict         
                    cmd = 'cp %s.fits %s.fits' % (fname_prepped, fname_simultaneous)
                    print(cmd)
                    if not glob(fname_simultaneous + '.fits'): os.system(cmd)

                    cmd = 'cp %s.sky.fits %s.sky.fits' % (fname_prepped, fname_simultaneous)
                    print(cmd)
                    if not glob(fname_simultaneous + '.sky.fits'): os.system(cmd)
                                                                                                                
                    if not glob(fname_simultaneous + '.fits') or not glob(fname_simultaneous + '.sky.fits'):
                        raise Exception
                                                                                                                        
                    extra_params['img%(imgNum)d_file' % info_dict ] = '%(file)s_%(SUFFIX)s.chip%(chip)d' % info_dict
                    extra_params['img%(imgNum)d_shift' % info_dict ] = '0 0'
                    extra_params['img%(imgNum)d_xform' % info_dict ] = '1 0 0'

                    info_params['img%(imgNum)d_instrument' % info_dict] =  INST                            
                    info_params['img%(imgNum)d_detector' % info_dict] = DETEC 
                    info_params['img%(imgNum)d_filt' % info_dict] = FILT 

                    orig = imroot + '/imaging/%s_%s.fits' % (info_dict['file'], info_dict['SUFFIX'])
                    masked = imroot + '/dolphot/%s_%s.fits' % (info_dict['file'], info_dict['SUFFIX'] )

                    info_params['img%(imgNum)d_orig' % info_dict] =  orig 
                    info_params['img%(imgNum)d_masked' % info_dict] =  masked 

    extra_params['Nimg'] = imgNum 

    if True:
        if imtype == 'subarray':
            output_dir = imroot + '/coadd/'
            #FLAG
            #ref_image_use = output_dir + FILT + 'glass_' + str(groupnum) + '_drz.fits'
            ref_image_use = output_dir + FILT + 'glass' + '_drz.fits'
        else:
            pass #ref_image_use = copy(ref_image)

    if imtype == 'subarray':
        if not glob(imdir_simultaneous + ref_image_use.split('/')[-1]):

            os.system('cp %s %s' % (ref_image_use, imdir_simultaneous) )

            os.chdir(imdir_simultaneous)

            cmd = dolphot_path + 'wfc3mask %s' % (ref_image_use.split('/')[-1]) 
            print(cmd)
            os.system(cmd)                                                                                           
            cmd = dolphot_path + 'splitgroups %s' % (ref_image_use.split('/')[-1])
            print(cmd) 
            os.system(cmd)                                                                                           
            for chip in [1]: 
                cmd = dolphot_path + 'calcsky %s.chip%s 15 35 4 2.25 2.00' % (ref_image_use.split('/')[-1].replace('.fits',''), chip) 
                print(cmd) 
                os.system(cmd)                                                                                       
            #os.system('
        #raw_input('here')

    else:
        if not glob(imdir_simultaneous + ref_image_use.replace('.fits','.chip1.fits').split('/')[-1]):      
            os.system('cp %s %s' % (ref_image_use.replace('.fits','.chip1.fits'), imdir_simultaneous) )
            os.system('cp %s %s' % (ref_image_use.replace('.fits','.chip1.sky.fits'), imdir_simultaneous) )

    if True:
        extra_params['img0_file'] = ref_image_use.replace('.fits','.chip1').split('/')[-1]  
        extra_params['img0_RAper'] = '4'
        extra_params['img0_RChi'] = '2.0'
        extra_params['img0_RSky'] = '15 35'
        extra_params['img0_RPSF'] = '15'
        #extra_params['img0_SigFind'] = '5'

    param_file = imdir_simultaneous + '/dolphot.params'

    os.chdir( imdir_simultaneous )
    
    ''' here using the recommended settings for a WFC3 UVIS registration image '''
    mk_param(INST, DETEC, param_file, extra_params, imtype)

    print(param_file)        
    os.system('pwd')        

    print(ref_image_use)

    print(imdir_simultaneous)
    os.system('pwd')

    cmd = dolphot_path + '/dolphot output -pdolphot.params' # photsec="0 1 3000 3000 4000 4000"' #% photsec
    print(cmd)
    os.system(cmd)

    #os.chdir( os.environ['icarus'] + '/processing/')

def sky2xy(img):
    ra = coord.Angle( sn_ra_me, unit=u.hour ) #pylint: disable = no-member
    ra_deg = ra.degree
                                                                                  
    dec = coord.Angle( sn_dec_me, unit=u.degree ) #pylint: disable = no-member
    dec_deg = dec.degree         
                                                                                  
    print(ra_deg, dec_deg)

    header = fits.open(img)[0].header

    naxis1 = header['NAXIS1']
    print(naxis1)
    naxis2 = header['NAXIS2']
    print(naxis2)

    try:
        w = WCS(header, fix=True)                                                      
        print('reference image', ref_image)
        print(sn_ra_me, sn_dec_me)
        print(ra_deg, dec_deg)
        
                                                                                      
        coords = coord.SkyCoord( ra_deg*u.degree, dec_deg*u.degree, equinox='J2000' ) #pylint: disable = no-member
        print(coords.ra, coords.dec)
                                                                                      
        print(w.wcs_world2pix(0., 0., 1))
        snx, sny = w.wcs_world2pix( ra_deg, dec_deg, 1, ra_dec_order=True)
        print('pix coords', snx, sny)
    except:
        cmd = 'sky2xy %s %f %f' % (img, ra_deg, dec_deg)
        output = subprocess.getoutput(cmd)
        print(output)

        #FLAG
        import re
        res = re.split('\s+', output.replace(' (off image)','')) #pylint: disable = anomalous-backslash-in-string
        snx = float(res[-2])
        sny = float(res[-1])

    return snx, sny, naxis1, naxis2 

def mk_param(instrument, detector, param_file, extra_params, imtype):

    #print(dolphot_params_all)
    f = open(param_file,'w')
    from copy import copy
    #for key in dolphot_params_all:
    #    f.write('%s = %s\n' % (key, str(dolphot_params_all[key])))
    if instrument == 'WFPC2':        
        dolphot_params_use = copy(dolphot_params_WFPC2)
    elif instrument == 'WFC3' and detector == 'UVIS':        
        dolphot_params_use = copy(dolphot_params_WFC3_UVIS)
    elif instrument == 'WFC3' and detector == 'IR':        
        dolphot_params_use = copy(dolphot_params_WFC3_IR)
    elif instrument == 'ACS' and detector == 'WFC':        
        dolphot_params_use = copy(dolphot_params_ACS_WFC)

    if imtype == 'subarray':
        dolphot_params_use['UseWCS'] = 2

    for key in dolphot_params_use:
        f.write('%s = %s\n' % (key, str(dolphot_params_use[key])))
    for key in extra_params:
        f.write('%s = %s\n' % (key, str(extra_params[key])))
    f.close()

def dolphot_force(special, apermag, force_same_mag, psfphot):
    #Only force_same_mag and psfphot change
    #special is false
    #apermag is false

    imdir_prepped = directory = imroot + '/dolphot_prepped/'
    imdir_simultaneous = directory = imroot + '/diffs/'

    output_dir = imroot + 'coadd/'
    ref_image_subarray = FILT + 'glass_drz'

    #FLAG
    os.chdir( imdir_simultaneous)

    xytfile = open('xytfile', 'w')

    objs = []

    #FLAG
    objNames = ['S1','S2','S3','S4','SX']

    for obj in objNames: 

        _, _, small_ra, small_dec = objCoords[obj]

        w = WCS(fits.open(ref_image)['SCI']) 

        import astropy.units as u

        import astropy.coordinates as coord
        ra = coord.Angle( small_ra, unit=u.hour) #pylint: disable = no-member
        ra_deg = ra.degree
                                                                                      
        dec = coord.Angle( small_dec, unit=u.degree) #pylint: disable = no-member
        dec_deg = dec.degree         
 

        ''' need to translate '''
        big_x, big_y = w.wcs_world2pix(small_ra, small_dec, 1, ra_dec_order=True) 
        print(big_x, big_y, obj)

        objs.append([obj, big_x, big_y])

        xytfile.write('0 1 %f %f 2 10\n' % (big_x, big_y))
    xytfile.close()

    if apermag:
        cmd = 'dolphot singlestar -pdolphot.params xytfile=xytfile usephot=output PSFPhot=0 Force1=1 SigFind=-99 Force1=1 SigFindMult=1.0 SigFinal=-99' 
    else:
        cmd = 'dolphot singlestar -pdolphot.params xytfile=xytfile usephot=output PSFPhot=%s Force1=1 FitSky=1 SigFind=-99 SigFindMult=1.0 SigFinal=-99 ' % psfphot

        if force_same_mag:
            cmd += ' and ForceSameMag=1'
        else:
            cmd += ' and ForceSameMag=0'

    #FLAG
    os.system(cmd)

    statinfo = os.stat('singlestar')
    #FLAG
    os.chdir(imroot + '/processing/')

    print(imdir_simultaneous +  'singlestar')            
    #a = dolphot_output( imdir_simultaneous +  'singlestar', groupnum)
    for objname,x,y in objs:

        if apermag:
            settings = 'forceaper'
        else:
            if force_same_mag:
                settings = 'force_force_same_mag_%d' % psfphot                   
            else:
                settings = 'force_no_force_same_mag_%d' % psfphot                    
            
        #a.info_single_object( x, y, objname, 'diff', settings, verbose=False)

    os.chdir(imroot + '/processing/')

