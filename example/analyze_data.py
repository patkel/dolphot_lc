import os
import functools

class dolphot_output:

    def __init__(self, file, param_file='dolphot.params'):
        self.file = file
        self.param_file_name = param_file
        self.pandas_parse()

    def pandas_parse(self):
        ''' run on the dolphot output file in the simultaneous folder ''' 

        import pandas, scipy
        columns_raw = [x[:-1] for x in open(self.file + '.columns','r').readlines()]
        trans_short = {'Object X position on reference image (or first image, if no reference)': 'x', 
                'Object Y position on reference image (or first image, if no reference)': 'y'
                }

        columns = []
        fmts = []
        for colName in columns_raw:
            colName = colName.replace(colName.split(' ')[0] + ' ','')
            print("colName = ", colName)

            if False:
                if len(colName.split('(')) > 1:    
                    colName = colName.split('(')[0]
                if colName[-1] == ' ': colName = colName[:-1]
                colName = colName.replace(' ','_').replace(',','').lower()
                colName = colName.replace('instrumental', 'instr').replace('measured', 'meas').replace('normalized','norm').replace('transformed', 'trans')
                print("colName = ", colName)
            if colName in columns: colName += '_prime'
            columns.append(colName)
            fmts.append( scipy.float64 )

        print("columns = ", columns)


        self.columns = columns
        self.tab = pandas.read_table(self.file, names=columns, sep='\s+', mangle_dupe_cols=True)
        self.data = list(zip(self.tab['Object X position on reference image (or first image, if no reference)'], self.tab['Object Y position on reference image (or first image, if no reference)'] ) )
        reg = open(self.file + '.reg', 'w')
        reg.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nphysical\n')
        print("self.data = ", self.data)
        print("list(self.data) = ", list(self.data))
        for i in range(len(list(self.data))):   ###range(len(self.data)):
            reg.write('circle(' + str(self.data[i][0]) + ',' + str(self.data[i][1]) + ',2)' + '\n')
        reg.close()

        ###directory = reduce(lambda x,y: x + '/' + y, self.file.split('/')[:-1] )
        directory = functools.reduce(lambda x,y: x + '/' + y, self.file.split('/')[:-1] )

        ''' read in params and then construct addstars with appropriate counts for each object '''
        param_file = directory + '/' + self.param_file_name #'/dolphot.params'

        print("param_file = ", param_file)
        lines = open( param_file , 'r').readlines()
        key_value_pairs = [a.replace('\n','').split(' = ') for a in lines]
        self.dolphot_params = dict( key_value_pairs )

        self.Nimg = int(float(self.dolphot_params['Nimg']))
        ''' parse the info file to get exposure times, instruments, and chips? '''
        ''' the problem is that if there are multiple instruments with the same filter, you have to figure out which the magnitude is '''

        ''' warmstart option doesn't create new info file '''
        ###info_file = reduce(lambda x,y: x + '/' + y, self.file.split('/')[:-1] ) + '/output.info'
        info_file = functools.reduce(lambda x,y: x + '/' + y, self.file.split('/')[:-1] ) + '/output.info'
        print("info_file = ", info_file)
        lines = open( info_file , 'r').readlines()
        filt_lines = list(filter(lambda x: x[0] == '*', lines))

        print("filt_lines = ", filt_lines)

        d = {}

        import string, re
        import pysynphot as S

        ''' make a list of all images '''
        for instru in ['WFPC2', 'ACS', 'WFC3']:
            print("instrument = ", instru)
            for i in range(len(filt_lines)):   ###range(len(filt_lines)):
                ###if string.find(filt_lines[i],'specific info') != -1:
                if str.find(filt_lines[i], 'specific info') != -1:
                    res = re.split('\s+',filt_lines[i])
                    instru_current = res[1].split('-')[0] 
                else:
                    res = re.split('\s+',filt_lines[i])
                    img_number = int(float(res[2][:-1]))
                    filt = res[3]                    
                    exptime = float(res[5].replace('\n',''))

                    fname = self.dolphot_params['img%d_file' % img_number]
                    cmd = 'gethead %s/%s.fits INSTRUME DETECTOR CCDCHIP' % (directory, fname)
                    import commands
                    print("cmd = ", cmd)

                    output = commands.getoutput(cmd)
                    res = re.split('\s+', output)
                    instru_name = res[0]
                    detector = res[1]
                    print("instru_name = ", instru_name)

                    if instru_name == 'WFPC2':
                        if float(detector) == 1: det_name = 'pc'
                        elif 2 <= float(detector)  <= 4: det_name = 'wfc'
                                                                                                                                                                                   
                        o = S.ObsBandpass('%s,%s,%s' % (instru_current, detector, filt))
                    elif instru_name == 'WFC3':
                        det_name = detector
                        if det_name == 'UVIS':
                            chip = res[2]
                            o = S.ObsBandpass('%s,%s%s,%s' % (instru_current, detector, chip, filt))
                        if det_name == 'IR':
                            print(instru_current, detector, filt)
                            o = S.ObsBandpass('%s,%s,%s' % (instru_current, detector, filt))
                    elif instru_name == 'ACS':                                                                    
                        det_name = detector
                        if det_name == 'WFC':
                            chip = res[2]
                            o = S.ObsBandpass('%s,%s%s,%s' % (instru_current, detector, chip, filt))

                    wave_pivot = o.pivot()
                    d[img_number] = {'filter': filt, 'exptime': exptime, 'instrument': instru_current, 'wave_pivot': wave_pivot, 'detector': detector, 'filename': fname}
                                                                                                                                                                                   
        print("d.keys() = ", d.keys())

        self.d = d

        filts_all = []

        ''' make a list of the filters in each set '''                                                                            
        for instru in ['WFPC2', 'ACS', 'WFC3']:
            pivot_dict = dict([ [a['filter'],a['wave_pivot']] for a in filter(lambda x: x['instrument'] == instru, d.values()) ])
            filts = list(set( pivot_dict.keys() ))
                                                                                                                                  
            filt_list = [ [pivot_dict[a],a] for a in filts]
            filt_list.sort()
                                                                                                                                  
            filts = [a[1] for a in filt_list]
            print("instru, filts = ", instru, filts)

            for filt in filts:
                filts_all.append( [instru, filt] )

        self.filts_all = filts_all

        if self.filts_all == []:
            ###instrument, detector, filt = get_group_info(self.groupnum)
            ###self.filts_all = [ [instru, filt] ]
            self.filts_all = [['WFC3', 'F110W']]


        print(" self.filts_all = ", self.filts_all)

        ''' parse filters.dat file '''
        ###self.zps = dolphot_zeropoints()

    def notsure(self):
        print(self.Nimg, self.dolphot_params)

        import string
        for im_num in range(1,self.Nimg+1):
            fname = self.dolphot_params['img%d_file' % im_num]
            key = 'Measured counts, %s (%f, 1000.0 sec)'
            dict_info = d[im_num] 
                                                                                              
            for col in columns:
                ''' find column with the same filename '''
                if string.find( col, fname ) != -1 and string.find( col, 'Measured counts') != -1: # and string.find( col, 'uncertainty') == -1:
                    zp = zps[ '%s_%s' % (dict_info['instrument'], dict_info['filter']) ] 
                    exptime = dict_info['exptime']
                    print(col, self.tab[col], 10.**((self.tab[col.replace('Measured counts','Instrumental VEGAMAG magnitude')] - zp )/-2.5) * exptime )
                    print(col, fname)
                    print(d[im_num])

        raw_input('here')
                                                                                                                                                 
        ''' so if we want to input the correct number of counts, we do 10.**( (mag - zp) / -2.5 ) * exptime '''
        ''' for each exposure measure the conversion between instrumental magnitude and counts '''
        ''' well it turns out that normalized count rate is just the instrument mag 10**(mag / -2.5 ) ''' 
        ''' then we need to get from magnitudes to counts '''
        ''' zeropoints are in 'filters.dat' file in each instrument directory '''

        ''' figure out exposure time '''

        ''' figure out location on reference image '''
                                                                                                                                                 
        ''' make fake sources '''


    def info_single_object(self, star_x, star_y, objname, diff_suffix, settings, save=True, typecode='obj', wcs=None, verbose=True):
        from scipy import spatial 
        from astropy.io import fits as pyfits
        import pylab                                                                                   
        import string
        print("...here?", verbose)
        if verbose: print("HERE RN", star_x, star_y, objname, diff_suffix)
        ###kdtree = spatial.KDTree(self.data)
        kdtree = spatial.KDTree(list(self.data))
        d, row_index = kdtree.query((star_x,star_y))

        print("here before")
        if verbose: print(d, row_index)

        print(d, row_index)
        print("here ")
        x_pos, y_pos = self.data[row_index]

        ''' extract info for each filter '''

        if verbose: print(self.columns)

        meas_dict = {}

        col_index = 11
        for instru, filt in self.filts_all:

            key = '%s_%s' % (instru,filt)

            meas_dict[key] = {} 

            for i in range(col_index, col_index + 13):
                print(self.columns[i].split(',')[0], self.tab[self.columns[i]][row_index])
                meas_dict[key][self.columns[i].split(',')[0]] = self.tab[self.columns[i]][row_index]                             
            if verbose: print(meas_dict)
            col_index += 13

        ''' now save photometry to database '''

        if verbose: print(meas_dict.keys())

        if verbose: print(self.filts_all)

        for col in self.columns[:11]:

            if True: #string.find(col,'F164N') != -1:
                if verbose: print(col, self.tab[col][row_index])

                for instru_filt in meas_dict:
                    short_col = ' '.join(col.split('(')[0].split(' ')[:3] )
                    if short_col[-1] == ' ': short_col = short_col[:-1]
                    meas_dict[instru_filt][short_col] = self.tab[col][row_index]
        curr_im = None

        self.ims_meas_dict = {}

        if True: #for im_num in range(1,self.Nimg+1):
            key = 'Measured counts, %s (%f, 1000.0 sec)'
            for col in self.columns:
                if True: #string.find( col, fname ) != -1 and (string.find( col, 'Normalized') != -1 or string.find( col, 'Instrumental') != -1 or  string.find( col, 'uncertainty') != -1):
                    if verbose: 
                        print(col, self.tab[col][row_index])

                    if True:
                        if col.find(',') != -1:
                            if col.split(',')[1].find('(') != -1:
                                full_im = col.split(',')[1].split('(')[0].replace(' ','').replace('\n','')
                                im = full_im #.split('.')[0]

                                print(col, full_im, im)

                                if im != curr_im: 
                                    curr_im = im
                                    directory = '/'.join(self.file.split('/')[:-1]) + '/'
                                    image_file = directory + full_im + '.fits'

                                    cmd = 'gethead %s MJD-OBS INSTRUME DETECTOR FILTER FILTER1 EXPTIME FILENAME' % (image_file) 
                                    ###import commands
                                    import subprocess
                                    output = subprocess.getoutput(cmd)

                                    import re
                                    res = re.split('\s+', output)
                                    print(res)
                                    mjd = res[0]
                                    instrument = res[1]
                                    detector = res[2]
                                    filt = res[3]
                                    exptime = res[4]
                                    im_use = res[5].split('.')[0]

                                    im_use = im_use.split('_')[0].upper() + '_' + im_use.split('_')[1].lower()

                                    db_im = {} # search_esa_icarus.imphotsql(self.groupnum, im_use, objname, typecode, 'dolphot', settings)

                                    print(im_use)

                                    db_im['mjd'] = mjd
                                    db_im['x'] = str(x_pos) 
                                    db_im['y'] = str(y_pos) 

                                    if wcs is not None:
                                        ra, dec = wcs.wcs_pix2world(x_pos, y_pos, 1)

                                        db_im['radeg'] = str(ra) 
                                        db_im['decdeg'] = str(dec)

                                    db_im['instrument'] = instrument
                                    db_im['detector'] = detector 
                                    db_im['filter'] = filt 
                                    db_im['exptime'] = exptime
                                    db_im['imagetype'] = diff_suffix

                                col_short = col.split(',')[0]

                                ''' add date, flag, etc, version, groupnum '''
                    else: pass
        self.meas_dict = meas_dict