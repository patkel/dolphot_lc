from drizzlepac import astrodrizzle
import os
import sys
import shutil

def download(items):
    for thing in items:
        os.system(f'wget {thing}')
        os.rename(thing.split('/Download/')[1].replace('/','%2F'), thing.split('/')[-1])

items = ['https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/jcdu51vrq_flc.fits',
         'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/jcdu51vtq_flc.fits',
         'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/jcdu51vxq_flc.fits',
         'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/jcdu51w0q_flc.fits']

cwd = os.getcwd()
print(cwd)

os.mkdir('ref')
os.mkdir('ref/ims')

combine_type = 'minmed'
ims_dir = f'{cwd}/ref/ims'

os.chdir(ims_dir)

download(items)

files = os.listdir(ims_dir)
cfiles = []

for file in files:
    cfiles.append(f'{ims_dir}/{file}')

os.chdir(f'{cwd}/ref')

astrodrizzle.AstroDrizzle(
    [f for f in cfiles], output=str('output'), final_wcs=True,
    driz_cr_corr=True, num_cores=48,
    combine_type=combine_type, skysub=False,
    skymethod='localmin', skystat='mode', build=True,
    preserve=False, final_rot=0)
    
shutil.rmtree(ims_dir)
os.rename('output_drc.fits', 'registration.fits')

files = os.listdir()
for file in files:
    if file != 'registration.fits':
        os.remove(file)

