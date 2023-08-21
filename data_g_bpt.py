
from __future__ import print_function
from astropy.io import fits

try:
    import pyfits
except:
    from astropy.io import fits as pyfits
from scipy import ndimage

import glob

from time import clock
from astropy.table import vstack, Table
import multiprocessing as mp
import xlrd
import os

from astropy.io import fits
from matplotlib import pyplot
from matplotlib import image
import matplotlib as mpl
from matplotlib import pyplot as plt

import numpy
from scipy import integrate
import math
import numpy as np

import time
from scipy.optimize import curve_fit
from ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.miles_util as lib
import ppxf as ppxf_package
from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.miles_util as lib
import astropy.units as u
import speclite.filters
from astropy.cosmology import WMAP9 as cosmo

a=xlrd.open_workbook('/pegasus/czh276/age_re/drp-v3_1_1.xls')# get the target names from excel
drp=fits.open('/pegasus/czh276/age_re/drpall-v3_1_1.fits') #drpall file
table=a.sheets()[0]
HYBpath="/pegasus/manga/dap/MPL-11/HYB10-MILESHC-MASTARHC2/"
VORpath="/pegasus/manga/dap/MPL-11/VOR10-MILESHC-MASTARHC2/"
SPXpath="/pegasus/manga/dap/MPL-11/SPX-MILESHC-MASTARHC2/"

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def channel_dictionary(hdu, ext):
    """
    Construct a dictionary of the channels in a MAPS file.
    """
    channel_dict = {}
    for k, v in hdu[ext].header.items():
        if k[0] == 'C':
            try:
                i = int(k[1:])-1
            except ValueError:
                continue
            channel_dict[v] = i
    return channel_dict


def channel_units(hdu, ext):
    """
    Construct an array with the channel units.
    """
    nchannels = 1 if len(hdu[ext].data.shape) == 2 else hdu[ext].data.shape[0]
    channel_units = numpy.empty(nchannels, dtype=object)
    for k, v in hdu[ext].header.items():
        if k[0] == 'U':
            try:
                i = int(k[1:])-1
            except ValueError:
                continue
            channel_units[i] = v.strip()
    return channel_units


def apply_index_dispersion_correction(indx, indxcorr, unit):
    """
    Apply a set of dispersion corrections.
    """
    if unit not in [ 'ang', 'mag' ]:
        raise ValueError('Unit must be mag or ang.')
    return indx * indxcorr if unit == 'ang' else indx + indxcorr
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def feeddata():
    L=len(table.col_values(0))
    out=[[],[],[],[],[],[],[],[],[],[],[]]
    m=0
    k=1
    print(L)
    while m<L:
        print(m)
        Tar=int(table.col_values(3)[m])
        if os.path.exists(HYBpath+str(drp[1].data[Tar][0])+'/'+str(drp[1].data[Tar][1])+'/'+'manga-'+drp[1].data[Tar][2]+'-MAPS-HYB10-MILESHC-MASTARHC2.fits.gz')==False:
            print('HYB'+drp[1].data[Tar][2])
            
            m+=1
            continue
        
        re=drp[1].data[Tar][87]#in arcsecond

        re=re/0.5 #half long axis in spaxels
        if table.col_values(2)[m]!='P' or re<=0:
            print(table.col_values(2)[m],re)
            m+=1
            continue
        
        hybmap=fits.open(HYBpath+str(drp[1].data[Tar][0])+'/'+str(drp[1].data[Tar][1])+'/'+'manga-'+drp[1].data[Tar][2]+'-MAPS-HYB10-MILESHC-MASTARHC2.fits.gz')
        
        hdu=hybmap.copy()
        emlc = channel_dictionary(hdu, 'EMLINE_GFLUX')
        mask_ext = hdu['EMLINE_GFLUX'].header['QUALDATA']
        ha_flux = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data[emlc['Ha-6564'],:,:].copy(), mask=hdu[mask_ext].data[emlc['Ha-6564'],:,:].copy() > 0)
        ha_flux_IVAR= numpy.ma.MaskedArray(hdu['EMLINE_GFLUX_IVAR'].data[emlc['Ha-6564'],:,:].copy(), mask=hdu[mask_ext].data[emlc['Ha-6564'],:,:].copy() > 0)
        ha_flux_err=1/numpy.sqrt(ha_flux_IVAR)
        OII1_flux = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data[emlc['OII-3727'],:,:].copy(), mask=hdu[mask_ext].data[emlc['OII-3727'],:,:].copy() > 0)
        OII1_flux_IVAR = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX_IVAR'].data[emlc['OII-3727'],:,:].copy(), mask=hdu[mask_ext].data[emlc['OII-3727'],:,:].copy() > 0)
        OII1_flux_err=1/numpy.sqrt(OII1_flux_IVAR)
        OII2_flux = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data[emlc['OII-3729'],:,:].copy(), mask=hdu[mask_ext].data[emlc['OII-3729'],:,:].copy() > 0)
        OII2_flux_IVAR = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX_IVAR'].data[emlc['OII-3729'],:,:].copy(), mask=hdu[mask_ext].data[emlc['OII-3729'],:,:].copy() > 0)
        OII2_flux_err=1/numpy.sqrt(OII2_flux_IVAR)
        OII_flux = OII1_flux + OII2_flux
        OII_flux_err = OII1_flux_err + OII2_flux_err
        hb_flux = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data[emlc['Hb-4862'],:,:].copy(), mask=hdu[mask_ext].data[emlc['Hb-4862'],:,:].copy() > 0)
        hb_flux_IVAR= numpy.ma.MaskedArray(hdu['EMLINE_GFLUX_IVAR'].data[emlc['Hb-4862'],:,:].copy(), mask=hdu[mask_ext].data[emlc['Hb-4862'],:,:].copy() > 0)
        hb_flux_err=1/numpy.sqrt(ha_flux_IVAR)

        OIII_flux = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data[emlc['OIII-5008'],:,:].copy(), mask=hdu[mask_ext].data[emlc['OIII-5008'],:,:].copy() > 0)
        OIII_flux_IVAR= numpy.ma.MaskedArray(hdu['EMLINE_GFLUX_IVAR'].data[emlc['OIII-5008'],:,:].copy(), mask=hdu[mask_ext].data[emlc['OIII-5008'],:,:].copy() > 0)
        OIII_flux_err=1/numpy.sqrt(OIII_flux_IVAR)
        
        SII1_flux = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data[emlc['SII-6718'],:,:].copy(), mask=hdu[mask_ext].data[emlc['SII-6718'],:,:].copy() > 0)
        SII1_flux_IVAR = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX_IVAR'].data[emlc['SII-6718'],:,:].copy(), mask=hdu[mask_ext].data[emlc['SII-6718'],:,:].copy() > 0)
        SII1_flux_err=1/numpy.sqrt(SII1_flux_IVAR)
        SII2_flux = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data[emlc['SII-6732'],:,:].copy(), mask=hdu[mask_ext].data[emlc['SII-6732'],:,:].copy() > 0)
        SII2_flux_IVAR = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX_IVAR'].data[emlc['SII-6732'],:,:].copy(), mask=hdu[mask_ext].data[emlc['SII-6732'],:,:].copy() > 0)
        SII2_flux_err=1/numpy.sqrt(SII2_flux_IVAR)
        SII_flux = SII1_flux + SII2_flux
        SII_flux_err = SII1_flux_err + SII2_flux_err
        
        NII_flux = numpy.ma.MaskedArray(hdu['EMLINE_GFLUX'].data[emlc['NII-6585'],:,:].copy(), mask=hdu[mask_ext].data[emlc['NII-6585'],:,:].copy() > 0)
        NII_flux_IVAR= numpy.ma.MaskedArray(hdu['EMLINE_GFLUX_IVAR'].data[emlc['NII-6585'],:,:].copy(), mask=hdu[mask_ext].data[emlc['NII-6585'],:,:].copy() > 0)
        NII_flux_err=1/numpy.sqrt(OIII_flux_IVAR)
        
        #MASK=fractional(l,Tar,0.2,float('nan'))
        #SPEC=avespec(spx,vormap,MASK,redshift)
        '''if SPEC=='error':
            print('noise invalid!')
            m+=1
            continue'''
        out[0].append(Tar)
        out[1].append(ha_flux)
        out[2].append(ha_flux_err)
        out[3].append(hb_flux)
        out[4].append(hb_flux_err)
        out[5].append(OIII_flux)
        out[6].append(OIII_flux_err)
        out[7].append(SII_flux)
        out[8].append(SII_flux_err)
        out[9].append(NII_flux)
        out[10].append(NII_flux_err)
        k+=1
        #HDU.close()
        hybmap.close()
        
        #spx.close()
        m+=1
    return(np.transpose(numpy.array(out)))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def output(out):
    
    total=0
    o2=0
    ha=0
    good=0
    l1=0
    l2=0
    L=len(out)
    print(L,'*')
    i=0
    while i<L:
        l1+=out[i][0]
        l2+=out[i][1]
        i+=1
    
    print(l1,l2,l2/float(l1))
    
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def run(data):
    L=len(data)
    dat=[]
    
    i=0
    while i<L:
        print(i,'of',L,'\n')
        Tar=int(data[i][0])
        plate_ifu=drp[1].data[Tar][2]
        ha_flux=data[i][1].data.flatten(order='C')
        ha_flux_err=data[i][2].data.flatten(order='C')
        ha_flux_mask=data[i][1].mask.astype(int).flatten(order='C')
        ha_ew=data[i][3].data.flatten(order='C')
        ha_ew_err=data[i][4].data.flatten(order='C')
        ha_ew_mask=data[i][3].mask.astype(int).flatten(order='C')
        OII_flux=data[i][5].data.flatten(order='C')
        OII_flux_err=data[i][6].data.flatten(order='C')
        OII_flux_mask=data[i][5].mask.astype(int).flatten(order='C')
        OII_ew=data[i][7].data.flatten(order='C')
        OII_ew_err=data[i][8].data.flatten(order='C')
        OII_ew_mask=data[i][7].mask.astype(int).flatten(order='C')
        NII_flux=data[i][9].data.flatten(order='C')
        NII_flux_err=data[i][10].data.flatten(order='C')
        NII_flux_mask=data[i][9].mask.astype(int).flatten(order='C')
        p_i=np.array([plate_ifu for j in range(0,ha_flux.shape[0])])
        subdat=Table()
        subdat['ha_flux']=ha_flux
        subdat['ha_flux_err']=ha_flux_err
        subdat['ha_flux_mask']=ha_flux_mask
        subdat['hb_flux']=ha_ew
        subdat['hb_flux_err']=ha_ew_err
        subdat['hb_flux_mask']=ha_ew_mask
        subdat['o3_flux']=OII_flux
        subdat['o3_flux_err']=OII_flux_err
        subdat['o3_flux_mask']=OII_flux_mask
        subdat['s2_flux']=OII_ew
        subdat['s2_flux_err']=OII_ew_err
        subdat['s2_flux_mask']=OII_ew_mask
        subdat['n2_flux']=NII_flux
        subdat['n2_flux_err']=NII_flux_err
        subdat['n2_flux_mask']=NII_flux_mask
        subdat['plate_ifu']=p_i
        dat.append(subdat)
        i+=1
        
    print(len(dat))
    datall=vstack(dat)
    datall.write('./data_g_bpt.fits')
    return(0)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
flag1=time.time()
#book=open('0711p1_3.txt','w')
data=feeddata()
flag4=time.time()
print(len(data))


print('start output')
run(data)
print('end output')
flag2=time.time()
Time=flag4-flag1
m=int(Time/60.)
s=Time-m*60
#out1=run(data)
print('load data       ',m,'min,',s,'s')
Time=flag2-flag4
m=int(Time/60.)
s=Time-m*60
#out1=run(data)
print('output data    ',m,'min,',s,'s')

'''
out=[]
i=0
while i< len(data):
    out.append(task(data[i]))
    i+=1
'''
'''data1=[]
L=len(out1)
i=0
while i<L:
    if out1[2]!='invalid1' and:
        data1.append(out1[i])
    i+=1
'''

flag3=time.time()


Time=flag3-flag1
m=int(Time/60.)
s=Time-m*60
#out1=run(data)
print('Total runn time' ,m,'min,',s,'s')
print('\n data size='+str(len(data))+'\n')
