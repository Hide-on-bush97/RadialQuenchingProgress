
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
def FIND(name):
    i=0
    while i<6779:
        if drp[1].data[i][2]==name:
            return i
            break
        i+=1
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
    out=[[],[],[],[],[],[],[],[],[]]
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
            print(table.col_values(2)[m])
            m+=1
            continue
        
        hybmap=fits.open(HYBpath+str(drp[1].data[Tar][0])+'/'+str(drp[1].data[Tar][1])+'/'+'manga-'+drp[1].data[Tar][2]+'-MAPS-HYB10-MILESHC-MASTARHC2.fits.gz')
        
        hdu=hybmap.copy()
        emlc = channel_dictionary(hdu, 'EMLINE_SFLUX')
        mask_ext = hdu['EMLINE_SFLUX'].header['QUALDATA']
        ha_flux = numpy.ma.MaskedArray(hdu['EMLINE_SFLUX'].data[emlc['Ha-6564'],:,:].copy(), mask=hdu[mask_ext].data[emlc['Ha-6564'],:,:].copy() > 0)
        ha_flux_IVAR= numpy.ma.MaskedArray(hdu['EMLINE_SFLUX_IVAR'].data[emlc['Ha-6564'],:,:].copy(), mask=hdu[mask_ext].data[emlc['Ha-6564'],:,:].copy() > 0)
        ha_flux_err=1/numpy.sqrt(ha_flux_IVAR)
        OII1_flux = numpy.ma.MaskedArray(hdu['EMLINE_SFLUX'].data[emlc['OIId-3728'],:,:].copy(), mask=hdu[mask_ext].data[emlc['OIId-3728'],:,:].copy() > 0)
        OII1_flux_IVAR = numpy.ma.MaskedArray(hdu['EMLINE_SFLUX_IVAR'].data[emlc['OIId-3728'],:,:].copy(), mask=hdu[mask_ext].data[emlc['OIId-3728'],:,:].copy() > 0)
        OII1_flux_err=1/numpy.sqrt(OII1_flux_IVAR)
        OII_flux = OII1_flux
        OII_flux_err = OII1_flux_err
        mask_ext = hdu['EMLINE_SEW'].header['QUALDATA']
        ha_ew_unmasked=hdu['EMLINE_SEW'].data[emlc['Ha-6564'],:,:].copy()
        #print(ha_ew_unmasked)
        l1=len(np.where(ha_ew_unmasked!=0)[0])
        
        ha_ew = numpy.ma.MaskedArray(hdu['EMLINE_SEW'].data[emlc['Ha-6564'],:,:].copy(), mask=hdu[mask_ext].data[emlc['Ha-6564'],:,:].copy() > 0)
        #print(ha_ew)
        l2=len(np.where(ha_ew>0)[0])
        
        ha_ew_IVAR= numpy.ma.MaskedArray(hdu['EMLINE_SEW_IVAR'].data[emlc['Ha-6564'],:,:].copy(), mask=hdu[mask_ext].data[emlc['Ha-6564'],:,:].copy() > 0)
        ha_ew_err=1/numpy.sqrt(ha_ew_IVAR)
        
        OII1_ew = numpy.ma.MaskedArray(hdu['EMLINE_SEW'].data[emlc['OIId-3728'],:,:].copy(), mask=hdu[mask_ext].data[emlc['OIId-3728'],:,:].copy() > 0)
        OII1_ew_IVAR = numpy.ma.MaskedArray(hdu['EMLINE_SEW_IVAR'].data[emlc['OIId-3728'],:,:].copy(), mask=hdu[mask_ext].data[emlc['OIId-3728'],:,:].copy() > 0)
        OII1_ew_err=1/numpy.sqrt(OII1_ew_IVAR)
        OII_ew = OII1_ew
        OII_ew_err = OII1_ew_err
        #MASK=fractional(l,Tar,0.2,float('nan'))
        #SPEC=avespec(spx,vormap,MASK,redshift)
        '''if SPEC=='error':
            print('noise invalid!')
            m+=1
            continue'''
        out[0].append(Tar)
        out[1].append(ha_flux)
        out[2].append(ha_flux_err)
        out[3].append(ha_ew)
        out[4].append(ha_ew_err)
        out[5].append(OII_flux)
        out[6].append(OII_flux_err)
        out[7].append(OII_ew)
        out[8].append(OII_ew_err)
        
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
        p_i=np.array([plate_ifu for j in range(0,ha_flux.shape[0])])
        subdat=Table()
        subdat['ha_flux']=ha_flux
        subdat['ha_flux_err']=ha_flux_err
        subdat['ha_flux_mask']=ha_flux_mask
        subdat['ha_ew']=ha_ew
        subdat['ha_ew_err']=ha_ew_err
        subdat['ha_ew_mask']=ha_ew_mask
        subdat['o2_flux']=OII_flux
        subdat['o2_flux_err']=OII_flux_err
        subdat['o2_flux_mask']=OII_flux_mask
        subdat['o2_ew']=OII_ew
        subdat['o2_ew_err']=OII_ew_err
        subdat['o2_ew_mask']=OII_ew_mask
        subdat['plate_ifu']=p_i
        dat.append(subdat)
        i+=1
        
    print(len(dat))
    datall=vstack(dat)
    datall.write('./data_s_ha.fits')
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
