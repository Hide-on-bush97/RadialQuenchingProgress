
from __future__ import print_function
from astropy.io import fits

try:
    import pyfits
except:
    from astropy.io import fits as pyfits
from scipy import ndimage

import glob

from time import clock

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
hdulist1=fits.open('./data_g_bpt.fits',memmap=False)
hdulist=fits.open('./data_g_ha.fits',memmap=False)
data=hdulist[1].data
data1=hdulist1[1].data

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def FIND(name):
    i=0
    while i<6779:
        if drp[1].data[i][2]==name:
            return i
            break
        i+=1
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def fractional(Tar,L,r):
    
    z=drp[1].data[Tar][77]
    objra=drp[1].data[Tar][11]
    objdec=drp[1].data[Tar][12]
    ifura=drp[1].data[Tar][15]
    ifudec=drp[1].data[Tar][16]
    ba=drp[1].data[Tar][84]
    re=drp[1].data[Tar][87]#in arcsecond

    phi=drp[1].data[Tar][85]
    mass=drp[1].data[Tar][83]
    sern=drp[1].data[Tar][94]
    
    l=L
    x=(objra-ifura)*3600*2*np.cos(objdec*np.pi/180)
    y=(objdec-ifudec)*3600*2
    x=-x+l/2.
    y=y+l/2.
    
    
    d=re*r*2# arcsecond -> pixel
    a=d
    b=d*ba  # this is Re remember we need 1.5 Re later
    theta=(phi+90)*3.1416/180
    A=pow(a,2)*pow(math.sin(theta),2)+pow(b,2)*pow(math.cos(theta),2)
    B=2*(pow(b,2)-pow(a,2))*math.cos(theta)*math.sin(theta)
    C=pow(a,2)*pow(math.cos(theta),2)+pow(b,2)*pow(math.sin(theta),2)
    f=-pow(a,2)*pow(b,2)
    
    X=np.array([[i for i in range(0,L)] for x0 in range(0,L)])
    Y=X.transpose()
    #print(np.concatenate([1,2,3,4],[1,2,3,4]))
    AX=[]
    AY=[]
    i=0
    while i<0.95:
        j=0
        while j<0.95:
            AX.append(X+j)
            AY.append(Y+i)
            
            j+=0.1
        i+=0.1
    cord=np.array([AX,AY])
    #cord=cord.transpose()
    out=np.zeros([100,L,L])
    out[(A*pow((cord[0]-x),2)+B*(cord[0]-x)*(cord[1]-y)+C*pow((cord[1]-y),2)<-f)]=1
    out=np.sum(out,axis=0)/100.
    return(out)
def fractionalR(Tar,L):
    
    z=drp[1].data[Tar][77]
    objra=drp[1].data[Tar][11]
    objdec=drp[1].data[Tar][12]
    ifura=drp[1].data[Tar][15]
    ifudec=drp[1].data[Tar][16]
    ba=drp[1].data[Tar][84]
    re=drp[1].data[Tar][87]#in arcsecond

    phi=drp[1].data[Tar][85]
    mass=drp[1].data[Tar][83]
    sern=drp[1].data[Tar][94]
    l=L
    x=(objra-ifura)*3600*2*np.cos(objdec*np.pi/180)
    y=(objdec-ifudec)*3600*2
    x=-x+l/2.
    y=y+l/2.
    
    
    d=re*2# arcsecond -> pixel
    a=d
    b=d*ba  # this is Re remember we need 1.5 Re later
    theta=(phi+90)*3.1416/180
    A=pow(a,2)*pow(math.sin(theta),2)+pow(b,2)*pow(math.cos(theta),2)
    B=2*(pow(b,2)-pow(a,2))*math.cos(theta)*math.sin(theta)
    C=pow(a,2)*pow(math.cos(theta),2)+pow(b,2)*pow(math.sin(theta),2)
    f=-pow(a,2)*pow(b,2)
    
    X=np.array([[i for i in range(0,L)] for x0 in range(0,L)])
    Y=X.transpose()
    #print(np.concatenate([1,2,3,4],[1,2,3,4]))
    
    cord=np.array([X,Y])
    #cord=cord.transpose()
    out=np.zeros([L,L])
    out=pow(-(A*pow((cord[0]-x),2)+B*(cord[0]-x)*(cord[1]-y)+C*pow((cord[1]-y),2))/f,0.5)
    return(out)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def sfdetBPT(o3,o2,ha,s2):
    L=len(ha)
    valid=numpy.zeros((int(L),int(L)))
    sf=numpy.zeros_like(valid)
    i=0
    while(i<L):
        j=0
        while j<L:
            if o3[j][i]>0 and o2[j][i]>0 and ha[j][i]>0 and s2[j][i]>0:
                valid[j][i]=1
                if numpy.log10(o3[j][i]/o2[j][i])<(2/(numpy.log10(s2[j][i]/ha[j][i])-0.65)+2.05) and numpy.log10(s2[j][i]/ha[j][i])<0.65:
                    sf[j][i]=1
            j+=1
        i+=1
    return(valid,sf)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def sfdetMassHa(hasb,mass):
    L=len(hasb)
    valid=numpy.zeros((int(L),int(L)))
    sf=numpy.zeros_like(valid)
    i=0
    while(i<L):
        j=0
        while j<L:
            if hasb[j][i]>0 and mass[j][i]>0:
                valid[j][i]=1
                if hasb[j][i]>mass[j][i]+28.5:
                    sf[j][i]=1
            j+=1
        i+=1
    return(valid,sf)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def percentage(valid,sf,array):
    
    out=numpy.sum(np.multiply(numpy.multiply(sf,array),valid))/float(numpy.sum(numpy.multiply(array,valid)))
    n=numpy.sum(numpy.multiply(array,valid))
    sigma=numpy.sqrt(out*(1-out)/float(n))
    return(out,n,sigma)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def testnring(nring,n):
    L=len(nring)
    i=n
    out=True
    while i<L-1:
        if nring[i]>0:
            out=False
            break
        i+=1
    return(out)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def fixpring(pring,nring,ntheory,rring,sigmaring):
    L=len(pring)
    badvaluemask=[]
    i=0
    while i<L:
        if math.isnan(pring[i])==True or nring[i]/float(ntheory[i])<0.5:
            badvaluemask.append(i)
        i+=1
    pring=numpy.delete(pring,badvaluemask)
    nring=numpy.delete(nring,badvaluemask)
    rring=numpy.delete(rring,badvaluemask)
    ntheory=numpy.delete(ntheory,badvaluemask)
    sigmaring=numpy.delete(sigmaring,badvaluemask)
    return(pring,nring,rring,sigmaring,ntheory)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def lnfit(x,y,yerr):
    
    fitfunc=lambda x,a,b:a*x+b
    popt,pcov=curve_fit(fitfunc,x,y,p0=(0.0,0.0),sigma=yerr) 
    
    perr = np.sqrt(np.diag(pcov))
    return(popt,perr)
                   
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''def qradius(valid,sf,Tar,p,jd):
    r1=0.20
    p1=1
    p2=1
    k1=1
    k2=1
    sig=0
    flag1=0
    flag2=0
    flag3=0
    out1=2.5
    out2=2.5
    out3=2.5
    ba=drp[1].data[Tar][84]
    re=drp[1].data[Tar][87]#in arcsecond
    re=re*2# in spaxel
    while r1<2.5:
        #array1,array2=fractional(len(sf),Tar,r1,jd)
        array1=fractional(len(sf),Tar,r1,jd)
        p2=p1
        
        sig2=sig
        p1,n,sig=percentage(valid,sf,array1)
        
        if p1-sig<p and flag1==0:
            out1=r1
            flag1=1
        if p1<p and flag2==0:
            out2=r1
            flag2=1
        if p1+sig<p and flag3==0:
            out3=r1
            flag3=1
        if p1+sig>p and p1+sig>p2+sig2:
            flag1=0
            flag2=0
            flag3=0
        if flag1==1 and flag2==1 and flag3==1:
            break
        r1+=0.01
        
    return(out1,out2,out3)'''
def qradiusAGN(valid,sf,Tar,p):
    re=drp[1].data[Tar][87]
    ba=drp[1].data[Tar][84]
    L=len(valid)
    rring1=np.array([0,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2]) # 1-49 pixels in radius, equals to 1-49 arcsconds in Diameter
    
    rring=[]# r/re in radius, equals to 1-49 pixels in radius
    pring=[]
    nring=[]
    sigmaring=[]
    ntheory=[]
    i=0
    while i<len(rring1)-1:
        rring.append(rring1[i+1])
        mask=fractional(Tar,L,rring1[i+1])-fractional(Tar,L,rring1[i])
        a,b,c=percentage(valid,sf,mask)
        pring.append(a)
        nring.append(b)
        sigmaring.append(c)
        ntheory.append(np.pi*pow(rring1[i+1],2)*ba*pow(2*re,2)-np.pi*pow(rring1[i],2)*ba*pow(2*re,2))
        i+=1
    pring,nring,rring,sigmaring,ntheory=fixpring(pring,nring,ntheory,rring,sigmaring)
    index1=np.where(pring>0.2)
    index2=np.where(pring>0.05)
    out='N'
    if len(index1[0])>0:
        #print(pring)
        out='Y'
    if len(index2[0])>0 and len(index1[0])==0:
        #print(pring)
        out='W'
    if len(index2[0])==0 and len(index1[0])>0:
        #print(pring)
        out='err'
    return(Tar,out,pring,rring)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def sfdetHaO2(haew,o2ew,hae,o2e):
    L=len(haew)
    valid1=numpy.zeros((int(L),int(L)))
    valid2=numpy.zeros((int(L),int(L)))
    '''valid3=numpy.zeros((int(L),int(L)))
    valid4=numpy.zeros((int(L),int(L)))
    valid5=numpy.zeros((int(L),int(L)))
    valid6=numpy.zeros((int(L),int(L)))'''
    
    i=0
    while(i<L):
        j=0
        while j<L:
            if math.isnan(haew[j][i])==False:
                valid1[j][i]=1
            '''if math.isnan(o2ew[j][i])==False:
                valid2[j][i]=1
            if math.isnan(o2ew[j][i])==False:# and math.isnan(haew[j][i])==False:
                valid1[j][i]=1'''
            if math.isnan(haew[j][i])==False and (haew[j][i]> 5*hae[j][i] or hae[j][i] <= 0.2):
                valid2[j][i]=1
            #if math.isnan(o2ew[j][i])==False and (o2ew[j][i]> o2e[j][i] or o2e[j][i] <= 1):
            #    valid2[j][i]=1
            #if math.isnan(o2ew[j][i])==False and math.isnan(haew[j][i])==False and (haew[j][i]> 5*hae[j][i] or hae[j][i] <= 0.2) and (o2ew[j][i]> o2e[j][i] or o2e[j][i] <= 1):
            #    valid6[j][i]=1
            j+=1
        i+=1
    return(valid1,valid2)
#valid1,valid2,
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


def TEST(spec):
    L=len(spec['flux'])
    k=0
    i=0
    while i<L:
        if (spec['noise'][i]<0 or spec['noise'][i]==float('inf') or spec['noise'][i]==float('nan')  ):
            k=1
            print(k,spec['noise'][i])
            break
        i+=1
    return(k)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def avespec(spx,vormap,mask,z):
    l=len(mask)
    k=0.0
    spec=spectrum(int(l/2),int(l/2),z,spx,vormap)
    wave=spec['wavelength']
    spec=np.zeros_like(spec['flux'])
    noise=np.zeros_like(spec)
    
    i=0
    while i<l:
        j=0
        while(j<l):
            if mask[i][j]>0:
                aaa=spectrum(i,j,z,spx,vormap)
                if TEST(aaa)==1:
                    j+=1
                    continue
                spec+=mask[i][j]*aaa['flux']
                noise+=mask[i][j]*aaa['noise']
                k+=1
            j+=1
        i+=1
    name=['flux','wavelength','noise']
    if k!=0:
        spec=spec/k
        noise=noise/k
        out=dict(zip(name,[spec,wave,noise]))
    if k==0:
        out='error'
    return out
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
    out=[]
    
    
    
    m=0
    k=1
    print(L)
    while m<L:
        Tar=int(table.col_values(3)[m])
        #agn,tc=ANGtypechange(Tar)
        re=drp[1].data[Tar][87]#in arcsecond
 #half long axis in spaxels
        if table.col_values(2)[m]!='P' or re<=0:
            #print(table.col_values(2)[m])
            m+=1
            continue
        
        out.append(Tar)
        
        m+=1
    return(np.transpose(numpy.array(out)))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def prob(x0,y0,sigmax,sigmay,k,b):
    return(1-norm.cdf((k*x0-y0+b)/np.sqrt(sigmay*sigmay+k*k*sigmax*sigmax)))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def task(Tar):
    
    #print(Tar,type(Tar)) 
    index=np.where(data['plate_ifu']==drp[1].data[Tar][2])
    l=int(np.sqrt(len(index[0])))
    ha_flux=np.reshape(data['ha_flux'][index],newshape=[l,l],order='C')
    ha_flux_mask=np.reshape(data['ha_flux_mask'][index],newshape=[l,l],order='C')
    ha_flux_ivar=np.reshape(data['ha_flux_ivar'][index],newshape=[l,l],order='C')

    ha_ew=np.reshape(data['ha_ew'][index],newshape=[l,l],order='C')
    ha_ew_mask=np.reshape(data['ha_ew_mask'][index],newshape=[l,l],order='C')
    ha_ew_ivar=np.reshape(data['ha_ew_ivar'][index],newshape=[l,l],order='C')

    o2_3727_flux=np.reshape(data['o2_3727_flux'][index],newshape=[l,l],order='C')
    o2_3727_flux_mask=np.reshape(data['o2_3727_flux_mask'][index],newshape=[l,l],order='C')
    o2_3727_flux_ivar=np.reshape(data['o2_3727_flux_ivar'][index],newshape=[l,l],order='C')

    o2_3729_flux=np.reshape(data['o2_3729_flux'][index],newshape=[l,l],order='C')
    o2_3729_flux_mask=np.reshape(data['o2_3729_flux_mask'][index],newshape=[l,l],order='C')
    o2_3729_flux_ivar=np.reshape(data['o2_3729_flux_ivar'][index],newshape=[l,l],order='C')
    
    o2_3727_ew=np.reshape(data['o2_3727_ew'][index],newshape=[l,l],order='C')
    o2_3727_ew_mask=np.reshape(data['o2_3727_ew_mask'][index],newshape=[l,l],order='C')
    o2_3727_ew_ivar=np.reshape(data['o2_3727_ew_ivar'][index],newshape=[l,l],order='C')

    o2_3729_ew=np.reshape(data['o2_3729_ew'][index],newshape=[l,l],order='C')
    o2_3729_ew_mask=np.reshape(data['o2_3729_ew_mask'][index],newshape=[l,l],order='C')
    o2_3729_ew_ivar=np.reshape(data['o2_3729_ew_ivar'][index],newshape=[l,l],order='C')
            
    mask=np.ones_like(ha_flux)
    mask[(ha_ew_mask & 1073741824==0)&(ha_ew_ivar != 0)
         &(o2_3727_ew_mask & 1073741824==0)&(o2_3727_ew_ivar != 0)
         &(o2_3729_ew_mask & 1073741824==0)&(o2_3729_ew_ivar != 0)]=0
    inverse_mask=np.ones_like(mask)-mask
    
    ha_flux=numpy.ma.MaskedArray(ha_flux,ha_flux_mask>0)
    ha_flux_err=1/np.sqrt(numpy.ma.MaskedArray(ha_flux_ivar,ha_flux_mask>0))

    ha_ew=numpy.ma.MaskedArray(ha_ew,mask>0)
    ha_ew_err=1/np.sqrt(numpy.ma.MaskedArray(ha_ew_ivar,mask>0))

    o2_ew=numpy.ma.MaskedArray(o2_3727_ew+o2_3729_ew,mask>0)
    o2_ew_err=1/np.sqrt(numpy.ma.MaskedArray(o2_3727_ew_ivar,mask>0))+1/np.sqrt(numpy.ma.MaskedArray(o2_3729_ew_ivar,mask>0))
    
    hb_flux=numpy.ma.MaskedArray(np.reshape(data1['hb_flux'][index],newshape=[l,l],order='C'),
                                         mask=np.reshape(data1['hb_flux_mask'][index],newshape=[l,l],order='C'))
            
    hb_flux_err=numpy.ma.MaskedArray(np.reshape(data1['hb_flux_err'][index],newshape=[l,l],order='C'),
                                         mask=np.reshape(data1['hb_flux_mask'][index],newshape=[l,l],order='C'))
    o3_flux=numpy.ma.MaskedArray(np.reshape(data1['o3_flux'][index],newshape=[l,l],order='C'),
                                         mask=np.reshape(data1['o3_flux_mask'][index],newshape=[l,l],order='C'))
            
    o3_flux_err=numpy.ma.MaskedArray(np.reshape(data1['o3_flux_err'][index],newshape=[l,l],order='C'),
                                         mask=np.reshape(data1['o3_flux_mask'][index],newshape=[l,l],order='C'))
    s2_flux=numpy.ma.MaskedArray(np.reshape(data1['s2_flux'][index],newshape=[l,l],order='C'),
                                         mask=np.reshape(data1['s2_flux_mask'][index],newshape=[l,l],order='C'))
            
    s2_flux_err=numpy.ma.MaskedArray(np.reshape(data1['s2_flux_err'][index],newshape=[l,l],order='C'),
                                         mask=np.reshape(data1['s2_flux_mask'][index],newshape=[l,l],order='C'))
    n2_flux=numpy.ma.MaskedArray(np.reshape(data1['n2_flux'][index],newshape=[l,l],order='C'),
                                         mask=np.reshape(data1['n2_flux_mask'][index],newshape=[l,l],order='C'))
            
    n2_flux_err=numpy.ma.MaskedArray(np.reshape(data1['n2_flux_err'][index],newshape=[l,l],order='C'),
                                         mask=np.reshape(data1['n2_flux_mask'][index],newshape=[l,l],order='C'))
            
    R=fractionalR(Tar,len(ha_flux))
    
    #metal,age=ppxf_population_example_sdss(SPEC,redshift)
    r3=np.log10(o3_flux/hb_flux)
    s2=np.log10(s2_flux/ha_flux)
    n2=np.log10(n2_flux/ha_flux)
    x1=0.63*n2+0.51*s2+0.59*r3
    y1=-0.63*n2+0.78*s2
    z1=-0.46*n2-0.37*s2+0.81*r3
    #__________________________________________________________
    #begin matrices for all
    index_5lines=np.where((ha_ew>3*ha_ew_err)
                          &(inverse_mask==1)    
                          &(ha_flux>3*ha_flux_err)
                          &(s2_flux>3*s2_flux_err)&(o3_flux>3*o3_flux_err)
                          &(hb_flux>3*hb_flux_err)&(n2_flux>3*n2_flux_err)
                          &(ha_flux.mask==0)
                          &(hb_flux.mask==0)&(o3_flux.mask==0)
                          &(s2_flux.mask==0)&(n2_flux.mask==0))
    varray_5lines=np.zeros_like(ha_ew.data)
    varray_5lines[index_5lines]=1
    
    index_2lines=np.where((ha_ew>3*ha_ew_err)
                          &(inverse_mask==1)    
                          &(ha_flux>3*ha_flux_err)
                          &(n2_flux>3*n2_flux_err)
                          &(ha_flux.mask==0)
                          &(n2_flux.mask==0))
    varray_2lines=np.zeros_like(ha_ew.data)
    varray_2lines[index_2lines]=1
    
    index_halines=np.where((ha_ew>3*ha_ew_err)
                           &(inverse_mask==1)    
                           &(ha_flux>3*ha_flux_err)
                           &(ha_flux.mask==0))
    varray_halines=np.zeros_like(ha_ew.data)
    varray_halines[index_halines]=1

    index_inverse_mask=np.where((inverse_mask==1))
    varray_inverse_mask=np.zeros_like(ha_ew.data)
    varray_inverse_mask[index_inverse_mask]=1
    '''plt.figure(figsize=(12,12))
    plt.subplot(221)
    plt.imshow(varray_5lines)
    plt.subplot(222)
    plt.imshow(varray_2lines)
    plt.subplot(223)
    plt.imshow(varray_halines)
    plt.subplot(224)
    plt.imshow(varray_inverse_mask)
    plt.show()'''
    #end matrices for all
    #__________________________________________________________
    #begin matrices for re<0.5
    index_5lines_5=np.where((ha_ew>3*ha_ew_err)
                            &(R<0.5)
                            &(inverse_mask==1)    
                            &(ha_flux>3*ha_flux_err)
                            &(s2_flux>3*s2_flux_err)&(o3_flux>3*o3_flux_err)
                            &(hb_flux>3*hb_flux_err)&(n2_flux>3*n2_flux_err)
                            &(ha_flux.mask==0)
                            &(hb_flux.mask==0)&(o3_flux.mask==0)
                            &(s2_flux.mask==0)&(n2_flux.mask==0))
    varray_5lines_5=np.zeros_like(ha_ew.data)
    varray_5lines_5[index_5lines_5]=1
    
    index_2lines_5=np.where((ha_ew>3*ha_ew_err)
                            &(R<0.5)
                            &(inverse_mask==1)    
                            &(ha_flux>3*ha_flux_err)
                            &(n2_flux>3*n2_flux_err)
                            &(ha_flux.mask==0)
                            &(n2_flux.mask==0))
    varray_2lines_5=np.zeros_like(ha_ew.data)
    varray_2lines_5[index_2lines_5]=1
    
    index_halines_5=np.where((ha_ew>3*ha_ew_err)
                             &(R<0.5)
                             &(inverse_mask==1)    
                             &(ha_flux>3*ha_flux_err)
                             &(ha_flux.mask==0))
    varray_halines_5=np.zeros_like(ha_ew.data)
    varray_halines_5[index_halines_5]=1

    index_inverse_mask_5=np.where((R<0.5)
                                 &(inverse_mask==1))
    varray_inverse_mask_5=np.zeros_like(ha_ew.data)
    varray_inverse_mask_5[index_inverse_mask_5]=1
    #end matrices for re<0.5
    #upper limit of n2 for region with valid ha and invalid n2
    #print(len(index_n2_upper[0]))
    index_n2_upper=np.where((varray_halines==1)
                            &(varray_2lines==0))
    
    n2_upper=3*n2_flux_err[index_n2_upper]
    ha_upper=ha_flux[index_n2_upper]
    ha_ew_upper=ha_ew[index_n2_upper]

    index_n2_upper=np.where((varray_halines_5==1)
                            &(varray_2lines_5==0))
    
    n2_upper_5=3*n2_flux_err[index_n2_upper]
    ha_upper_5=ha_flux[index_n2_upper]
    ha_ew_upper_5=ha_ew[index_n2_upper]
    #print(len(index_n2_upper[0]))
    #________________________________________________
    varray10=inverse_mask[R<1]
    varray5=inverse_mask[R<0.5]
    v=np.ones_like(ha_ew.data)
    v10=v[R<1]
    v5=v[R<0.5]
    #________________________________________________
    #begin agn selection with whan as suupliment, ALL
    agn=np.zeros_like(ha_ew.data)
    index=np.where((ha_ew>3)&(x1>-1.51*y1*y1-0.355*y1+0.002)
                   &(varray_5lines==1))
    
    index_sub=np.where((ha_ew>3)&(n2>-0.4)
                       &(varray_5lines==0)&(varray_2lines==1))
    
    #test begin
    test1=np.zeros_like(ha_ew.data)
    test2=np.zeros_like(ha_ew.data)
    test1[index]=1
    test2[index_sub]=1
    test_index=np.where((test1==1)&(test2==1))
    '''plt.figure(figsize=(12,12))
    plt.subplot(221)
    plt.imshow(varray_5lines)
    plt.subplot(222)
    plt.imshow(varray_2lines)
    plt.subplot(223)
    plt.imshow(test1)
    plt.subplot(224)
    plt.imshow(test2)
    plt.show()'''
    if len(test_index[0])!=0:
        
        print(len(test_index[0]),'err')
    #test end
    agn[index]=1
    agn[index_sub]=1
    out=qradiusAGN(inverse_mask,agn,Tar,0.1)
    # end agn selection with whan as suupliment, ALL
    #_________________________________________________________________
    #begin agn selection with whan as suupliment, re<0.5
    agn_5=np.zeros_like(ha_ew.data)
    index=np.where((ha_ew>3)&(x1>-1.51*y1*y1-0.355*y1+0.002)
                   &(varray_5lines_5==1))
    
    index_sub=np.where((ha_ew>3)&(n2>-0.4)
                       &(varray_5lines_5==0)&(varray_2lines_5==1))
    
    agn_5[index]=1
    agn_5[index_sub]=1
    out_5=qradiusAGN(inverse_mask,agn_5,Tar,0.1)
    # end agn selection with whan as suupliment, ALL
    #metal,age=ppxf_population_example_sdss(SPEC,redshift)
    
    #print(len(index1[0])/float(len(index2[0])),out)
    
    return(out_5,out,
           np.sum(varray5),np.sum(varray10),np.sum(v5),np.sum(v10),
           np.sum(varray_5lines),np.sum(varray_2lines),np.sum(varray_halines),np.sum(varray_inverse_mask),
           np.sum(varray_5lines_5),np.sum(varray_2lines_5),np.sum(varray_halines_5),np.sum(varray_inverse_mask_5),
           n2_upper,ha_upper,ha_ew_upper,
           n2_upper_5,ha_upper_5,ha_ew_upper_5)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def output(out):
    
    L=len(out)
    i=0
    v5=0
    v10=0
    va5=0
    va10=0
    lines5=0
    lines2=0
    haline=0
    inverse_mask=0
    lines5_5=0
    lines2_5=0
    haline_5=0
    inverse_mask_5=0
    n2_u=np.array([])
    ha_u=np.array([])
    ha_ew_u=np.array([])
    n2_u_5=np.array([])
    ha_u_5=np.array([])
    ha_ew_u_5=np.array([])
    book=open('AGN0613ha3R5.txt','w')
    plt.figure(figsize=(20,20))
    while i<L:
        v5+=out[i][2]
        v10+=out[i][3]
        va5+=out[i][4]
        va10+=out[i][5]
        lines5+=out[i][6]
        lines2+=out[i][7]
        haline+=out[i][8]
        inverse_mask+=out[i][9]
        lines5_5+=out[i][10]
        lines2_5+=out[i][11]
        haline_5+=out[i][12]
        inverse_mask_5+=out[i][13]
        n2=out[i][14]
        ha=out[i][15]
        ha_ew=out[i][16]
        n2_5=out[i][17]
        ha_5=out[i][18]
        ha_ew_5=out[i][19]
        #print(len(n2),len(n2_5),inverse_mask)
        
        n2_u=np.concatenate((n2_u,n2))
        ha_u=np.concatenate((ha_u,ha))
        ha_ew_u=np.concatenate((ha_ew_u,ha_ew))
        n2_u_5=np.concatenate((n2_u_5,n2_5))
        ha_u_5=np.concatenate((ha_u_5,ha_5))
        ha_ew_u_5=np.concatenate((ha_ew_u_5,ha_ew_5))
        book.write(str(int(out[i][0][0]))+'*'+str(out[i][0][1])+'\n')
        if out[i][0][1]=='N' and (out[i][1][1]=='Y' or out[i][1][1]=='W'):
            print(out[i][0][1],out[i][1][1])
            plt.plot(out[i][1][3],out[i][1][2])
        i+=1
    print('\n0.5re valid array coverage:', round(v5/float(va5),4),
          '\n1re valid array coverage:', round(v10/float(va10),4),
          '\n5lines valid fraction', round(lines5/float(inverse_mask),4),
          '\n2lines valid fraction', round(lines2/float(inverse_mask),4),
          '\nha-line valid fraction', round(haline/float(inverse_mask),4),
          '\n5lines valid fraction, <0.5re', round(lines5_5/float(inverse_mask_5),4),
          '\n2lines valid fraction, <0.5re', round(lines2_5/float(inverse_mask_5),4),
          '\nha-line valid fraction, <0.5re', round(haline_5/float(inverse_mask_5),4),
          '\n')
    plt.xlim(0,1.5)
    plt.ylim(0,1)
    plt.xlabel('r/Re',fontsize=35)
    plt.ylabel('AGN percentage',fontsize=35)
    plt.tick_params(labelsize=35)
    plt.tight_layout()
    plt.savefig('0613.png',dpi=400)
    plt.close()
    book.close()
    plt.figure(figsize=(20,20))
    plt.scatter(np.log10(n2_u/ha_u),np.log10(ha_ew_u),c='blue',s=1)
    plt.xlim(-1,0.6)
    plt.ylim(-1,2.3)
    plt.xlabel('log10(n2/ha)',fontsize=35)
    plt.ylabel('log10(Ha_ew)',fontsize=35)
    xxxx=np.array([i/10. for i in range(-10,6)])
    yyyy=xxxx*0+0.477
    yyyy2=xxxx*0+0.778
    yyyy3=np.array([i/10. for i in range(-10,25)])
    xxxx2=yyyy3*0-0.4
    xxxx3=yyyy3*0
    plt.plot(xxxx,yyyy,'r')
    plt.plot(xxxx,yyyy2,'r')
    plt.plot(xxxx2,yyyy3,'r')
    plt.tick_params(labelsize=35)
    plt.tight_layout()
    plt.savefig('0613_2.png',dpi=400)
    plt.close()
    plt.figure(figsize=(20,20))
    plt.scatter(np.log10(n2_u_5/ha_u_5),np.log10(ha_ew_u_5),c='blue',s=1)
    plt.title('re<0.5',fontsize=30)
    plt.xlim(-1,0.6)
    plt.ylim(-1,2.3)
    plt.xlabel('log10(n2/ha)',fontsize=35)
    plt.ylabel('log10(Ha_ew)',fontsize=35)
    xxxx=np.array([i/10. for i in range(-10,6)])
    yyyy=xxxx*0+0.477
    yyyy2=xxxx*0+0.778
    yyyy3=np.array([i/10. for i in range(-10,25)])
    xxxx2=yyyy3*0-0.4
    xxxx3=yyyy3*0
    plt.plot(xxxx,yyyy,'r')
    plt.plot(xxxx,yyyy2,'r')
    plt.plot(xxxx2,yyyy3,'r')
    plt.tick_params(labelsize=35)
    plt.tight_layout()
    plt.savefig('0613_3.png',dpi=400)
    plt.close()
    return(0)
    
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def run(data):
    L=len(data)
    out=[]
    i=0
    while i<L:
        out.append(task(data[i]))
        i+=1
    return(out)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
flag1=time.time()
#book=open('0711p1_3.txt','w')
Tar=feeddata()
flag4=time.time()
print(len(Tar))
cores = mp.cpu_count()


pool = mp.Pool(processes=cores)

out1=pool.map(task,Tar)

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
print('process data    ',m,'min,',s,'s')

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
output(out1)
flag3=time.time()
Time=flag3-flag2
m=int(Time/60.)
s=Time-m*60
print('output data    ',m,'min,',s,'s')

Time=flag3-flag1
m=int(Time/60.)
s=Time-m*60
#out1=run(data)
print('Total runn time' ,m,'min,',s,'s')
print('\n data size='+str(len(out1))+'\n')
