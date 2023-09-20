
from __future__ import print_function
from astropy.io import fits
from scipy.stats import norm
from matplotlib.patches import Ellipse
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

from matplotlib import pyplot as plt
from astropy.table import vstack, Table
import numpy
from scipy import integrate
import math
import numpy as np

import time
import astropy.units as u
import speclite.filters
from astropy.cosmology import WMAP9 as cosmo

a=xlrd.open_workbook('/pegasus/czh276/age_re/Tarlist0615.xls')# get the target names from excel
drp=fits.open('/pegasus/czh276/age_re/drpall-v3_1_1.fits') #drpall file
table=a.sheets()[0]
hdulist1=fits.open('./data_g_bpt.fits',memmap=False)
hdulist=fits.open('./data_g_ha.fits',memmap=False)
data=hdulist[1].data
data1=hdulist1[1].data
dat=fits.open('./bgdata.fits')
x_bg=dat[1].data['Ha']
y_bg=dat[1].data['OII']


index=np.where(5*x_bg-7>y_bg)
x_bg_q=x_bg[index]
y_bg_q=y_bg[index]
print(x_bg.shape,x_bg_q.shape)

    
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
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def percentage(valid,sf,array,mask2):
    
    out=np.sum(sf*array*valid*mask2)/float(np.sum(array*valid*mask2))
    n=np.sum(array*valid*mask2)
    sigma=np.sqrt(out*(1-out)/float(n))
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
#    print(L)
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
def fixpringAGN(pring,nring,ntheory,rring,sigmaring):
    L=len(pring)
#    print(L)
    badvaluemask=[]
    i=0
    while i<L:
        if math.isnan(pring[i])==True:
            badvaluemask.append(i)
        i+=1
    pring=numpy.delete(pring,badvaluemask)
    nring=numpy.delete(nring,badvaluemask)
    rring=numpy.delete(rring,badvaluemask)
    ntheory=numpy.delete(ntheory,badvaluemask)
    sigmaring=numpy.delete(sigmaring,badvaluemask)
    return(pring,nring,rring,sigmaring,ntheory)

                   
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def qradiusR(valid,sf,Tar,p):
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
        rring.append(np.sqrt((pow(rring1[i+1],2)+pow(rring1[i],2))/2.))
        mask=fractional(Tar,L,rring1[i+1])-fractional(Tar,L,rring1[i])
        a,b,c=percentage(valid,sf,mask)
        pring.append(a)
        nring.append(b)
        sigmaring.append(c)
        ntheory.append(np.pi*pow(rring1[i+1],2)*ba*pow(2*re,2)-np.pi*pow(rring1[i],2)*ba*pow(2*re,2))
        i+=1
    pring,nring,rring,sigmaring,ntheory=fixpring(pring,nring,ntheory,rring,sigmaring)
    
    
    return(Tar,pring,nring,rring,sigmaring,ntheory)
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
        rring.append(np.sqrt((pow(rring1[i+1],2)+pow(rring1[i],2))/2.))
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
def qradius(valid,sf,Tar,p,mask2):
    
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
        rring.append(np.sqrt((pow(rring1[i+1],2)+pow(rring1[i],2))/2.))
        mask=fractional(Tar,L,rring1[i+1])-fractional(Tar,L,rring1[i])
        a,b,c=percentage(valid,sf,mask,mask2)
        pring.append(a)
        nring.append(b)
        sigmaring.append(c)
        ntheory.append(np.pi*pow(rring1[i+1],2)*ba*pow(2*re,2)-np.pi*pow(rring1[i],2)*ba*pow(2*re,2))
        i+=1
    pring,nring,rring,sigmaring,ntheory=fixpring(pring,nring,ntheory,rring,sigmaring)
    pmid=pring
#    print(L)
    out1=2.5
    out2=2.5
    L=len(pring)
    print(pring,len(pring))
#    print(pring)
#    print(L)
    if L>=2:
        i=0
        while i<L:
            if i!=L-1:
                if pmid[i]<0.5 and pmid[i+1]>0.5:
                    out1=(0.5-pmid[i])*(rring[i]-rring[i+1])/(pmid[i]-pmid[i+1])+rring[i]
                    break
            i+=1
        i=0
        while i<L:
            if i!=L-1:
                if pmid[i]>0.5 and pmid[i+1]<0.5:
                    out2=(0.5-pmid[i])*(rring[i]-rring[i+1])/(pmid[i]-pmid[i+1])+rring[i]
                    break
            i+=1

    if L==0 or L==1:
        out1='invalid'
        out2='invalid'
#        print('test')
    if out1==2.5 and out2==2.5:
        out1=pmid[0]
        out2='invalid'
    if L==0:
        out3='invalid'
        out4='invalid'
        out5='invalid'
        out6='invalid'
    if L!=0:
        out3=pring[0]
        out4=rring[0]
        out5=pring[L-1]
        out6=rring[L-1]
#        print('test')
    return(Tar,out1,out2,out3,out4,out5,out6)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
    out=[[],[]]
    
    
    
    m=0
    k=1
    print(L)
    while m<L:
        Tar=int(table.col_values(0)[m])
        AGN=table.col_values(1)[m]
        re=drp[1].data[Tar][87]#in arcsecond
 #half long axis in spaxels
        if re<=0:
            
            print(table.col_values(2)[m],'Impossible',re)
            
            m+=1
            continue
        
        out[0].append(Tar)
        out[1].append(AGN)
        m+=1
    return(np.transpose(numpy.array(out)))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def prob1(x0,y0,dx,dy):

    
    #print(np.amax(x))
    PEH=np.sum(np.exp(-((x_bg_q-x0)**2)/(2*dx**2)-((y_bg_q-y0)**2)/(2*dy**2)),axis=0)
    tot=np.sum(np.exp(-((x_bg-x0)**2)/(2*dx**2)-((y_bg-y0)**2)/(2*dy**2)),axis=0)
    if tot<0.00005:
        if y0>5*x0-7:
            out=0
        else:
            out=1
    else:
        out=PEH/tot
    #print(PEH,tot,x0,y0,out)
    return(out)
def prob(x0,y0,sigmax,sigmay,k,b):
    x=x0.flatten(order='C')
    y=y0.flatten(order='C')
    dx=sigmax.flatten(order='C')
    dy=sigmay.flatten(order='C')
    l=len(x0)
    #print(int(l/2-3))
    res=[]
    for i in range(len(x)):
        res.append(prob1(x[i],y[i],dx[i],dy[i]))
    res=np.array(res)
    out=np.reshape(res,newshape=[int(np.sqrt(len(x))),int(np.sqrt(len(x)))],order='C')
    #print(x0.shape,x0[int(l/2-3):int(l/2+3),int(l/2-3):int(l/2+3)])
    #print(y0.shape,y0[int(l/2-3):int(l/2+3),int(l/2-3):int(l/2+3)])
    #print(out.shape,out[int(l/2-3):int(l/2+3),int(l/2-3):int(l/2+3)])
    return(out)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def task(datainput):
    Tar=int(datainput[0])
    AGN=datainput[1]
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
    
    mask2=np.ones_like(inverse_mask)
    sf=prob(ha_ew,o2_ew,ha_ew_err,o2_ew_err,5,-7)
    sf=mx = np.ma.masked_invalid(sf)
    #print('SF',sf.shape,sf[int(l/2-3):int(l/2+3),int(l/2-3):int(l/2+3)])
    if AGN=='Y':
        #begin matrices for re<0.5
        index_5lines_5=np.where((ha_ew>3*ha_ew_err)
                                #&(R<0.5)
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
                                #&(R<0.5)
                                &(inverse_mask==1)    
                                &(ha_flux>3*ha_flux_err)
                                &(n2_flux>3*n2_flux_err)
                                &(ha_flux.mask==0)
                                &(n2_flux.mask==0))
        varray_2lines_5=np.zeros_like(ha_ew.data)
        varray_2lines_5[index_2lines_5]=1
    
        index_halines_5=np.where((ha_ew>3*ha_ew_err)
                                 #&(R<0.5)
                                 &(inverse_mask==1)    
                                 &(ha_flux>3*ha_flux_err)
                                 &(ha_flux.mask==0))
        varray_halines_5=np.zeros_like(ha_ew.data)
        varray_halines_5[index_halines_5]=1

        index_inverse_mask_5=np.where((inverse_mask==1))
                                      #&(R<0.5))
        varray_inverse_mask_5=np.zeros_like(ha_ew.data)
        varray_inverse_mask_5[index_inverse_mask_5]=1
    #end matrices for re<0.5
        
        index=np.where((ha_ew>3)&(x1>-1.51*y1*y1-0.355*y1+0.002)
                   &(varray_5lines_5==1))
    
        index_sub=np.where((ha_ew>3)&(n2>-0.4)
                       &(varray_5lines_5==0)&(varray_2lines_5==1))
    
        
        sf[index]=1
        sf[index_sub]=1
        print('treat agn as Q\n')
    if AGN=='W':
        #begin matrices for re<0.5
        index_5lines_5=np.where((ha_ew>3*ha_ew_err)
                                #&(R<0.5)
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
                                #&(R<0.5)
                                &(inverse_mask==1)    
                                &(ha_flux>3*ha_flux_err)
                                &(n2_flux>3*n2_flux_err)
                                &(ha_flux.mask==0)
                                &(n2_flux.mask==0))
        varray_2lines_5=np.zeros_like(ha_ew.data)
        varray_2lines_5[index_2lines_5]=1
    
        index_halines_5=np.where((ha_ew>3*ha_ew_err)
                                 #&(R<0.5)
                                 &(inverse_mask==1)    
                                 &(ha_flux>3*ha_flux_err)
                                 &(ha_flux.mask==0))
        varray_halines_5=np.zeros_like(ha_ew.data)
        varray_halines_5[index_halines_5]=1

        index_inverse_mask_5=np.where((inverse_mask==1))
                                      #&(R<0.5))
        varray_inverse_mask_5=np.zeros_like(ha_ew.data)
        varray_inverse_mask_5[index_inverse_mask_5]=1
    #end matrices for re<0.5
        
        index=np.where((ha_ew>3)&(x1>-1.51*y1*y1-0.355*y1+0.002)
                   &(varray_5lines_5==1))
    
        index_sub=np.where((ha_ew>3)&(n2>-0.4)
                       &(varray_5lines_5==0)&(varray_2lines_5==1))
    
        
        mask2[index]=0
        mask2[index_sub]=0
        print('treat agn as null\n')
    
    out=qradius(inverse_mask,sf,Tar,0.50,mask2)

    #book.write(str(out[0])+'*'+str(out[1])+'*'+str(out[2])+'*'+str(out[3])+'*'+str(out[4])+'\n')
    '''if r2 > 0.2 and r2 < 1.5:
        out=[Tar,r1,r2,r3,l]'''
    return(out)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def image1(Tar,ax,L):
    plate=drp[1].data[Tar][0]
    ifu=drp[1].data[Tar][1]
    objra=drp[1].data[Tar][11]
    objdec=drp[1].data[Tar][12]
    ifura=drp[1].data[Tar][15]
    ifudec=drp[1].data[Tar][16]
    ba=drp[1].data[Tar][84]
    re=drp[1].data[Tar][87]#in arcsecond
    
    dra=objra-ifura
    ddec=objdec-ifudec
    
    phi=drp[1].data[Tar][85]
    
    #im=plt.imread('/pegasus/czh276/age_re/images/'+str(plate)+'/images/'+str(ifu)+'.png')
    #print(type(im),len(im))
    #print(phi,plate,ifu)
    print(objra, objdec,re)
    circle1= Ellipse((L/2.+dra*3600*2, L/2.+ddec*3600*2),width=4*re*ba, height=4*re,angle=phi, color='red', fill=False)
    
    circle2= Ellipse((L/2.+dra*3600*2, L/2.+ddec*3600*2),width=10*re*ba, height=10*re,angle=phi, color='red', fill=False)
    
    
    ax.add_artist(circle1)
    ax.add_artist(circle2)

def image2(Tar,qr,age,ax,color):
    plate=drp[1].data[Tar][0]
    ifu=drp[1].data[Tar][1]
    objra=drp[1].data[Tar][11]
    objdec=drp[1].data[Tar][12]
    ifura=drp[1].data[Tar][15]
    ifudec=drp[1].data[Tar][16]
    ba=drp[1].data[Tar][84]
    re=drp[1].data[Tar][87]#in arcsecond
    
    dra=objra-ifura
    ddec=objdec-ifudec
    
    phi=drp[1].data[Tar][85]
    if os.path.exists('/pegasus/czh276/age_re/images/'+str(plate)+'/images/'+str(ifu)+'.png')==False:
        im=np.zeros((100,100))
    else:
        
        im=plt.imread('/pegasus/czh276/age_re/images/'+str(plate)+'/images/'+str(ifu)+'.png')
    #print(type(im),len(im))
    #print(phi,plate,ifu)
    L=len(im)
    circle1= Ellipse((L/2.-dra*3600*L/50.0, L/2.-ddec*3600*L/50.0),width=L/25.0*re*ba, height=L/25.0*re,angle=-phi, color='red', fill=False)
    
    circle3= Ellipse((L/2.-dra*3600*L/50.0, L/2.-ddec*3600*L/50.0),width=L/25.0*re*ba*2.5, height=L/25.0*re*2.5,angle=-phi, color=color, fill=False)
    
    ax.imshow(im)
    ax.add_artist(circle1)
    
    ax.add_artist(circle3)
    
    #ax.text(500,100,str(round(age,2)),color='white')
    

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def outputR(out):
    print(out[0][0][3],out[0][0][5],out[0][0][2])
    i=0
    while i<len(out):
        print('test1')
        plt.figure(figsize=(18,12))
        Tar=out[i][0][0]
        pring=out[i][0][1]
        nring=out[i][0][2]
        rring=out[i][0][3]
        sigmaring=out[i][0][4]
        ntheory=out[i][0][5]
        print('test3')
        plt.subplot(231)
        plt.errorbar(rring,pring,yerr=sigmaring)
        plt.xlabel('R/Re',fontsize=25)
        plt.ylabel('Q fraction',fontsize=25)
        plt.tick_params(labelsize=25)
        plt.subplot(232)
        plt.plot(rring,nring/ntheory)
        plt.xlabel('R/Re',fontsize=25)
        plt.ylabel('Percentage',fontsize=25)
        plt.tick_params(labelsize=25)
        ax=plt.subplot(233)
        cb=plt.imshow(out[i][1])
        cb.set_label('quiscent probability')
        plt.colorbar(cb)
        plt.xlim(0,len(out[i][1])-1)
        plt.ylim(0,len(out[i][1])-1)
        print('test3')
        image1(Tar,ax,len(out[i][1]))
        plt.xlabel('X',fontsize=25)
        plt.ylabel('Y',fontsize=25)
        plt.tick_params(labelsize=25)
        ax=plt.subplot(234)
        print('test4')
        image2(Tar,0,0,ax,'red')
        plt.tick_params(labelsize=25)
        plt.subplot(235)
        print('test5')
        cb=plt.imshow(out[i][3])
        plt.xlim(0,len(out[i][1])-1)
        plt.ylim(0,len(out[i][1])-1)
        plt.colorbar(cb)
        plt.xlabel('ha_ew',fontsize=25)
        plt.tick_params(labelsize=25)
        plt.tight_layout()
        plt.savefig('./pics0830/D/'+str(drp[1].data[Tar][2])+'.png',dpi=400)
        plt.close()
        i+=1

def output(out):
    
    
    L=len(out)
    print('output begin, total='+str(int(L))+'\n')
    tar=np.array([])
    r_up=np.array([])
    r_down=np.array([])
    p_in=np.array([])
    r_in=np.array([])
    p_out=np.array([])
    r_out=np.array([])
    
    i=0
    while i<L:
        tar=np.append(tar,out[i][0])
        r_up=np.append(r_up,out[i][1])
        r_down=np.append(r_down,out[i][2])
        p_in=np.append(p_in,out[i][3])
        r_in=np.append(r_in,out[i][4])
        p_out=np.append(p_out,out[i][5])
        r_out=np.append(r_out,out[i][6])
        
        i+=1
    dat=Table()
    dat['tar']=tar
    dat['r_up']=r_up
    dat['r_down']=r_down
    dat['p_in']=p_in
    dat['r_in']=r_in
    dat['p_out']=p_out
    dat['r_out']=r_out
    dat.write('./data_age_re_as_q_0613_final.fits')
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
print(len(data))
cores = mp.cpu_count()


pool = mp.Pool(processes=8)

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
print('begin output')
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
