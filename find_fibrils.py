#!/usr/bin/python

import glob
import sys
import matplotlib.pyplot as plt
from scipy import fftpack
from scipy import misc
import numpy as np
import pylab as py
import struct
from numpy import *
import math
import os

vers=0.3

outdir = 'output'
cutoff = 1.4



#------------- READ MRCs--------------
# from David Stokes
class mrc_image:
    def __init__(self,filename):
        self.numbytes1=56           # 56 long ints
        self.numbytes2=80*10          # 10 header lines of 80 chars each
        self.filename=filename

def read(self):
    input_image=open(self.filename,'rb')
    self.header1=input_image.read(self.numbytes1*4)
    self.header2=input_image.read(self.numbytes2)
    byte_pattern='=' + 'l'*self.numbytes1   #'=' required to get machine independent standard size
    self.dim=struct.unpack(byte_pattern,self.header1)[:3]   #(dimx,dimy,dimz)
    self.imagetype=struct.unpack(byte_pattern,self.header1)[3]  #0: 8-bit signed, 1:16-bit signed, 2: 32-bit float, 6: unsigned 16-bit (non-std)
    if (self.imagetype == 0):
        imtype='b'
    elif (self.imagetype == 1):
        imtype='h'
    elif (self.imagetype == 2):
        imtype='f4'
    elif (self.imagetype == 6):
        imtype='H'
    else:
        imtype='unknown'   #should put a fail here
    input_image_dimension=(self.dim[1],self.dim[0])  #2D images assumed
    self.image_data=fromfile(file=input_image,dtype=imtype,count=self.dim[0]*self.dim[1]).reshape(input_image_dimension)
    input_image.close()
    return(self.image_data)

def bin_ndarray(ndarray, new_shape, operation='mean'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.
    """
    operation = operation.lower()
    if not operation in ['sum', 'mean']:
        raise ValueError("Operation not supported.")
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d,c in zip(new_shape,
                                                  ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        op = getattr(ndarray, operation)
        ndarray = op(-1*(i+1))
    return ndarray

def get_PS(image_in,outname):
    mrcim = mrc_image(image_in)
    image = read(mrcim)
    #image = plt.imread(image_in)           # for reading from tiffs instead of mrcs

    binned = bin_ndarray(image,new_shape=(image.shape[0]/2,image.shape[1]/2),operation='mean')

    # fourier transform of the image.
    F1 = fftpack.fft2(binned)
    #F1 = fftpack.fft2(image)

    # make pretty - put low res in center
    F2 = fftpack.fftshift( F1 )
     
    # Calculate 2D power spectrum
    ps2D = np.abs( F2 )**2
    
    return (ps2D)


def get_radialpoint(angle,radius):
    x = round(radius*math.cos(math.radians(angle)))
    y = round(radius*math.sin(math.radians(angle)))
    return(x,y)


def get_data(file):
    filename = file.split('.')[0]
    the_PS = get_PS(file,'{0}_fullPS.png'.format(filename))
    cx,cy = (int(the_PS.shape[0]/2)+1,int(the_PS.shape[1]/2)+1)
    return(the_PS,cx,cy)


def calc_dist(psdata,pxsize,resolution):
    """ given a resolution calculate its frequency
    on the PS
    """
    nyquist = 2*pxsize
    distance = ((np.shape(psdata)[0]*.5)*nyquist)/resolution
    return(distance)

def do_radii_calcs(thedata,cx,cy):
    allradii= []
    allamyloid = []
    list1 = []
    
    amyloid_data=[]
    amyloid_repeat=int(calc_dist(thedata,1.07,4.8))
    amyloid_radaii = range(amyloid_repeat-60,amyloid_repeat+60)
    
    for angle in (angles):
        for i in radii:
            pt= get_radialpoint(angle,i)
            list1.append(thedata.item(int(cx+pt[0]),int(cy+pt[1])))  
            thedata[int(cx+pt[0]),int(cy+pt[1])] = thedata[int(cx+pt[0]),int(cy+pt[1])]*0.25
        
        #for i in amyloid_radaii:                                   # code for getting 4.8A data as well
        #    pt= get_radialpoint(angle,i)                           # doesn't work very well - too much noise
        #    #list1.append(np.log10(middle.item(cx+pt[0],cy+pt[1])))   #removing this speeds it up.
        #    amyloid_data.append(thedata.item(cx+pt[0],cy+pt[1]))  
        #    thedata[cx+pt[0],cy+pt[1]] = 1
        
        allradii.append(np.mean(list1))
        #allamyloid.append(np.mean(amyloid_data))
        list1= []
        amyloid_data = []
    
    x = np.array(allradii)
    a = np.array(allamyloid)
    
    return(x,a)


def moving_average(a, n) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n



angles = range(0,1441)
angs = []
n=0
for i in angles:
    angs.append(n)
    n+=0.25
angles = angs

#remove 5 degrees around 0, 90, 180, and 270 degree sampling because of edge artifacts
del angles[1420:1441]
del angles[1060:1101]
del angles[700:741]
del angles[340:381]
del angles[0:21]
radii = range(40,90)


cwd = os.getcwd()

if os.path.isdir('{0}/{1}'.format(cwd,outdir)) == False:
    os.system('mkdir {0}/output'.format(cwd))


files = sys.argv[1:]
fibrilsout = open('{0}/{1}/probably_fibrils.txt'.format(cwd,outdir),'w')
emptyout = open('{0}/{1}/probably_not_fibrils.txt'.format(cwd,outdir),'w')




runcount = float(0)
fibcount = 0
emptycount = 0
pthresh = 10
print('running on {} files'.format(len(files)))
for i in files:
    if int((runcount/len(files))*100) > pthresh:
        print('{0}% done'.format(pthresh))
        pthresh +=20
    count = []
    filename = i.split('/')[-1].split('.')[0]
    ps_chunk,cx,cy=get_data(i)
    lowres,amyloid_peak = do_radii_calcs(ps_chunk,cx,cy)
    lowres = moving_average(lowres,50)
    normed = [j/min(lowres) for j in lowres]
    cutoffvals = np.full((len(lowres),1),cutoff)

    #amyloid_peak=moving_average(amyloid_peak,50)                       # Code for checking at 4.8 Angstroms as well, 
    #amyloid_normed = [j/min(amyloid_peak) for j in amyloid_peak]       # nice idea, but doesn't seem to give good results
    #amyloid_shift = [j+360 for j in np.arange(len(amyloid_normed))]    # too much noise
    
    for j in normed:
        if j > cutoff:
            count.append(j)
    totalscore = np.sum(count)
    if totalscore > 1:
        fibrilsout.write('{0}\n'.format(filename))
        fibcount +=1
    else:
        emptyout.write('{0}\n'.format(filename))
        emptycount +=1
    plt.plot(np.arange(len(normed)),normed,label='{0}_lores'.format(i.split('.')[0]))
    plt.plot(np.arange(len(normed)),cutoffvals)
    #plt.plot(amyloid_shift,amyloid_normed,label='{0}_4.8'.format(i.split('.')[0]))#
    plt.legend(loc='best',fontsize='xx-small')
    plt.savefig("{2}/{1}/{0}_output.png".format(filename,outdir,cwd))
    plt.close()
    
    #fig,ax = plt.subplots()
    #ax.imshow(np.log10(ps_chunk[cx-100:cx+100,cy-100:cy+100]),cmap='Greys')       # print out the PS with the regions sampled                
    #ax.axis('off')
    #fig.savefig("{2}/{1}/{0}_PS.png".format(filename,outdir,cwd))
    #plt.close()
    #
    #fig,ax = plt.subplots()
    #mrcim = mrc_image(i)
    #image = read(mrcim)
    #ax.imshow(image,cmap='Greys') 
    #ax.axis('off')
    #fig.savefig("{2}/{1}/{0}_img.png".format(filename,outdir,cwd))
    #plt.close()
    
    runcount +=1
print('--\n{0} predicted fibrils\n{1} predicted empty'.format(fibcount,emptycount))
