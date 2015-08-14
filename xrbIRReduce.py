import glob, sys, os
import numpy as np 
import pandas as pd
from astropy.io import fits
import alipy

def redWrap(color):
	fileListing=pd.read_csv('HrawList.csv')
	grDate=fileListing.groupby('Date')

	for name, group in grDate:
		fileList=group['file'].values

		flat=nearestFlat(group['JulianDate'].values, color)
		subFlat(fileList,flat)

		align(map(lambda x: 'flat_sky_'+x[2:], fileList))

		alignedFiles=[]
		for aligned in map(lambda x: 'flat_sky_'+x[2:-5]+'_gregister.fits',fileList):
			if os.path.isfile(aligned):
				alignedFiles.append(aligned)
		if len(alignedFiles) >= 3:
			combine(alignedFiles)
		else:
			print "not enough files were aligned"
		#combine(map(lambda x: 'flat_sky_'+x[2:-5]+'_gregister.fits',fileList))
	return



def subFlat(fileList, flatFile):
    n=len(fileList)
    dataArray=np.zeros((n,512,512))
    headerinfo=[]
    
    hdulist=fits.open(flatFile)
    flatData=hdulist[0].data
    hdulist.close()
    flatData[flatData < 1.0]=1.0
    
    norm_flat_data=flatData/np.mean(flatData)
    #hdu=fits.PrimaryHDU(norm_flat_data)
    #hdu.writeto('norm_flat_py.fits')
    
    for i in xrange(0,n):
        hdulist=fits.open(fileList[i])
        dataArray[i]=hdulist[0].data
        headerinfo.append(hdulist[0].header)
        
    sky=np.median(dataArray, axis=0)
    
    #hdu=fits.PrimaryHDU(sky)
    #hdu.writeto(str(date)+'sky_py.fits')
    
    skySub=np.zeros((n,512,512))
    skySubFlat=np.zeros((n,512,512))
    
    for i in xrange(0,n):
        skySub[i]=dataArray[i]-sky
        hdu=fits.PrimaryHDU(skySub[i])
        #hdu.header=headerinfo[i]
        #hdu.writeto('sky_'+fileList[i])
        
        skySubFlat[i]=skySub[i]/norm_flat_data
        hdu=fits.PrimaryHDU(skySubFlat[i])
        hdu.header=headerinfo[i]
        hdu.writeto('flat_sky_'+fileList[i][2:])
        
    return

def align(fileList):
	ref_image=fileList[0]

	#get a listing of the path for each image for a given filter

	#everything below here i copied from the alipy demo http://obswww.unige.ch/~tewes/alipy/tutorial.html
	identifications = alipy.ident.run(ref_image, fileList, visu=False)
	# That's it !
	# Put visu=True to get visualizations in form of png files (nice but much slower)
	# On multi-extension data, you will want to specify the hdu (see API doc).

	# The output is a list of Identification objects, which contain the transforms :
	for id in identifications: # list of the same length as images_to_align.
		if id.ok == True: # i.e., if it worked

			print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
			# id.trans is a alipy.star.SimpleTransform object. Instead of printing it out as a string,
			# you can directly access its parameters :
			#print id.trans.v # the raw data, [r*cos(theta)  r*sin(theta)  r*shift_x  r*shift_y]
			#print id.trans.matrixform()
			#print id.trans.inverse() # this returns a new SimpleTransform object

		else:
			print "%20s : no transformation found !" % (id.ukn.name)

	# Minimal example of how to align images :
	outputshape = alipy.align.shape(ref_image)
	# This is simply a tuple (width, height)... you could specify any other shape.
	for id in identifications:
	    if id.ok == True:
	        # Variant 2, using geomap/gregister, correcting also for distortions :
	            alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars,verbose=False,
	                                      shape=outputshape, makepng=False, outdir=".")

	return


def combine(fileList):
    n=len(fileList)
    dataArray=np.zeros((n,512,512))
    #headerinfo=[]
    
    #date=fileList[0][5:10]
    for i in xrange(0,n):
        hdulist=fits.open(fileList[i])
        dataArray[i]=hdulist[0].data
        #headerinfo.append(hdulist[0].header)
        
    combined=np.sum(dataArray, axis=0)
    hdu=fits.PrimaryHDU(combined)
    #hdu.header=headerinfo[i]
    hdu.writeto('final_'+fileList[0].split('_')[-2].split('.')[0]+'.fits')
    return

def nearestFlat(JDlist,color):
	'''
	given the dfview, find the median observation time, and using this time
	identify the jflat that was observed closest in time
	return a string that gives the path to this flat
	'''
	#find the median observation time of the dfView
	meanTime=np.nanmean(JDlist)
	#load up the flat df
	flatDF=pd.read_csv(color.upper()+'flatList.csv')
	#subtract the median time form all the flat times, take absolute value, 
	#find index of min, use this index to grab the file name with this observation date
	return flatDF.file[np.argmin(np.abs(flatDF.JulianDate.values.astype(float) - meanTime))]


