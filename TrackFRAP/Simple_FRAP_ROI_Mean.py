#Script inspired by
#https://imagej.net/Analyze_FRAP_movies_with_a_Jython_script

# How to use this script...
# 1) Load a movie
# 2) Draw and add to ROI manager
#	- ROI ref
#	- ROI background
# 3) Launch script and find intensities values written in chosen folder


##########################################################################################
# 1) Importations
##########################################################################################
import java.awt.Color as Color
from ij import WindowManager
from ij.plugin.frame import RoiManager
from ij.process import ImageStatistics
from ij.measure import Measurements
from ij import IJ
from ij.measure import CurveFitter
from ij.gui import Plot
from ij.gui import PlotWindow
import math
import os.path
from ij.io import DirectoryChooser
from ij.gui import Roi

##########################################################################################
# 2) Image handling and parameters definition
##########################################################################################
imp			= WindowManager.getCurrentImage()	# Get open movie
title		= imp.getTitle()					# Get title of movie
n_slices 	= imp.getNFrames()					# Define number of frames max
stack       = imp.getImageStack()				#
calibration = imp.getCalibration()				#get calibration for measurements

##########################################################################################
# 3) ROI management
##########################################################################################
roi_manager = RoiManager.getInstance()		# Get ROIs
roi_list    = roi_manager.getRoisAsArray()	# Store them in array of ROI

roi_FRAP    = roi_list[0];	# We assume 1st ROI is FRAP ROI
roi_REF	    = roi_list[1];	# We assume 2nd ROI is REF ROI
roi_BACK    = roi_list[2];	# We assume 3rd ROI is BACK ROI

##########################################################################################
# 4) Create bigger ROI for FRAP intensities measurement (because mito moves)
##########################################################################################
n = 3 #prefere an odd number...							#Define multiplying factor 

X = roi_FRAP.getBounds().x								#Get X the x coordinate of the FRAP ROI
Y = roi_FRAP.getBounds().y								#Get Y the y coordinate of the FRAP ROI
W = roi_FRAP.getBounds().width							#Get W width of FRAP ROI
H = roi_FRAP.getBounds().height							#Get H height of FRAP ROI

if n==2:
   MyRoi = Roi(X-W/2, Y-H/2, 2*W, 2*H)
else:
   MyRoi = Roi(X-((n-1)/2)*W, Y-((n-1)/2)*H, n*W, n*H)	#Create new ROI bigger/smaller of factor n and centered around FRAP ROI

roi_manager.addRoi(MyRoi)
##########################################################################################
# 5) Collecting intensities over time
##########################################################################################
If = []  # Definition of array that will contain intensities of FRAP ROI
Ir = []  # Definition of array that will contain intensities of REF ROI
Ib = []  # Definition of array that will contain intensities of BACK ROI

# Loop over each slice of the stack
for i in range(0, n_slices):

    # Get the current slice 
    ip = stack.getProcessor(i+1)
    
    # FRAP intensities
    ip.setRoi(MyRoi)															#Define ROI												
    stats = ImageStatistics.getStatistics(ip, Measurements.MEAN, calibration);	#Calculate mean										
    If.append(stats.mean)														#Store calculated mean in If
   
    # non-FRAPed area (ref)
    ip.setRoi(roi_REF)
    stats = ImageStatistics.getStatistics(ip, Measurements.MEAN, calibration);
    Ir.append(stats.mean)
    
    # Do the same for background
    ip.setRoi(roi_BACK)
    stats = ImageStatistics.getStatistics(ip, Measurements.MEAN, calibration);
    Ib.append(stats.mean)

##########################################################################################
# 6) Create .txt file with results (format FRAP Analyzer)
##########################################################################################
#dc = DirectoryChooser("Choose directory for saving results")	# Open file manager to chose a directory
#srcDir = dc.getDirectory()										# Get chosen directory 
srcDir = "Z:/SHARED/Timo/LSM700/Emilie/Analysis"
date="20180828"

completepath = os.path.join(srcDir, title+"_"+date+"_FreeDiff_intensities.txt")
MyFile = open(completepath, "w")
MyFile.write("Frame \t Intensity_FRAP \t Intensity_REF \t Intensity_BACK \r\n")
for o in range(0,n_slices):
	MyFile.write( str(o+1)+ "\t" + str(If[o]) + "\t" + str(Ir[o]) +  "\t" +str( Ib[o])+ "\r\n")
MyFile.close()
print("Intensities files created")
##########################################################################################
# 7) Create .txt file with METADATA (for reproduction of results)
##########################################################################################
completepath = os.path.join(srcDir, title+"_"+date+"_FreeDiff_metadata.txt")
MyFile2 = open(completepath, "w")
MyFile2.write("METADATA FILE LINKED TO THE FILE .txt NAMED:"+title+"_"+date+"_FreeDiff"+"\r\n")
MyFile2.write("\r\n")
MyFile2.write("Chosen ROIs for measurements (rectangular, in pixels): \r\n")
MyFile2.write("FRAP ROI: upper-left corner coordinates (x="+ str(roi_FRAP.getBounds().x)+",y="+ str(roi_FRAP.getBounds().y)+"); width="+str(roi_FRAP.getBounds().width)+"; Height="+str(roi_FRAP.getBounds().height)+" \r\n")
MyFile2.write("Measured frap ROI: upper-left corner coordinates (x="+ str(MyRoi.getBounds().x)+",y="+ str(MyRoi.getBounds().y)+"); width="+str(MyRoi.getBounds().width)+"; Height="+str(MyRoi.getBounds().height)+" \r\n")
MyFile2.write("Reference ROI: upper-left corner coordinates (x="+ str(roi_REF.getBounds().x)+",y="+ str(roi_REF.getBounds().y)+"); width="+str(roi_REF.getBounds().width)+"; Height="+str(roi_REF.getBounds().height)+" \r\n")
MyFile2.write("BACKGROUND ROI: upper-left corner coordinates (x="+ str(roi_BACK.getBounds().x)+",y="+ str(roi_BACK.getBounds().y)+"); width="+str(roi_BACK.getBounds().width)+"; Height="+str(roi_BACK.getBounds().height)+" \r\n")
MyFile2.close()
print("Metadata file created")