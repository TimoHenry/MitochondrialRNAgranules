# FRAP analysis for moving spots

Mitochondrial RNA granules (MRGs) are moving structures within mitochondria organelles inside living cells. The TrackFRAP script was developed to automatized the data collection of Fluorescence Recovery After Photobleaching (FRAP) experiments videos by tracking the intensity recovery of MRGs after FRAP experiment. This script could be extended to other applications. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

## Installing

Download directly from the git repository. You can clone it directly in Fiji Marcos folder to access the macros from the Fiji menu.

```

https://c4science.ch/source/TrackFRAP.git
ssh://git@c4science.ch/source/TrackFRAP.git

```

## Goals and input/output formats



### Open_Bioformats.ijm input/output
	
- **What it does**: open the file with bioformats metadata
- **Input**: .lsm format FRAP movie
- **Output**: ImagePlus created in imageJ/Fiji with FRAP area defined in ROI manager

### TrackFRAP.groovy input/output

- **What it does**: Track the bleached moving spot and measure the recovery
- **Input**: a FRAP .lsm format movie
- **Output**: TrackFRAP output is a .txt file with measurement of time and mean intensity of calculated ROI tracking spots.	The TrackFRAP output.txt file is compatible with FRAP analyser software from ActinSim (http://actinsim.uni.lu/eng/Downloads/FRAPAnalyser).

### Plane_timing.ijm input/output

- **What it does**: print timing of the video on the frames
- **Input**: a .lsm FRAP movie
- **Output**: RGB copy of the input with timings printed on frames


### Simple_FRAP_ROI_Mean.py input/output

- **What it does**: measures recovery of the defined ROI
- **Input**: .lsm format FRAP movie
- **Output**: .txt file

### Data_Analysis.ipynb

- **What it does**: fit frap recovery data to exponential recovery
- **Input**: .txt files generated with trackFRAP of Simple_ROI_MEAN.py macros
- **Output**: graphs and half times and mobile fraction











## Example of use

In the folder named example_data you will find 3 FRAP videos that can be used for the following tutorial on how to use TrackFRAP analysis. The macro uses an input the .lsm extension format created when taking FRAP videos with LSM Zeiss microscopes.

First you will need to open the video you want to analyse using the Open BioFormats macro. Just run the macro, chose your file. The video will open as an ImagePlus object, the log will congratulate you and the bleach area will be available from the Fiji ROI manager.

{F9936356, layout=center, size=full, alt="output"} 


Once this step is done you will need to prepare two things because computation of the recovery relies on three measurements: the frap area recovery measurement, some references measurements, and the background measurement. Thanks to the Open Bioformats macro, the bleaching area is already in the ROI manager. So, first, you will need to add to the ROI manager the ROIs that englobe or are near the references you also want to track. The shape has no importance because the trackFRAP macro will take the closest spot and adapt the ROI for measurements. Then, after all the references, you will need to add to the ROI manager the background ROI. Here the shape of the ROI you set is important because it will represent the measured area overt time. Names have no influence.

{F9936893, layout=center, size=full, alt="preparation"} 

You are now ready to run TrackFRAP. Run it and enter the wanted parameters: 
- **Spot radius**: radius in um of the spot you want to track. 
- **Spot Threshold**: higher threshold reduces number of detected spots. You can put it higher if the spots are really bright, but put it as low as possible if the spots are not that bright (minimum is 0)
- **Spot quality**: Similar to threshold, can be a bit higher if
- **Minimum track duration**:
- **Max Frame gap**: maximum gap of frames authorized
- **Frame gap closing**: maximum distance in um allowed between two frames of the same tracking
- **Tracking diameter**: Diameter of the ROI created for measurement around the spot detected
- **How many references? **: how many ROI for references did you put on the ROI manager? How many references do you take into account?
- **Output folder**: put here the path of the folder in which you want your .txt file to be created
- **Label**: label that will be used in the name of the –txt file to recognize it, can be further used to sort the files


Because TrackFRAP is based on TrackMate plugin, you may want to understand how the plugins and the different parameters act on detection by playing with the TrackMate plugin directly. The next image show examples of parameters that can be used with our image12.lsm that work with TrackFRAP. 

{F9937326, layout=center, size=full, alt="parameters”}

Next, you see that TrackFRAP will upload in the ROI manager the tracking ROIs, and by selecting them you can assess of the tracking acuity. The results window will display the measurements. Some warnings may display on the log window.

{F9937237, layout=center, size=full, alt="output"} 

{F9986951, layout=center, size=thumb, alt="video"} 

The final output is a .txt that will look like the one below. Separations are made by tabs.

{F9987017, layout=center, size=full, alt=".txt output"} 

When you have all the .txt files created for an experiment you can gather them in a folder and use the python jupyter notebook Data_Analysis.ipynb file to analyze them, plot and fit to the desired recovery function. 


## Authors

* **Timo REY** timo.rey@epfl.ch
* **Emilie CUILLERY** emilie.cuillery@epfl.ch
* **Olivier BURRI**

https://secure.phabricator.com/book/phabricator/article/remarkup/

https://c4science.ch/file/




