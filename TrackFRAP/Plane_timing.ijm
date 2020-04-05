//*********************************************************************************************
// Opening
//*********************************************************************************************
//chose a file via a file manager
my_file = File.openDialog("Choose a file");

//Use BioFormats to open file and import metadata
run("Bio-Formats Importer", "open=my_file autoscale color_mode=Grayscale display_rois rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
title = getTitle() //get image name
//close("DuplicataOf"+title); //in case one was already open 

//for better vizualization in log if several opening happen...
print("\r\n")
print("************* BIO FORMAT IMPORT ************")
print("Image named "+title+" was open, congratulations \r\n")

//*********************************************************************************************
// ROI managing
//*********************************************************************************************
//Select FRAP ROI and rename it
roiManager("Select", 0); 
roiManager("Rename", "ROI BLEACH");
//Select second ROI imported by boiformat and delete it
roiManager("Select", 1);
roiManager("Delete");
//deselect all ROI
roiManager("Deselect");

//*********************************************************************************************
// create copy with time stamping (BEWARE, created image is RGB = cannot be analyzed)
//*********************************************************************************************

run("Duplicate...", "title=DuplicataOf"+title+" duplicate"); 


selectWindow("DuplicataOf"+title);
run("Green"); //Set Green LUT to stack
	
run("Bio-Formats Macro Extensions"); 
	
Ext.setId(my_file);
Ext.getSeriesCount(seriesCount);
print("Series Count = " + seriesCount);
print("\r\n");

//Get metadata for each series and print them to log window
// Based on : https://www.researchgate.net/post/How_to_use_the_time_information_stored_in_the_frames_of_a_video_in_ImageJ_instead_of_individually_labeling_the_stacks
for (s=0; s<seriesCount; s++) {
	Ext.setSeries(s);
	Ext.getSizeX(sizeX);
	Ext.getSizeY(sizeY);
	Ext.getSizeZ(sizeZ);
	Ext.getSizeC(sizeC);
	Ext.getSizeT(sizeT);
	print("Series #" + s + ": image resolution is " + sizeX + " x " + sizeY + " pixels");
	print("Focal plane count = " + sizeZ);
	print("Channel count = " + sizeC);
	print("Time point count = " + sizeT);

	Ext.getImageCount(imageCount);
	print("Plane count: " + imageCount);

	creationDate = "";
	Ext.getImageCreationDate(creationDate);
	print("Creation date: " + creationDate);

run("RGB Color");
run("Colors...", "foreground=red background=green selection=orange");

deltaT = newArray(imageCount); //new array size of number of images in stack

for (no=0; no<imageCount; no++) {
	Ext.getPlaneTimingDeltaT(deltaT[no], no);
	if (deltaT[no] == deltaT[no]) { // not NaN
		slice=no+1;
		stamp=""+deltaT[no]+"(s)";
		xpos=sizeX-55;
		ypos=sizeY-5;
		para="format=Text starting=0 interval=1 x="+xpos+" y="+ypos+" font=10 text="+stamp+" range="+slice+"-"+slice;
		run("Label...", para);
	}
}
print(deltaT[1] + " seconds/frame");
}

//*********************************************************************************************
// Close original
//*********************************************************************************************
close(title);
