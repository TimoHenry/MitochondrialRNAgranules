//*********************************************************************************************
// Cleaning
//*********************************************************************************************
roiManager("reset"); 

//*********************************************************************************************
// Opening
//*********************************************************************************************
//chose a file via a file manager
my_file = File.openDialog("Choose a file");

//Use BioFormats to open file and import metadata
run("Bio-Formats Importer", "open=my_file autoscale color_mode=Grayscale display_rois rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
title = getTitle() //get image name

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
//Select second ROI imported by bioformat and delete it
roiManager("Select", 1);
roiManager("Delete");
//deselect all ROI
roiManager("Deselect");

//*********************************************************************************************
// LUT selection for better visualisation
//*********************************************************************************************
//Run glow LUT on analysis stack for better visualization 
selectWindow(title);
run("Fire");
roiManager("Select", 0); // select ROI FRAP (why not)