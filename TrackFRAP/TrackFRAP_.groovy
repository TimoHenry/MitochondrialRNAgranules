//***********************************************************************************************************
//Static variables
//***********************************************************************************************************
#@ImagePlus imp	//Get current image
#@RoiManager rm	//Get Roi manager

//Get TrackMate parameters 
#@Double ( label="Spot Radius (um)", value=0.2) spotRadius
#@Double ( label="Spot Threshold", value=0.1) spotThr
#@Double ( label="Spot Quality", value=0) spotQ
#@Integer ( label="Minimum Track Duration (frames)", value=50 ) minFrames
#@Integer ( label="Max Frame Gap (frames)", value=3 ) frameGap
#@Double ( label="Frame Gap closing (um)", value=0.8 ) gapClosing
#@Double ( label="Tracking diameter (pixels)", value=10 ) size

//Get number of references defined in ROI manager
#@int ( label="How many references ?", value=3 ) Nb_ref

//Get data saving parameters
#@string ( label="Outuput folder", value="Z:/SHARED/Timo/LSM700/Emilie/Analysis" ) save_path
#@string ( label="Label (YearMonthDay_else)", value="YearMonthDay" ) DateOfToday

//***********************************************************************************************************
// Librairies import
//***********************************************************************************************************
import ij.IJ
import ij.ImagePlus
import ij.gui.Roi
import ij.plugin.*
import ij.process.FloatPolygon
import ij.measure.Calibration
import ij.gui.OvalRoi
import ij.gui.PointRoi
import ij.plugin.frame.RoiManager

import ij.process.ImageProcessor
import ij.process.ImageStatistics as IS

import fiji.plugin.trackmate.Model
import fiji.plugin.trackmate.Logger
import fiji.plugin.trackmate.Settings
import fiji.plugin.trackmate.TrackMate
import fiji.plugin.trackmate.SelectionModel

import fiji.plugin.trackmate.detection.LogDetectorFactory
import fiji.plugin.trackmate.tracking.sparselap.SparseLAPTrackerFactory
import fiji.plugin.trackmate.tracking.LAPUtils
import fiji.plugin.trackmate.features.FeatureFilter
import fiji.plugin.trackmate.features.track.TrackDurationAnalyzer
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer
import fiji.plugin.trackmate.Spot
import fiji.plugin.trackmate.action.CloseGapsByLinearInterpolationAction

import net.sf.ij_plugins.clustering.*

import java.io.File
import java.io.IOException
import javax.swing.JFileChooser

//For MetaParser Class 
import loci.common.DateTools
import loci.common.services.ServiceFactory
import loci.formats.ImageReader
import loci.formats.services.OMEXMLService
import ome.units.quantity.Time
import ome.units.UNITS

//***********************************************************************************************************
//Gather image information
//***********************************************************************************************************
def N			= imp.getNFrames();
def FileName	= imp.getTitle();
def infos 		= imp.getOriginalFileInfo();
File id 		= new File(infos.directory + FileName);
def timestamps 	= MetaParser.getTimeStamps(id);

//***********************************************************************************************************
//Roi management
//***********************************************************************************************************
// Get Frap ROI
rm.select(0);
rm.runCommand("Rename", "ROI BLEACH");
def roi = rm.getRoi(0);

//Get Ref(s) ROI(s)
def roi_ref = [];
for(k=0; k<Nb_ref; k++){
rm.select(k+1);
rm.runCommand("Rename", "REFERENCE"+(k+1).toString());
roi_ref[k] = rm.getRoi(k+1);
}

//Get Background ROI
rm.select(Nb_ref+1);
rm.runCommand("Rename", "BACKGROUND");
def roi_back = rm.getRoi(Nb_ref+1);

//***********************************************************************************************************
//MEASUREMENTS 
//***********************************************************************************************************

//##################################
//##### BACKGROUND measurements#####
//##################################
def Ib = []; //definition of the arry that will contain collected background intensities 
IJ.run("Set Measurements...", "mean min display redirect=None decimal=3"); //set measurement for measure command

rm.reset();												//reset roi manager before entering for loop
for(i=1; i<N+1; i++){									//For each slice of the stack (N slices, first index is 1)
	imp.setSlice(i);									//set slice
	rm.addRoi(roi_back);								//set ROI
	IJ.run(imp, "Measure", ""); 						//measure
	stats = imp.getProcessor().getStatistics();			//get measures
	Ib[i-1] = stats.mean; //starts at 0					//store measures in Ib (first index is 0)
}

//##################################
//### FRAP recovery measurements ###
//##################################
frap_roi_array = []; //definition of the ROIs after the tracking
If = []; //definition of the array that will contain collected FRAP intensities
IJ.run("Set Measurements...", "mean min display redirect=None decimal=3");

rm.reset();
rm.addRoi(roi);

def tf = new TrackFRAP(imp);
tf.findRois(size); //before default value was 11, now is dynamic
tf.trackRois(spotRadius, spotThr, spotQ, frameGap, gapClosing, minFrames);

for(i=1; i<N+1; i++){
	rm.select(i);
	IJ.run(imp, "Measure", ""); 
	stats = imp.getProcessor().getStatistics();
	If[i-1] = stats.mean;
	frap_roi_array[i] = rm.getRoi(i);
}
//##################################
//######## REF measurements ########
//##################################
def Ir_all  = new double[Nb_ref][N]; //definition of the array that will contain collected REF intensities FOR EACH REF
def ref_roi_array  = new Roi[Nb_ref][N]; //definition of the ROIs after the tracking FOR EACH REF
Ir=[]; //definition of the array that will contain collected mean of all references
sum=0; //definition of variable sum


// Cleanup roimanager, keep only first roi
for(k=0; k<Nb_ref; k++){
	rm.reset();
	rm.addRoi(roi_ref[k]);

	def tf_ref = new TrackFRAP(imp)

	tf_ref.findRois(size);
	tf_ref.trackRois(spotRadius, spotThr, spotQ, frameGap, gapClosing, minFrames);

	// Data collection 
	for(i=0; i<N; i++){
		rm.select(i+1);
		IJ.run(imp, "Measure", ""); 
		stats = imp.getProcessor().getStatistics();
		Ir_all[k][i] = stats.mean; //starts at 0
		ref_roi_array[k][i] = rm.getRoi(i+1)
	}
}

// Compute the mean over all references
for(i=0; i<N; i++){
	for(k=0; k<Nb_ref; k++){
		sum = sum + Ir_all[k][i];
	}
	Ir[i] = sum / Nb_ref;
	sum=0;
}

//********************************************************************************
//Class/methods/function definition 
//********************************************************************************

// TrackFRAP class ***************************************************************
public class TrackFRAP{

	private ImagePlus frap
	private RoiManager rm
	private ArrayList<Roi> frap_rois
	private start_points = []
	private Calibration cal
	private Model model
	private double size
	
	/**
	 * Constructor needs currently open image and start ROI
	 */ 
	public TrackFRAP(ImagePlus imp) {
		this.frap = imp
		this.cal = imp.getCalibration()
		this.rm = RoiManager.getInstance()
		this.frap_rois = rm.getRoisAsArray() as List
		
	}

	/**
	 * This puts the ROIs into the ROI manager
	 */
	public void findRois( double size ) {
		rm.reset()
		
		this.size = size
		frap_rois.eachWithIndex{ roi, idx ->

			def x = roi.getBounds().x + roi.getBounds().width / 2
			def y = roi.getBounds().y + roi.getBounds().height / 2

			// Add calibrated points
			start_points.add([x:cal.getX(x), y:cal.getX(y)])
			
			def aRoi = new OvalRoi(x-this.size/2, y-this.size/2, this.size, this.size)
			aRoi.setName("FRAP Region #"+IJ.pad(idx+1,2))
			rm.addRoi(aRoi)
		}
	}

	public void trackRois(double spotRadius, double spotThr, double spotQ, int frameGap, double gapClosing, int minFrames) {
		// Track them using TrackMate
		model = new Model();
		model.setLogger(Logger.IJ_LOGGER);
		def settings = new Settings();

		frap.killRoi();
		settings.setFrom(frap);
       
		settings.detectorFactory = new LogDetectorFactory();
		settings.detectorSettings = [ 
    		'DO_SUBPIXEL_LOCALIZATION': true,
    		'RADIUS': (double) spotRadius,
    		'TARGET_CHANNEL': (int) 0,
    		'THRESHOLD': (double) spotThr,
    		'DO_MEDIAN_FILTERING': false
		];

		//Configure spot filters - Classical filter on quality
		def filter1 = new FeatureFilter('QUALITY', spotQ, true);
		settings.addSpotFilter(filter1);

		// Configure tracker
		settings.trackerFactory =  new SparseLAPTrackerFactory() ;
		settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap(); // almost good enough
		settings.trackerSettings['ALLOW_TRACK_SPLITTING']   = false;
		settings.trackerSettings['ALLOW_TRACK_MERGING']     = false;
		settings.trackerSettings['MAX_FRAME_GAP']           = (int) frameGap;
		settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE']= (double) gapClosing;
		settings.trackerSettings['LINKING_MAX_DISTANCE']    = (double) gapClosing;
		
		settings.addTrackAnalyzer(new TrackDurationAnalyzer());
		def filter2 = new FeatureFilter('TRACK_DURATION', minFrames* cal.frameInterval, true);
		settings.addTrackFilter(filter2);
		

		def trackmate = new TrackMate(model, settings);
		
		if (trackmate.checkInput()) trackmate.process()

		// Close gaps by linear interpolation
		def cgblia = new CloseGapsByLinearInterpolationAction();
		cgblia.execute( trackmate );


		def selectionModel = new SelectionModel(model);
		def displayer =  new HyperStackDisplayer(model, selectionModel, frap);

		// find closest points to ROIs in frame one
		def spots = model.getSpots();
		
		def trackIds = new ArrayList<Integer>();

		start_points.eachWithIndex{ point, i -> 
			
			def ax = point['x']
			def ay = point['y']

			def closest = spots.getClosestSpot(new Spot(ax, ay, 0, 5.0, 10), 0, true);
			
			def trackID = model.getTrackModel().trackIDOf(closest);
			
			trackIds.add(trackID);
			trackToRois(trackID, "FRAP Region", i, this.size);
		}

		// now repeat for all tracks that ARE NOT the previous ones
		model.getTrackModel().trackIDs(true).each{ id ->
			if(!trackIds.contains(id)) {
				trackToRois(id, "Other Region", id, 2 * cal.getRawX(spotRadius));
			}
		}

		model.getLogger().log(model.toString())
	}

	private void trackToRois(int trackID, String name, int id, double width) {
		// get all the spots
		def trackSpots = model.getTrackModel().trackSpots(trackID);
		def sorted = new ArrayList< Spot >( trackSpots );
		def comparator = Spot.frameComparator;
       	Collections.sort(sorted, comparator);
		
			
		// for each spot, get the mean intensty in the original ROI size
		sorted.each{
			// Get the position in XY T
			print(it.getFeatures());
			def x = cal.getRawX(it.getFeature("POSITION_X"));
			def y = cal.getRawY(it.getFeature("POSITION_Y"));
			def t = Math.round(it.getFeature("POSITION_T") / cal.frameInterval);
			
			def r = new OvalRoi(x-width/2,y-width/2,width,width);
			r.setPosition((int) t+1);
			r.setName(name+" #"+(id+1)+" frame "+IJ.pad((int)t,2));
			rm.addRoi(r);
		}
	}
}


// Metaparser class ***************************************************************
// Can contain tools to parse metadata nicely
//https://github.com/openmicroscopy/bioformats/blob/develop/components/formats-gpl/utils/PrintTimestamps.java
public class MetaParser {

	public static List<Double> getTimeStamps(File image_file) {
		return getTimeStamps(image_file, 0)
	}
	
	// This returns the timestamps for the first channel and first z of each timepoint for the given sereis in the given file
	public static List<Double> getTimeStamps(File image_file, Integer series_number) {

		// the Bioformats way for getting metadata from a file
		def service = new ServiceFactory().getInstance(OMEXMLService.class)
		def meta = service.createOMEXMLMetadata()

		def reader = new ImageReader()
		reader.setMetadataStore(meta)
		reader.setId(image_file.getAbsolutePath())

		// Get the data from this series only (useful for things like LIF files)
		reader.setSeries(series_number)
		
		def timestamps = []
		// Get the total number of planes
    	def planeCount = meta.getPlaneCount(series_number)
    	(0..planeCount-1).each{plane ->
    		Time deltaT = meta.getPlaneDeltaT(series_number, plane)
      		if (deltaT != null) {
				// convert plane ZCT coordinates into image plane index
				def z = meta.getPlaneTheZ(series_number, plane).getValue().intValue()
				def c = meta.getPlaneTheC(series_number, plane).getValue().intValue()
				def t = meta.getPlaneTheT(series_number, plane).getValue().intValue()
				if (z == 0 && c == 0) {
					timestamps.add(deltaT.value(UNITS.SECOND).doubleValue())
				}
      		}
    	}
    	
    	return timestamps
	}
}


//***********************************************************************************************************
// Show only Rois of interest to follow tracking
//***********************************************************************************************************
rm.reset()

// 1) Add FRAP roi
rm.addRoi(roi)
// 2) Add tracking of frap
for(i=1; i<frap_roi_array.size(); i++){
rm.addRoi(frap_roi_array[i])}
// 3) Add ref ROI + tracking of ref (for each ref)
for(k=0; k<Nb_ref; k++){ 
	rm.addRoi(roi_ref[k])
	for(i=0; i<N; i++){
	rm.addRoi(ref_roi_array[k][i])
	}
}
// 4) Add background Roi
rm.addRoi(roi_back)

// For better visibility of trajectories we remove the overlay
IJ.run("Remove Overlay", "");

//***********************************************************************************************************
// Create .txt file (output format compatible with FRAPanalyzer)
//***********************************************************************************************************
def myFile = new File(save_path, FileName+DateOfToday+"_intensities_"+".txt");
myFile.text = '' // empty file in the case if already exxisted, otherwise the measurements add
	
myFile << "Frame \t Time \t Intensity_FRAP \t Intensity_REF \t Intensity_BACK \r\n" ;
for(o=0; o<N; o++) {
	//print("ligne numéro: "+o)
	myFile <<  o+1 + "\t" + timestamps[o] + "\t" + If[o] + "\t" + Ir[o] +  "\t" + Ib[o] + "\r\n" ;
};
	
myFile.createNewFile();

//***********************************************************************************************************
// Create .txt file (output format compatible with FRAPanalyzer) with METADATA
//***********************************************************************************************************
def myFile2 = new File(save_path, FileName+DateOfToday+"_metadata_"+".txt");
myFile2.text = '' // empty file in the case if already exxisted, otherwise the measurements add

myFile2 << "METADATA FILE LINKED TO THE FILE .txt NAMED: " + myFile.getName() + "\r\n" + "\r\n"

myFile2 << "TrackFRAP PARAMETERS:" + "\r\n" 
myFile2 << "Spot Radius (um): "+ spotRadius + "\r\n"
myFile2 << "Spot Threshold: "+ spotThr + "\r\n"
myFile2 << "Spot Quality: "+ spotQ + "\r\n"
myFile2 << "Minimum Track Duration (frames): " + minFrames + "\r\n"
myFile2 << "Max Frame Gap (frames): " + frameGap + "\r\n"
myFile2 << "Frame Gap closing (um): " + gapClosing + "\r\n"
myFile2 << "Tracking diameter (pixels): " + size + "\r\n"
myFile2 << "\r\n"

myFile2 << "Initial ROIs" + "\r\n"
myFile2 << "Roi are presented as ROI java objects. Rectangle is the shape, x and y are coordinates of the up left-hand corner of the rectangle and the width and the height of the rectangle" + "\r\n"
myFile2 <<  "FRAP ROI: " + roi +  "\r\n"
myFile2 << "How many references ?: "+ Nb_ref + "\r\n"
for(r=0; r<Nb_ref; r++) {	
	myFile2 << "Reference number " + (r+1) +": " + roi_ref[r] + "\r\n"
};
myFile2 <<  "Background ROI: " + roi_back +  "\r\n"

myFile2.createNewFile();