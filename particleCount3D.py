'''
Author: Laurent Guerard
Group: IMCF
Email: laurent.guerard@unibas.ch
Creation Date: Wednesday, 24th July 2019 3:37:12 pm
-----
Last Modified: Thursday, 25th July 2019 2:58:55 pm
Modified By: Laurent Guerard
-----
HISTORY:
Date        	By	Comments
------------	---	---------------------------------------------------------

25-07-192019	LG	Added a "Clear Outside" to work only on ROI

'''


# ****************************************************************************
# *                                  Imports                                 *
# ****************************************************************************

from fiji.plugin.trackmate.detection import LogDetector
from net.imglib2.img.display.imagej import ImageJFunctions

import os
import sys
import csv
import glob
from itertools import izip

from ij import IJ, ImagePlus, ImageStack, WindowManager as wm
from ij.plugin import Duplicator
from ij.plugin.frame import RoiManager
from ij.gui import PointRoi, WaitForUserDialog
from ij.measure import ResultsTable
from ij.process import ImageConverter

# Bioformats imports
from loci.plugins import BF
from loci.plugins.in import ImporterOptions

# 3DSuite imports
from mcib3d.geom import Objects3DPopulation
from mcib3d.image3d import ImageInt, ImageHandler
from mcib3d.image3d.regionGrowing import Watershed3D
from mcib3d.image3d.IterativeThresholding import TrackThreshold

# MorpholibJ imports
from inra.ijpb.morphology import Strel3D
from inra.ijpb.morphology import Morphology

# import mcib3d.image3d.regionGrowing.Watershed3D;

# ****************************************************************************
# *                                 Variables                                *
# ****************************************************************************

# ############################# #
# TRACKMATE DETECTION VARIABLES #
# ############################# #

# Radius for the spots
radius     = 0.4
# Threshold
threshold  = 500
# Do Subpixel detection
doSubpixel = True
# Apply median before detection
doMedian   = False

# ################################ #
# ITERATIVE THRESHOLDING VARIABLES #
# ################################ #

# minimum volume
volMin  = 30
# maximum volume
volMax  = 20000
# min contrast
minCont = 0
# check threshold every step
step    = 10
# minimum threshold to start aka background noise
thmin   = 0

# ****************************************************************************
# *                         FILTER OBJECTS TOUCHING Z                        *
# ****************************************************************************

filter_objects_touching_z = False


# ****************************************************************************
# *                                 Functions                                *
# ****************************************************************************


def checkForFiles(filepath):
    """Check if files are there no matter the extension

    Arguments:
        filepath {string} -- Path and name to check if exists

    Returns:
        bool -- Returns true if exists otherwise returns false
    """
    for filepath_object in glob.glob(filepath):
        if os.path.isfile(filepath_object):
            return True

    return False


def extractChannel(imp, nChannel, nFrame):
    """Extract one channel from a hyperstack

    Arguments:
        imp {imagePlus} -- ImagePlus of the hyperstack
        nChannel {int}  -- number of the channel of interest
        nFrame {int}    -- number of the frame in case

    Returns:
        {imagePlus} -- Stack with only the channel of interest
    """
    stack = imp.getImageStack()
    ch = ImageStack(imp.width, imp.height)
    for i in range(1, imp.getNSlices() + 1):
        index = imp.getStackIndex(nChannel, i, nFrame)
        ch.addSlice(str(i), stack.getProcessor(index))
    return ImagePlus("Channel " + str(nChannel), ch)


def cellDetection3D(implus, rad, thresh, subpix, med):
    """Function to detect the cells in 3D using TrackMate

    Arguments:
        implus {imagePlus} -- ImagePlus of the image to use for detection
        rad    {int}       -- Radius of the cell to detect, half the diameter
        thresh {int}       -- Intensity threshold for the detection
        subpix {bool}      -- Option for subpixel detection
        med {bool}         -- Option for median filter before detection

    Returns:
        cellCount {int} -- Number of cells found
    """
    dim = implus.getDimensions()

    implus2 = implus.duplicate()

    # Set the parameters for LogDetector
    img           = ImageJFunctions.wrap(implus2)
    interval      = img
    cal           = impWTH.getCalibration()
    # calibration = [round(cal.pixelWidth,3), round(cal.pixelHeight,3), round(cal.pixelDepth,3)]
    calibration   = [cal.pixelWidth, cal.pixelHeight, cal.pixelDepth]

    radius     = rad  # the radius is half the diameter
    threshold  = thresh
    doSubpixel = subpix
    doMedian   = med

    # print cal.pixelDepth

    # Setup spot detector (see http://javadoc.imagej.net/Fiji/fiji/plugin/trackmate/detection/LogDetector.html)
    #
    # public LogDetector(RandomAccessible<T> img,
    #            Interval interval,
    #            double[] calibration,
    #            double radius,
    #            double threshold,
    #            boolean doSubPixelLocalization,
    #            boolean doMedianFilter)

    detector = LogDetector(img, interval, calibration, radius,
                           threshold, doSubpixel, doMedian)

    # Start processing and display the results
    if detector.process():
        # Get the list of peaks found
        peaks = detector.getResult()
        # print str(len(peaks)), "peaks were found."

        # Add points to ROI manager
        rm = RoiManager.getInstance()
        if not rm:
            rm = RoiManager()
        rm.reset()

        # Loop through all the peak that were found
        for peak in peaks:
            # Print the current coordinates
            # print peak.getDoublePosition(0), peak.getDoublePosition(1), peak.getDoublePosition(2)
            # Add the current peak to the Roi manager
            roi = PointRoi(peak.getDoublePosition(0) / cal.pixelWidth,
                           peak.getDoublePosition(1) / cal.pixelHeight)
            # print peak.getDoublePosition(2)/cal.pixelDepth
            roi.setPosition(
                int(round(peak.getDoublePosition(2) / cal.pixelDepth)) + 1)
            rm.addRoi(roi)
        # Show all ROIs on the image
        # rm.runCommand(imp2, "Show All")
        cellCount = rm.getCount()
        # Close the duplicate
        implus2.changes = False
        implus2.close()
    else:
        print "The detector could not process the data."

    # create an empty black (all intensities are 0) hyperstack with the same dimensions as the original
    markerimage = IJ.createImage(
        "Markerimage", "8-bit black", dim[0], dim[1], dim[2], dim[3], dim[4])

    # set the forground color to 255 for the fill command
    # resulting image will be binary
    IJ.setForegroundColor(255, 255, 255)
    for roi in range(rm.getCount()):
        # TODO: requires A (any) image window to be shown. not clean...
        rm.select(roi)
        rm.runCommand(markerimage, "Fill")
    # markerimage.show()
    IJ.run(markerimage, "Select None", "")
    IJ.run(implus, "Select None", "")
    return markerimage

    # print(type(markerimage))


####################################################################


# ****************************************************************************
# *                                 Main Code                                *
# ****************************************************************************

# Get image and some info
imp      = wm.getCurrentImage()
if imp is None:
    IJ.log("No image open")
    sys.exit("No image open")
origInfo = imp.getOriginalFileInfo()
dir      = origInfo.directory
if dir is None:
    dir = os.path.expanduser("~/Desktop")
fileName = origInfo.fileName
imp.show()

# Check for existing filename
outFileName     = fileName
testName        = outFileName + ".*"
outFullPathTest = os.path.join(dir, testName)
count = 1

# Save a new file for each ROI
while checkForFiles(outFullPathTest):
    # print(outFullPathTest)
    outFileName     = fileName + "_" + str(count)
    testName        = outFileName + ".*"
    outFullPathTest = os.path.join(dir, testName)
    count           = count + 1

outFullPath = os.path.join(dir, outFileName)
# print(outFileName)
# sys.exit()

rm = RoiManager.getInstance()
if not rm:
    rm = RoiManager()
rm.reset()

# Awaiting ROI or quits after 5 tries
count = 0
while(imp.getRoi() is None):
    WaitForUserDialog(
        "Draw the region of interest and press OK").show()
    if count == 5:
        sys.exit("Too many clicks without ROI")
    else:
        count = count + 1

rm.addRoi(imp.getRoi())
rm.runCommand("Save", outFullPath + ".zip")
rm.select(imp, 0)

# Apply morphological filters
# TODO: Could be improved
#imp2 = Duplicator().crop(imp)
imp2 = imp.duplicate()
rm.select(imp2, 0)
IJ.run(imp2, "Clear Outside", "stack");

# IJ.run("Morphological Filters (3D)",
#        "operation=[White Top Hat] element=Ball x-radius=0.8 y-radius=0.8 z-radius=4")
# impWTH = wm.getCurrentImage()

# filter with Morphological opening to get rif of the rings
# create structuring element (ball of x,y,z-radius in px)
strel = Strel3D.Shape.BALL.fromRadiusList(1, 1, 4)

# apply morphological opening filter to input image

imStWTH = Morphology.whiteTopHat(imp2.getImageStack(), strel)

impWTH  = ImagePlus("WTH results", imStWTH)
# assign correct calibration
impWTH.setCalibration(imp2.getCalibration())
# impWTH.show()

# sys.exit()
# Not necessary because otherwise some particles are not whole
# IJ.run("Restore Selection", "")
# IJ.run(impWTH, "Clear Outside", "stack")

# Get the marker image with the peaks of cells using TrackMate
MarkerImage = cellDetection3D(
    impWTH, radius, threshold, doSubpixel, doMedian)

# Converts to 16 bit for the thresholding
imp16bit = imp2.duplicate()
ImageConverter(imp16bit).convertToGray16()
imp16bit.show()

# IJ.run(imp16bit, "16-bit", "")
# markerImp.show()

# IJ.selectWindow(imp16bit.getTitle())
# IJ.run("3D Iterative Thresholding",
#        "min_vol_pix=100 max_vol_pix=2000 min_threshold=0 min_contrast=0 criteria_method=MSER threshold_method=KMEANS segment_results=All value_method=10")
# impThresh = wm.getCurrentImage()

# init the iterative thresholding
# various methods
IT = TrackThreshold(volMin, volMax, minCont, step, step, thmin)
# check threshold every step
tmethod = TrackThreshold.THRESHOLD_METHOD_KMEANS
IT.setMethodThreshold(tmethod)
# favours the edges
cri = TrackThreshold.CRITERIA_METHOD_MSER
# cri = TrackThreshold.CRITERIA_METHOD_MIN_ELONGATION
IT.setCriteriaMethod(cri)
# segment objects
impThresh = IT.segment(imp16bit, True)
impThresh.show()


# Returns the dimensions of this image (width, height, nChannels, nSlices, nFrames) as a 5 element int array.
# dims = impThresh.getDimensions()

# Duplicate the image and crop to only keep the first channel
# run(ImagePlus imp, int firstC, int lastC, int firstZ, int lastZ, int firstT, int lastT)
# IJ.run("Duplicate...", "duplicate channels=1")
impThreshChnl1 = extractChannel(impThresh, 1, 1)

# impThreshChnl1.setTitle("impThreshChnl1")
# IJ.run("3D Watershed",
#        "seeds_threshold=0 image_threshold=0 image=impThreshChnl1 seeds=Markerimage radius=0")
# impWatershed = wm.getCurrentImage()

water        = Watershed3D(
    impThreshChnl1.getImageStack(), MarkerImage.getImageStack(), 0, 0)
impWatershed = water.getWatershedImage3D().getImagePlus()
# plus3.show()
impWatershed.setCalibration(imp2.getCalibration())
impWatershed.show()

# IJ.run("3D Manager Options", "volume integrated_density mean_grey_value feret distance_between_centers=10 distance_max_contact=1.80 drawing=Contour")
# IJ.run("3D Manager", "")
# Ext.Manager3D_AddImage()
# Ext.Manager3D_Measure()
# Ext.Manager3D_SaveResult("M", origInfo.directory + origInfo.fileName + "_ROI1_Results.csv")

# get current image, an image with labelled objects
# plus=WindowManager.getCurrentImage()
# wrap ImagePlus into 3D suite image format
img         = ImageInt.wrap(impWatershed)
# create a population of 3D objects
pop         = Objects3DPopulation(img)
nb          = pop.getNbObjects()
unit = imp2.getCalibration().getUnits()
# print(nb)
IHimp2      = ImageHandler.wrap(imp2)
volList     = []
meanIntList = []
feretList   = []


# loop over the objects
for i in range(0, nb):
    obj = pop.getObject(i)
    if(obj.touchBorders(img, filter_objects_touching_z)):
        continue
    volList.append(obj.getVolumeUnit())
    # print(obj.getIntegratedDensity())

    # Measure volume unit
    # print(obj.getMeasure(2))

    # Measure mean intensity
    meanIntList.append(obj.getPixMeanValue(IHimp2))

    # Measure feret
    feretList.append(obj.getFeret())


outCSV = outFullPath + ".csv"

# print(volList)
with open(outCSV, 'wb') as f:
    writer = csv.writer(f)
    writer.writerow(
        ["Volume (" + unit + " cube), mean Intensity, feret diameter (" + unit + ")"])
    writer.writerows(izip(volList, meanIntList, feretList))

# Add the result to the array
# cellCountImp.append(cellCountImpValue)
# Add the results to the results table
# rt.addValue('File name', basename)
# rt.addValue('Number of cells found', cellCountImpValue)
# IJ.log('' + str(cellCountImpValue) + ' particles found')
# currentFile = currentFile + 1
# Close the image
# imp.changes = False
# imp.close()


# Save the results
# rt.save(os.path.join(folder, 'Results.csv'))

IJ.log('Results were saved as ' + outCSV)
IJ.log('###########################')
IJ.log('Script done')
IJ.log('###########################')
