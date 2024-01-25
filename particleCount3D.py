# ─── SCRIPT PARAMETERS ──────────────────────────────────────────────────────────

#@ Double(label="Minimum size of spots", value=0.05) min_volume_spots
#@ OpService ops
#@ RoiManager rm

# ─── IMPORTS ────────────────────────────────────────────────────────────────────

from net.imglib2.img.display.imagej import ImageJFunctions

import os
import sys
import csv
import glob
from itertools import izip

from ij import IJ, ImagePlus, ImageStack, WindowManager as wm
from ij.plugin import Duplicator
from ij.plugin.frame import RoiManager
from ij.gui import WaitForUserDialog
from ij.measure import ResultsTable
from ij.process import ImageConverter

# Bioformats imports
from loci.plugins import BF
from loci.plugins.in import ImporterOptions

# 3DSuite imports
from mcib3d.geom import Objects3DPopulation
from mcib3d.image3d import ImageInt, ImageHandler
from mcib3d.image3d.segment import Segment3DImage



# ─── VARIABLES ──────────────────────────────────────────────────────────────────

# H-watershed settings
hMin_cells                   = 150
segmentation_threshold_cells = 400
peak_flooding_cells          = 100
output_mask_cells            = True
allow_split_cells            = True

# min_volume_spots = 0.05

filter_objects_touching_z = False

# ─── FUNCTIONS ──────────────────────────────────────────────────────────────────

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


# ─── MAIN CODE ──────────────────────────────────────────────────────────────────

# Get image and some info
imp      = wm.getCurrentImage()
unit = imp.getCalibration().getUnits()
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
out_full_pathTest = os.path.join(dir, testName)
count = 1

# Save a new file for each ROI
while checkForFiles(out_full_pathTest):
    # print(out_full_pathTest)
    outFileName     = fileName + "_" + str(count)
    testName        = outFileName + ".*"
    out_full_pathTest = os.path.join(dir, testName)
    count           = count + 1

out_full_path = os.path.join(dir, outFileName)
# print(outFileName)
# sys.exit()

# rm.reset()
# print(rm.getCount())

# Awaiting ROI or quits after 5 tries
count = 0
while((imp.getRoi() is None) and (rm.getCount() == 0 ) ):
    WaitForUserDialog(
        "Draw the region of interest and press OK").show()
    if count == 5:
        sys.exit("Too many clicks without ROI")
    else:
        count = count + 1

if (rm.getCount() == 1):
	rm.select(imp, 0)
	region_roi = rm.getRoi(0)
else:
	rm.reset()
	region_roi = imp.getRoi()
	rm.addRoi(region_roi)
	rm.select(imp, 0)

rm.runCommand("Save", out_full_path + ".zip")


imp_cropped = Duplicator().run(imp, 1, 1, 1, imp.getNSlices(), 1, 1)
IJ.run(imp_cropped, "Median...", "radius=2 stack")
# imp_cropped.show()
IH_cropped  = ImageHandler.wrap(imp_cropped)

imp_watershed = ops.run("H_Watershed", imp_cropped, hMin_cells,
                segmentation_threshold_cells, peak_flooding_cells,
                output_mask_cells, allow_split_cells)
imp_watershed.setTitle("Watershed")
imp_watershed.setCalibration(imp.getCalibration())

imp_spots_segmented    = Segment3DImage(imp_watershed, 1, 65535)
imp_spots_segmented.segment()
stack_spots_segmented  = imp_spots_segmented.getLabelledObjectsStack()
spots_label            = ImagePlus("Spots_label", stack_spots_segmented)
spots_label.setCalibration(imp.getCalibration())
img_spots              = ImageInt.wrap(spots_label)
pop_spots              = Objects3DPopulation(img_spots)
nb_spots               = pop_spots.getNbObjects()
nb_spots_before        = nb_spots
spots_obj_to_remove    = []
spots_names            = []
spots_volumes          = []
spots_ferets           = []
spots_mean_intensities = []

for i in range(0, nb_spots):
    obj = pop_spots.getObject(i)

    if obj.getVolumeUnit() < min_volume_spots:
        spots_obj_to_remove.append(obj)
        continue
    elif obj.touchBorders(img_spots, filter_objects_touching_z):
        spots_obj_to_remove.append(obj)
        continue
    else:
        spots_names.append(obj.getName())
        spots_volumes.append(obj.getVolumeUnit())
        spots_ferets.append(obj.getFeret())
        spots_mean_intensities.append(obj.getPixMeanValue(IH_cropped))
        obj.translate(region_roi.getBounds().x, region_roi.getBounds().y, 0)

for obj in spots_obj_to_remove:
    pop_spots.removeObject(obj)

pop_spots.saveObjects(out_full_path + "_3DROIs.zip")

out_CSV = out_full_path + ".csv"

# print(volList)
with open(out_CSV, 'wb') as f:
    writer = csv.writer(f)
    writer.writerow(
        ["Object ID, Volume (" + unit + " cube), mean Intensity, " +
        "feret diameter (" + unit + ")"])
    writer.writerows(izip(spots_names, spots_volumes, spots_mean_intensities,
        spots_ferets))


IJ.log('Results were saved as ' + out_CSV)
IJ.log('###########################')
IJ.log('Script done')
IJ.log('###########################')
