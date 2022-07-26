[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6801606.svg)](https://doi.org/10.5281/zenodo.6801606)

# particle-count-3d

## Description

This script allows to detect 3D particles inside a cell selected by the user and then gives the volume, the mean intensity and the feret diameter of these particles.

## Requirements

The script is made in Python and requires [Fiji](http://fiji.sc/) to be ran. On top of this, multiple update sites need to be activated following [this guide](https://imagej.net/update-sites/#following-an-update-site): 
* 3DImageSuite
* MorpholibJ

Once activated, just drag and drop the script in the main Fiji window and click on the RUN button.

As Fiji is operating system independant, this can be run on Windows, Mac and Linux. No specific hardware is necessary and script should take a couple of minutes to finish.

## Run script

### Input

The script requires an image opened in Fiji on which to run. If no ROI are available, will prompt the user to create one. 

### Runtime 

Once selected, a white top hat filter will be applied to homogenize the background and TrackMate is then used to create a marker image containing 3D seeds for the particles.
An iterative threshold is applied on the original image to segment the particles and the 3D seeds are used to correct for wrongly separated objects.

The script will then loop through all the 3D objects, filtering out the ones touching the borders in X and Y, and measure volume, mean intensity and the feret diameter, saving the results in a CSV.

### Output

The script will save the drawn ROI for reproducibility (will increment number for saving if already ran on the same image) and save the measures in a CSV with the same name as the ROI.