import os
import cv2
import numpy as np
import matplotlib.pyplot as plt
import sys
import concurrent.futures

# Add the path of the image_processing_util package to the system path
sys.path.append('.../image_processing_util')

# Import the ImgProcUtils module from the image_processing_util package
# this custom module contains 
from ImgProcUtils import *


####### To load all images for the identification of a common threshold for each channel ####
# Define the directory to search for .tiff files
Dir = ".../IFAnalysis"

# A function to load a single image
def load_image(file_path):
    image = cv2.imread(file_path, cv2.IMREAD_UNCHANGED)
    return np.array(image)

# Collect all .tiff file paths
tiff_files = []
for root, dirs, files in os.walk(Dir):
    for file in files:
        if file.endswith(".tiff") and "_thrBin" not in file:
            tiff_files.append(os.path.join(root, file))

# Use a ThreadPoolExecutor to load the images in parallel
with concurrent.futures.ThreadPoolExecutor() as executor:
    image_arrays = list(executor.map(load_image, tiff_files))

# Print the list of image arrays
len(image_arrays)

# Count the positive signals
OutDir = ".../Results"
PositiveCount(I_list=image_arrays, lower_threshold=[18849, 17497, 10486], Filelist=tiff_files, OutsputDir=OutDir, FileName="PositiveCount")
PositiveIntensity(I_list=image_arrays, lower_threshold=[18849, 17497, 10486], Filelist=tiff_files, OutsputDir=OutDir, FileName="PositiveIntensity")




