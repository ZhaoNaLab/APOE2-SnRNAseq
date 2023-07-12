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
        if file.endswith(".tiff"):
            tiff_files.append(os.path.join(root, file))

# Use a ThreadPoolExecutor to load the images in parallel
with concurrent.futures.ThreadPoolExecutor() as executor:
    image_arrays = list(executor.map(load_image, tiff_files))

# Print the list of image arrays
len(image_arrays)
# print(image_arrays[0])

## Let's take a look at the data distribution of all the images for each channel
InDir = ".../IFAnalysis"
OutDir = ".../IFAnalysis/Results"

# SampleChannelIntensity(input_dir = InDir, output_dir = OutDir, N = 20, file_name = "IntensityDist1")

## Let's do this for all images
channelDensity(I_list=image_arrays, output_dir=OutDir, file_name="IntensityDistAll")
channelIntensityMedian(I_list=image_arrays, output_dir=OutDir, file_name="IntensityMedian")
