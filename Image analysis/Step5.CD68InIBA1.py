import os
import cv2
import numpy as np
import matplotlib.pyplot as plt
import sys
import concurrent.futures
from scipy import ndimage
from scipy.ndimage import median_filter
import pandas as pd

# Add the path of the image_processing_util package to the system path
sys.path.append('.../image_processing_util')

# Import the ImgProcUtils module from the image_processing_util package
# this custom module contains 
from ImgProcUtils import *

# Define the directory to search for .tiff files
Dir = ".../IFAnalysis"

# Threshold values
threshold_red = 10486
threshold_green = 17497  # Set your threshold for green channel

# A function to process a single image
def process_image(file_path, threshold_red, threshold_green, depesckle_size=5):
    # Load image
    image = cv2.imread(file_path, cv2.IMREAD_UNCHANGED)

    # Threshold the red channel and despeckle
    img_red = image[:, :, 2]
    img_red_thres = np.where(img_red > threshold_red, np.iinfo(img_red.dtype).max, 0)
    img_red_thres_despeckle = median_filter(img_red_thres, size=depesckle_size)

    # Fill holes
    img_red_thres_despeckle_filled = ndimage.binary_fill_holes(img_red_thres_despeckle).astype(int)
    img_red_thres_despeckle_filled = img_red_thres_despeckle_filled * np.iinfo(img_red.dtype).max

    # Threshold the green channel
    img_green = image[:, :, 1]
    img_green_thres = np.where(img_green > threshold_green, img_green, 0)

    # Apply the red threshold mask to the green channel
    img_green_masked = np.where(img_red_thres_despeckle_filled, img_green_thres, 0)

    # Count the number of non-zero pixels and calculate the sum of pixel intensities
    non_zero_pixels = np.count_nonzero(img_green_masked)
    total_intensity = np.sum(img_green_masked)

    return non_zero_pixels, total_intensity

# Collect all .tiff file paths
tiff_files = []
for root, dirs, files in os.walk(Dir):
    for file in files:
        if file.endswith(".tiff") and "_thrBin" not in file:
            tiff_files.append(os.path.join(root, file))

# Use a ThreadPoolExecutor to process the images in parallel
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_image, tiff_files, [threshold_red]*len(tiff_files), [threshold_green]*len(tiff_files)))

# Extract Sample ID and Image name from fileList
sample_ids = [file.split("/")[-2].split(",")[0] for file in tiff_files]
image_names = [file.split("/")[-1].split(".")[0] for file in tiff_files]

# Create DataFrame
df = pd.DataFrame(results, columns=['PixelNumbers', 'TotalIntensities'])
df.insert(0, "ImageName", image_names)
df.insert(0, "SampleID", sample_ids)

# Save to CSV
df.to_csv(".../Results/CD68InIBA1.csv", index=False)
