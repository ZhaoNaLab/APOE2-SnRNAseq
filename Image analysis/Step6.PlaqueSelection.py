import os
import cv2
import numpy as np
import pandas as pd
import concurrent.futures
from skimage import measure
import logging

logging.basicConfig(filename='image_processing.log', level=logging.INFO, 
                    format='%(asctime)s:%(levelname)s:%(message)s')

def collect_tiff_files(directory):
    """
    Collects all .tiff file paths in the directory.

    Args:
    directory (str): The directory to search for .tiff files.

    Returns:
    list: A list of file paths to .tiff files in the directory.
    """
    tiff_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".tiff") and "_thrBin" not in file:
                tiff_files.append(os.path.join(root, file))
    return tiff_files

def load_image(file_path):
    """
    Loads a single image.

    Args:
    file_path (str): The file path to the image.

    Returns:
    numpy.ndarray: The loaded image as a numpy array.
    """
    try:
        image = cv2.imread(file_path, cv2.IMREAD_UNCHANGED)
        return np.array(image)
    except Exception as e:
        logging.error(f"Error loading image {file_path}: {e}")
        return None

def process_image(image, threshold, struct_elem_erod=3, struct_elem_dilate=3, dilate_N=1, erod_N=1):
    """
    Processes an image, applying thresholding, erosion, hole filling, and median filtering.
    Returns the median filtered image and its region properties.

    Args:
    image (numpy.ndarray): The image to process.
    threshold (int): The threshold value for the blue channel.
    struct_elem_erod (int): The size of the structuring element for erosion.
    struct_elem_dilate (int): The size of the structuring element for dilation.
    dilate_N (int): The number of times to dilate the image.
    erod_N (int): The number of times to erode the image.

    Returns:
    tuple: A tuple containing the median filtered image and its region properties.
    """
    try:
        blue_channel = image[:, :, 0]
        _, thresh = cv2.threshold(blue_channel, threshold, np.iinfo(blue_channel.dtype).max, cv2.THRESH_BINARY)
        struct_elem_erode = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (struct_elem_erod, struct_elem_erod))
        struct_elem_dilate = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (struct_elem_dilate, struct_elem_dilate))
        dilated = cv2.dilate(thresh, struct_elem_dilate, iterations = dilate_N)
        eroded = cv2.erode(dilated, struct_elem_erode, iterations = erod_N)

        filled = cv2.morphologyEx(eroded, cv2.MORPH_CLOSE, struct_elem_erode)
        median_filtered = cv2.medianBlur(filled, 5)
        labels = measure.label(median_filtered, connectivity=2)
        regions = measure.regionprops(labels)
        return median_filtered, regions
    except Exception as e:
        logging.error(f"Error processing image: {e}")
        return None, None

def save_processed_image(image, file_path):
    """
    Saves the processed image to a .tiff file with "_processed" added to the end of the original image name.

    Args:
    image (numpy.ndarray): The processed image to save.
    file_path (str): The file path to the original image.
    """
    try:
        outDir = os.path.dirname(file_path)
        outName = os.path.splitext(os.path.basename(file_path))[0] + "_processedBlue.tiff"
        cv2.imwrite(os.path.join(outDir, outName), image)
    except Exception as e:
        logging.error(f"Error saving image {file_path}: {e}")

def main():
    """
    The main function that processes all .tiff files in a directory.
    """
    Dir = ".../IFAnalysis"
    threshold = 18849
    struct_elem_dilate = 3
    struct_elem_erod = 3
    dilate_N = 10
    erod_N = 3

    tiff_files = collect_tiff_files(Dir)
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        image_arrays = list(executor.map(load_image, tiff_files))
    
    for idx, image in enumerate(image_arrays):
        if image is not None:  # Only process if image was loaded successfully
            median_filtered, regions = process_image(image, threshold, struct_elem_dilate, struct_elem_erod, dilate_N, erod_N)
            if median_filtered is not None and regions is not None:  # Only proceed if image was processed successfully
                save_processed_image(median_filtered, tiff_files[idx])
    
if __name__ == "__main__":
    main()
