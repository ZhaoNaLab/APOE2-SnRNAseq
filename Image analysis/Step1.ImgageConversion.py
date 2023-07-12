import xml.etree.ElementTree as ET
from aicsimageio import AICSImage
import numpy as np
from tifffile import imwrite
import os
from glob import glob
import shutil

from aicspylibczi import CziFile

def CZIToRBG(image, channel):
    """
    Process a single color channel of an image.
    """
    data = image.get_image_data("ZYX", C=channel, S=0, T=0)
    data = np.squeeze(data)
    return data

def convert_image(input_path, output_path):
    try:
        # Load the image
        img = AICSImage(input_path)
        Czifile = CziFile(input_path)
    except Exception as e:
        print(f"Error loading image '{input_path}': {e}")
        return

    try:
        # Get the Excitation Wavelengths
        root = Czifile.meta
        wavelengths = []
        for channel in root.findall(".//Dimensions/Channels/Channel"):
            excitation_wavelength = channel.find("./ExcitationWavelength")
            if excitation_wavelength is not None:
                excitation_wavelength = float(excitation_wavelength.text)
            else:
                excitation_wavelength = None
            wavelengths.append(excitation_wavelength)

        # Order the channels based on the wavelengths
        wavelength_order = [543, 488.00000000000006, 405.00000000000006]
        sorted_channels = sorted(list(range(len(wavelengths))), key=lambda i: wavelength_order.index(wavelengths[i]))

        # Process each color channel based on the order
        channels = [CZIToRBG(img, i) for i in sorted_channels]

        # Merge the grayscale images into a three-channel image
        rgb_image = np.dstack(channels)
    except Exception as e:
        print(f"Error processing image '{input_path}': {e}")
        return

    try:
        # Save the image to tiff file
        imwrite(output_path, rgb_image)
    except Exception as e:
        print(f"Error saving image '{output_path}': {e}")
        return

# Define input and output directories
input_dir = "Path/to/czi files"
output_dir = ".../IFAnalysis"

# Find all .czi files in the input directory and its subdirectories
input_files = glob(os.path.join(input_dir, "**", "*.czi"), recursive=True)

for input_file in input_files:
    # Define the output file path
    relative_path = os.path.relpath(input_file, input_dir)
    output_file = os.path.splitext(relative_path)[0] + ".tiff"
    output_path = os.path.join(output_dir, output_file)

    # Create the output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Convert and save the image
    convert_image(input_file, output_path)

