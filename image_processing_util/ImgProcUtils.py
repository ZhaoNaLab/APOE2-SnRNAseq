import cv2
import numpy as np
import seaborn as sns
import random
import numpy as np
import matplotlib.pyplot as plt
import skfmm
from skimage import filters
import os
import pathlib
from matplotlib.backends.backend_pdf import PdfPages
from concurrent.futures import ProcessPoolExecutor
import tifffile as tiff
from scipy import ndimage as ndi
from scipy.ndimage import median_filter, distance_transform_edt
import pandas as pd
from skimage import measure
import numpy as np


# Function to load image and convert to numpy array
def load_image(file):
    image = cv2.imread(str(file), cv2.IMREAD_UNCHANGED)
    return np.array(image)

def centerImage(image):
    """
    Given an input image, put it in the center of a new image of given
    width and height. Assume the border is all 0's. Assume image is 8-bit.
    This functions expand the border of the image to the new image size.
    """

    # Get original image dimensions
    height, width = image.shape[:2]
    
    # Calculate new dimensions, 20% larger than the original ones
    new_height, new_width = int(height * 1.2), int(width * 1.2)

    # Create new image, filled with zeros (black)
    new_image = np.zeros((new_height, new_width, 3), dtype=np.uint8)

    # Calculate coordinates of top left corner of original image inside new image
    x = (new_width - width) // 2
    y = (new_height - height) // 2

    # Place the original image into the center of new image
    new_image[y:y+height, x:x+width] = image

    return new_image

def extremeSelect(I, upper_threshold, lower_threshold, struct_elem_width, struct_elem_height, min_area, max_area):
    # Define structuring element for erosion
    struct_elem = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (struct_elem_width, struct_elem_height))

    # Erode the image
    I_eroded = cv2.erode(I, struct_elem)

    # Apply thresholding and area filtering
    I_eroded_filt = ((I_eroded >= upper_threshold) * 1).astype(np.uint8)
    labeled_image, _ = measure.label(I_eroded_filt, return_num=True)
    props = measure.regionprops(labeled_image)
    
    # Vectorized area filtering
    region_labels = [region.label for region in props if region.area < min_area or region.area > max_area]
    for label in region_labels:
        I_eroded_filt[labeled_image == label] = 0

    # Find the coordinates of extreme points
    coords = np.column_stack(np.where(I_eroded_filt > 0))

    # Apply thresholding and select points connected to the extreme points
    mask = I > lower_threshold
    labeled_mask, _ = ndi.label(mask)
    I_extreme_filt_bin = np.zeros_like(I, dtype=bool)
    for coord in coords:
        I_extreme_filt_bin |= ndi.binary_fill_holes(labeled_mask == labeled_mask[coord[0], coord[1]])

    return I_extreme_filt_bin


def generateColors(num_colors, seed=None):
    """Generates a set of discrete colors."""
    if seed is not None:
        random.seed(seed)
    colors = ['#%06x' % random.randint(0, 0xFFFFFF) for _ in range(num_colors)]
    return colors

# Function to calculate and visualize channel density
def channelDensity(I_list, output_dir, file_name):
    # Convert the list to a numpy array and reshape it for easier access
    I_list_arr = np.array(I_list).reshape((-1, *I_list[0].shape))
    
    # Get a list of colors to use
    colors = generateColors(num_colors=len(I_list), seed=12535)

    fig, axs = plt.subplots(3, 1, figsize=(10, 8))

    for idx, channel_name in enumerate(['R', 'G', 'B']):
        for i, color in enumerate(colors):
            sns.kdeplot(I_list_arr[i, :, :, idx].ravel(), color=color, alpha=0.7, label=f'Image {i+1}', ax=axs[idx])
        axs[idx].set_title(channel_name)
        axs[idx].legend()

    plt.tight_layout()
    
    # Save the figure as a PDF
    with PdfPages(os.path.join(output_dir, file_name + '.pdf')) as pdf:
        pdf.savefig(fig)
    plt.close(fig)

# Function to sample directory and files
def SampleChannelIntensity(input_dir, output_dir, N, file_name):
    # Define the directory to search for subdirectories
    directory = pathlib.Path(input_dir)

    # Recursively search the directory and subdirectories for image files
    subdirs = [f for f in directory.glob("**/") if any(f.glob("*.tiff"))]

    # Randomly select N subdirectories
    random_dirs = random.sample(subdirs, N)

    # From these N directories, randomly select 1 tiff file each
    random_files = [random.choice(list(d.glob("*.tiff"))) for d in random_dirs]

    # Load all the images using parallel processing
    with ProcessPoolExecutor() as executor:
        image_arrays = list(executor.map(load_image, random_files))

    # Calculate and visualize channel density
    channelDensity(image_arrays, output_dir, file_name)

import pandas as pd

def PositiveCount(I_list, lower_threshold, Filelist, OutsputDir, FileName):
    """
    I_list : list of images, each element corresponding an M x N x 3 numpy array (RGB image of size M x N)
    lower_threshold : a vector of smallest pixel intensity levels that are considered "extreme"
    
    This function returns a dataframe with the number of pixels that surpass lower_threshold, per image, for each channel.
    """

    # Extract Sample ID and Image name from fileList
    sample_ids = [file.split("/")[-2].split(",")[0] for file in Filelist]
    image_names = [file.split("/")[-1].split(".")[0] for file in Filelist]

    # Stack the image list into a 4D array of shape (M, N, number of images, 3)
    I_list_mat = np.stack(I_list, axis=2)
    
    # Broadcasting to apply thresholds across channels
    # Resulting extreme_counts has shape (number of images, 3)
    positive_counts = np.sum(I_list_mat >= np.array(lower_threshold)[None, None, None, :], axis=(0, 1))
    
    # Create a DataFrame
    df = pd.DataFrame(positive_counts, columns=['BlueCount', 'GreenCount', 'RedCount'])
    df.insert(0, "ImageName", image_names)
    df.insert(0, "SampleID", sample_ids)

    # Save the DataFrame as a CSV
    df.to_csv(os.path.join(OutsputDir, FileName + '.csv'), index=False)
    
    return df


# Example usage
# Assuming you have a list of images named image_list and a lower_threshold list [r_threshold, g_threshold, b_threshold]
# channelExtremes(image_list, [r_threshold, g_threshold, b_threshold])

def PositiveIntensity(I_list, lower_threshold, Filelist, OutsputDir, FileName):
    """
    I_list : list of images, each element corresponding an M x N x 3 numpy array (RGB image of size M x N)
    lower_threshold : a vector of smallest pixel intensity levels that are considered "extreme"
    
    This function outputs a plot of the total intensity of pixels that surpass lower_threshold, per image, for each channel.
    It returns a table (numpy array) of total intensity of such pixels.
    """

    # Extract Sample ID and Image name from fileList
    sample_ids = [file.split("/")[-2].split(",")[0] for file in Filelist]
    image_names = [file.split("/")[-1].split(".")[0] for file in Filelist]

    # Stack the image list into a 4D array of shape (M, N, number of images, 3)
    I_list_mat = np.stack(I_list, axis=2)
    
    # Broadcasting to apply thresholds across channels
    # Resulting extreme_counts has shape (number of images, 3)
    positive_intensity = np.sum(np.where(I_list_mat >= np.array(lower_threshold)[None, None, None, :], I_list_mat, 0), axis=(0, 1))
    
    # Create a DataFrame
    df = pd.DataFrame(positive_intensity, columns=['BlueIntensity', 'GreenIntensity', 'RedIntensity'])
    df.insert(0, "ImageName", image_names)
    df.insert(0, "SampleID", sample_ids)

    # Save the DataFrame as a CSV
    df.to_csv(os.path.join(OutsputDir, FileName + '.csv'), index=False)
    
    return df

# A function to quantify the number and intensity of pixels that surpass a threshold
def GreenInRed(file_path, threshold):
    """
    Extracts the green channel from an image where the red channel is above a certain threshold.

    Args:
        file_path (str): The path to the image file.
        threshold (int): The threshold value for the red channel.

    Returns:
        Tuple[int, int]: A tuple containing the number of non-zero pixels and the sum of pixel intensities
            in the green channel of the masked image.
    """
    # Load image
    image = cv2.imread(file_path, cv2.IMREAD_UNCHANGED)

    # Threshold the red channel and despeckle
    img_red = image[:, :, 2]
    img_red_thres = np.where(img_red > threshold, np.iinfo(img_red.dtype).max, 0)
    img_red_thres_despeckle = median_filter(img_red_thres, size=5)

    # Fill holes
    img_red_thres_despeckle_filled = ndi.binary_fill_holes(img_red_thres_despeckle).astype(int)
    img_red_thres_despeckle_filled = img_red_thres_despeckle_filled * np.iinfo(img_red.dtype).max

    # Apply the mask
    img_green_masked = np.where(img_red_thres_despeckle_filled, image[:, :, 1], 0)

    # Count the number of non-zero pixels and calculate the sum of pixel intensities
    non_zero_pixels = np.count_nonzero(img_green_masked)
    total_intensity = np.sum(img_green_masked)

    return non_zero_pixels, total_intensity

def channelIntensityMean(I_list, output_dir, file_name):
    # I_list : list of images, assumed to be a list of numpy arrays of some length,
    # each element corresponding an M X N x 3 matrix (RGB image of size M x N)
    # output : a plot of the mean value for each of the 3 channels, grouped by
    # image. If we got the wrong channel, we'll see it in this plot.

    I_list_mean = [np.mean(x.reshape(-1, 3), axis=0) for x in I_list]
    I_list_mat_mean = np.array(I_list_mean).T

    image_ids = range(len(I_list))

    plt.subplot(3,1,1)
    plt.bar(image_ids, np.log1p(I_list_mat_mean[0,:]), tick_label=image_ids)
    plt.title('R')
    plt.xlabel('Image ID')

    plt.subplot(3,1,2)
    plt.bar(image_ids, np.log1p(I_list_mat_mean[1,:]), tick_label=image_ids)
    plt.title('G')
    plt.xlabel('Image ID')

    plt.subplot(3,1,3)
    plt.bar(image_ids, np.log1p(I_list_mat_mean[2,:]), tick_label=image_ids)
    plt.title('B')
    plt.xlabel('Image ID')

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{file_name}.pdf")
    plt.show()


def channelIntensityMedian(I_list, output_dir, file_name):
    # I_list : list of images, assumed to be a list of numpy arrays of some length,
    # each element corresponding an M X N x 3 matrix (RGB image of size M x N)
    # output : a plot of the median value for each of the 3 channels, grouped by
    # image. If we got the wrong channel, we'll see it in this plot.

    I_list_median = [np.median(x.reshape(-1, 3), axis=0) for x in I_list]
    I_list_mat_median = np.array(I_list_median).T

    image_ids = range(len(I_list))

    plt.subplot(3,1,1)
    plt.bar(image_ids, np.log1p(I_list_mat_median[0,:]), tick_label=image_ids)
    plt.title('R')
    plt.xlabel('Image ID')

    plt.subplot(3,1,2)
    plt.bar(image_ids, np.log1p(I_list_mat_median[1,:]), tick_label=image_ids)
    plt.title('G')
    plt.xlabel('Image ID')

    plt.subplot(3,1,3)
    plt.bar(image_ids, np.log1p(I_list_mat_median[2,:]), tick_label=image_ids)
    plt.title('B')
    plt.xlabel('Image ID')

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{file_name}.pdf")
    plt.show()

def Darkthreshold(I, D):
    """
    Given the dark region D, compute threshold of I using
    the max pixel intensity of I within D.
    Use the threshold to return a thresholded image.
    """
    T = np.max(I[D==1])
    _, I_thr = cv2.threshold(I, T, np.iinfo(I.dtype).max, cv2.THRESH_BINARY)
    return I_thr, T

def correctFixedDarkPatternNoise(I):
    """
    Remove dark-pattern noise from 3-channel images.
    """
    I = cv2.normalize(I, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    _, I_1 = cv2.threshold(I[:,:,0], 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)
    _, I_2 = cv2.threshold(I[:,:,1], 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)
    _, I_3 = cv2.threshold(I[:,:,2], 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)

    fixed_dark_pattern = np.logical_not(I_1)
    fixed_dark_pattern[np.logical_or(I_2, I_3)] = 0

    I_thr = np.zeros_like(I)
    T = np.zeros((3,))
    I_thr[:,:,0], T[0] = Darkthreshold(I[:,:,0], fixed_dark_pattern)
    I_thr[:,:,1], T[1] = Darkthreshold(I[:,:,1], fixed_dark_pattern)
    I_thr[:,:,2], T[2] = Darkthreshold(I[:,:,2], fixed_dark_pattern)
    
    return I_thr, T

def createBlurredRegions(I, sigma):
    """
    Select for regions of image that appear to have some connected features.
    I : input image (uint8). Recommended to remove uninformative, noisy stuff from I before using createBlurredRegions.
    sigma : degree of spatial averaging via 2-D symmetric Gaussian filter.
    """
    if sigma == 0:  # no blurring desired.
        _, I_blurred_bin = cv2.threshold(I, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)
    else:
        # I = np.where(I > 0, 255, 0).astype(np.uint8)  # Assume all non-zero pixels of I informative.
        I_blurred = cv2.GaussianBlur(I, (0, 0), sigma)  # Apply Gaussian filter
        _, I_blurred_bin = cv2.threshold(I_blurred, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)  # Binarize

    return I_blurred_bin


def createDistContourLevels(image, levels):
    """
    Create distance contours given binary image.
    Given levels, categorize contours into regions 0, 1, ..., length(levels) - 1
    Example:
    levels = [10, 46, 82]
    level 0 : 0 - 10
    level 1 : 10 - 46
    level 2 : 46 - 82
    level 3 : > 82
    """
    # Calculate the Euclidean distance transform of the image
    dist_transform = distance_transform_edt(image)
    
    # Initialize the contour image with zeros
    contour_image = np.zeros_like(dist_transform)
    
    # Loop through the levels in reverse order and update the contour image
    for index, level in enumerate(reversed(levels)):
        current_level = len(levels) - index
        if current_level == len(levels):
            contour_image[dist_transform > level] = current_level
        contour_image[dist_transform <= level] = current_level - 1

    return contour_image

def createDistGeodesicContourLevels(I, mask, levels):
    # Create distance contours given binary image I
    # Given levels, categorize contours into regions 0, 1, ..., len(levels) - 1
    # Example:
    # levels = [10, 46, 82]
    # level 0 : 0 - 10
    # level 1 : 10 - 46
    # level 2 : 46 - 82
    # level 3 : > 82

    # Convert mask to boolean where nonzero is True

    # Compute distance transform
    Id = skfmm.distance(mask)

    I_contours = np.zeros(Id.shape)
    L = len(levels)
    for l in reversed(levels):
        if L == len(levels):
            I_contours[Id > l] = L
        if L > 0:
            I_contours[Id <= l] = L - 1
        L = L - 1
    return I_contours


def createImageMontage(image_list, nrow, ncol, wspace=0.1, hspace=0.1):
    # Assume image_list is a list of image paths
    images = [plt.imread(img) for img in image_list]

    # Get the size of the first image
    size = images[0].shape

    # Check if all images are the same size
    for img in images:
        if img.shape != size:
            raise ValueError('Not all images have the same size.')

    # Compute number of images
    n = len(images)

    # Create a montage
    fig, axs = plt.subplots(nrow, ncol)

    for i, ax in enumerate(axs.flat):
        if i < n:
            # Display image if index is less than number of images
            ax.imshow(images[i])
        else:
            # Otherwise hide axes
            ax.axis('off')

        # Remove axis ticks
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    # Remove whitespace between subplots
    plt.tight_layout()
    plt.subplots_adjust(wspace=wspace, hspace=hspace)
    plt.show()

def estimateThresholds(images):
    """
    Given a list of images, estimate thresholds for each image using the
    correct_fixed_dark_pattern_noise function, and return the average threshold
    across all images. The thresholds and intermediate results are also returned.
    """
    thresholded_images, thresholds = map(correct_fixed_dark_pattern_noise, images)
    intermediate_results = thresholded_images
    avg_threshold = np.mean(thresholds)
    return avg_threshold, thresholds, intermediate_results


def removeObstacle(I_list):
    I_obstacle = I_list[0] > 0
    for i in range(1, len(I_list)):
        I_obstacle = I_obstacle * (I_list[i] > 0)
    for i in range(len(I_list)):
        I_list[i][I_obstacle > 0] = 0
    return I_list

def sharedThreshold(I_list, fileList, OutDir):
    # Concatenate images along the horizontal axis
    I_list_cat = np.concatenate(I_list, axis=1)

    # Separate RGB channels
    I_list_cat_1 = I_list_cat[:,:,0]
    I_list_cat_2 = I_list_cat[:,:,1]
    I_list_cat_3 = I_list_cat[:,:,2]

    # Determine Otsu's threshold for each channel
    T1 = filters.threshold_otsu(I_list_cat_1)
    T2 = filters.threshold_otsu(I_list_cat_2)
    T3 = filters.threshold_otsu(I_list_cat_3)
    
    # Build the result string
    result = f"The threshold for blue channel is {T1};\n"
    result += f"The threshold for green channel is {T2};\n"
    result += f"The threshold for red channel is {T3};"

    I_thr_list = []
    int_max = np.iinfo(I_list[0].dtype).max

    # Apply the threshold to each image
    for i, I in enumerate(I_list):
        I_thr = np.zeros_like(I)
        I_thr[:,:,0] = np.where(I[:,:,0] > T1, int_max, 0) # Blue channel
        I_thr[:,:,1] = np.where(I[:,:,1] > T2, int_max, 0) # Green channel
        I_thr[:,:,2] = np.where(I[:,:,2] > T3, int_max, 0) # Red channel
        I_thr_list.append(I_thr)

        # Save the thresholded image with a different name
        file_name = os.path.splitext(os.path.basename(fileList[i]))[0]
        new_file_name = file_name + '_thrBin.tiff'
        new_file_path = os.path.join(os.path.dirname(fileList[i]), new_file_name)
        # cv2.imwrite(new_file_path, I_thr)
        tiff.imwrite(new_file_path, I_thr)

    # Save the result string as a text file
    out_dir = OutDir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_file_path = os.path.join(out_dir, "SharedThreshold.txt")
    with open(out_file_path, "w") as f:
        f.write(result)

    return result

def tanhThreshold(I, m, b):
    """
    Apply a tanh thresholding operation to an input image.

    Parameters
    ----------
    I : numpy.ndarray
        The input image.

    m : float
        The gain, which approximately represents the slope of the linear region at b.

    b : float
        The intercept, which is subtracted from the rescaled image.

    Returns
    -------
    numpy.ndarray
        The thresholded image.
    """
    I = I.astype(float)
    I_rescaled = (I - I.min()) / (I.max() - I.min())
    I_thr = np.tanh(m * (I_rescaled - b)) * 0.5 + 0.5
    I_thr_rescaled = I_thr * (I.max() - I.min()) + I.min()
    return I_thr_rescaled.astype(I.dtype)


def threshold(I, T_min, T_max):
    I_thr = np.interp(I, (T_min, T_max), (0, np.iinfo(I.dtype).max)).astype(I.dtype)
    I_thr[I < T_min] = 0
    return I_thr
