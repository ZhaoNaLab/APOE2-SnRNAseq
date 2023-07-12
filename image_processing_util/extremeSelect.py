import cv2
import numpy as np
from skimage import measure, morphology
from scipy import ndimage as ndi

def extreme_select(I, upper_threshold, lower_threshold, struct_elem_size, min_area, max_area):
    # Define structuring element for erosion
    struct_elem = cv2.getStructuringElement(cv2.MORPH_RECT, (struct_elem_size, struct_elem_size))

    # Erode the image
    I_eroded = cv2.erode(I, struct_elem)

    # Apply thresholding and area filtering
    I_eroded_filt = ((I_eroded >= upper_threshold) * 1).astype(np.uint8)
    labeled_image, _ = measure.label(I_eroded_filt, return_num=True)
    props = measure.regionprops(labeled_image)
    for region in props:
        if region.area < min_area or region.area > max_area:
            I_eroded_filt[labeled_image == region.label] = 0

    # Find the coordinates of extreme points
    coords = np.column_stack(np.where(I_eroded_filt > 0))

    # Apply thresholding and select points connected to the extreme points
    mask = I > lower_threshold
    labeled_mask, _ = ndi.label(mask)
    I_extreme_filt_bin = np.zeros_like(I, dtype=bool)
    for coord in coords:
        I_extreme_filt_bin |= ndi.binary_fill_holes(labeled_mask == labeled_mask[coord[0], coord[1]])

    return I_extreme_filt_bin