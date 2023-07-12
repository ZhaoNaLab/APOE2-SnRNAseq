import cv2
import numpy as np

def create_blurred_regions(I, sigma):
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
