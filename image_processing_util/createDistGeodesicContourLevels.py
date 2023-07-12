import cv2
import numpy as np
import skfmm

def create_dist_geodesic_contour_levels(I, mask, levels):
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
