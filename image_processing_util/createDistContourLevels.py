import numpy as np
from scipy.ndimage import distance_transform_edt

def create_dist_contour_levels(image, levels):
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
