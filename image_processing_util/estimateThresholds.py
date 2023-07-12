import numpy as np

def estimate_thresholds(images):
    """
    Given a list of images, estimate thresholds for each image using the
    correct_fixed_dark_pattern_noise function, and return the average threshold
    across all images. The thresholds and intermediate results are also returned.
    """
    thresholded_images, thresholds = map(correct_fixed_dark_pattern_noise, images)
    intermediate_results = thresholded_images
    avg_threshold = np.mean(thresholds)
    return avg_threshold, thresholds, intermediate_results