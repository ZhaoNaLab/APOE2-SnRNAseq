import numpy as np
import matplotlib.pyplot as plt

def channelExtremes(I_list, lower_threshold):
    """
    I_list : list of images, each element corresponding an M x N x 3 numpy array (RGB image of size M x N)
    lower_threshold : a vector of smallest pixel intensity levels that are considered "extreme"
    
    This function outputs a plot of the number of pixels that surpass lower_threshold, per image, for each channel.
    It returns a table (numpy array) of number of such pixels.
    """
    
    # Stack the image list into a 4D array of shape (M, N, number of images, 3)
    I_list_mat = np.stack(I_list, axis=2)
    
    # Broadcasting to apply thresholds across channels
    # Resulting extreme_counts has shape (number of images, 3)
    extreme_counts = np.sum(I_list_mat >= np.array(lower_threshold)[None, None, None, :], axis=(0, 1))
    
    # Plotting results
    plt.figure(figsize=(10, 15))

    channels = ['Red', 'Green', 'Blue']
    for i, channel in enumerate(channels):
        plt.subplot(3, 1, i+1)
        plt.plot(extreme_counts[:, i], 'o')
        plt.title(f'{channel} channel, intensity >= {lower_threshold[i]}')

    plt.tight_layout()
    plt.show()
    
    return extreme_counts

# Example usage
# Assuming you have a list of images named image_list and a lower_threshold list [r_threshold, g_threshold, b_threshold]
# channelExtremes(image_list, [r_threshold, g_threshold, b_threshold])
