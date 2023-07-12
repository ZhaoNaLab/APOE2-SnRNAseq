import numpy as np
import matplotlib.pyplot as plt

def channelIntensityMean(I_list):
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
    plt.show()