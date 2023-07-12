from skimage import filters
import numpy as np

def shared_threshold(I_list):
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
    T = [T1, T2, T3]

    I_thr_list = []
    int_max = np.iinfo(I_list[0].dtype).max

    # Apply the threshold to each image
    for I in I_list:
        I_thr = np.zeros_like(I)
        I_thr[:,:,0] = np.where(I[:,:,0] > T1, int_max, 0)
        I_thr[:,:,1] = np.where(I[:,:,1] > T2, int_max, 0)
        I_thr[:,:,2] = np.where(I[:,:,2] > T3, int_max, 0)
        I_thr_list.append(I_thr)
        
    return T, I_thr_list

