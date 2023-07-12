import cv2
import numpy as np
import matplotlib.pyplot as plt

def dark_threshold(I, D):
    """
    Given the dark region D, compute threshold of I using
    the max pixel intensity of I within D.
    Use the threshold to return a thresholded image.
    """
    T = np.max(I[D==1])
    _, I_thr = cv2.threshold(I, T, np.iinfo(I.dtype).max, cv2.THRESH_BINARY)
    return I_thr, T

def correct_fixed_dark_pattern_noise(I):
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
    I_thr[:,:,0], T[0] = dark_threshold(I[:,:,0], fixed_dark_pattern)
    I_thr[:,:,1], T[1] = dark_threshold(I[:,:,1], fixed_dark_pattern)
    I_thr[:,:,2], T[2] = dark_threshold(I[:,:,2], fixed_dark_pattern)
    
    return I_thr, T