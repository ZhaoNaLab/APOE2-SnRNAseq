import numpy as np

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

