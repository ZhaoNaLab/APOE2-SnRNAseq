import numpy as np

def centerImage(image):
    """
    Given an input image, put it in the center of a new image of given
    width and height. Assume the border is all 0's. Assume image is 8-bit.
    This functions expand the border of the image to the new image size.
    """

    # Get original image dimensions
    height, width = image.shape[:2]
    
    # Calculate new dimensions, 20% larger than the original ones
    new_height, new_width = int(height * 1.2), int(width * 1.2)

    # Create new image, filled with zeros (black)
    new_image = np.zeros((new_height, new_width, 3), dtype=np.uint8)

    # Calculate coordinates of top left corner of original image inside new image
    x = (new_width - width) // 2
    y = (new_height - height) // 2

    # Place the original image into the center of new image
    new_image[y:y+height, x:x+width] = image

    return new_image
