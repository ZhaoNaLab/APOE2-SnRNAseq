def create_image_montage(image_list, nrow, ncol, wspace=0.1, hspace=0.1):
    # Assume image_list is a list of image paths
    images = [plt.imread(img) for img in image_list]

    # Get the size of the first image
    size = images[0].shape

    # Check if all images are the same size
    for img in images:
        if img.shape != size:
            raise ValueError('Not all images have the same size.')

    # Compute number of images
    n = len(images)

    # Create a montage
    fig, axs = plt.subplots(nrow, ncol)

    for i, ax in enumerate(axs.flat):
        if i < n:
            # Display image if index is less than number of images
            ax.imshow(images[i])
        else:
            # Otherwise hide axes
            ax.axis('off')

        # Remove axis ticks
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    # Remove whitespace between subplots
    plt.tight_layout()
    plt.subplots_adjust(wspace=wspace, hspace=hspace)
    plt.show()


## More complicated functions if needed
# from PIL import Image
# import numpy as np
# import matplotlib.pyplot as plt

# def create_image_montage(*args, **kwargs):
#     if len(args) == 0:
#         raise ValueError('Nothing provided to montage')

#     if isinstance(args[0], list): # if the first argument is a list
#         lengths = [len(arg) for arg in args if isinstance(arg, list)]
#         if len(set(lengths)) > 1:
#             raise ValueError('Not all inputs have the same length.')
        
#         sizes = [np.array(Image.open(arg[0]).size) for arg in args if isinstance(arg, list)]
#         if len(set(tuple(size) for size in sizes)) > 1:
#             raise ValueError('Not all lists have same sizes of images.')
        
#         for i in range(lengths[0]):
#             img_list = [arg[i] for arg in args if isinstance(arg, list)]
#             montage(img_list, **kwargs)
#     else:
#         sizes = [np.array(Image.open(arg).size) for arg in args]
#         if len(set(tuple(size) for size in sizes)) > 1:
#             raise ValueError('Not all inputs have the same size.')

#         montage(args, **kwargs)

# def montage(image_list, **kwargs):
#     images = [Image.open(img) for img in image_list]
    
#     kwargs.setdefault('BackgroundColor', 'black')
#     kwargs.setdefault('BorderSize', (0, 0))
#     kwargs.setdefault('Size', (1, len(images)))

#     fig, axs = plt.subplots(*kwargs['Size'])
    
#     for ax, img in zip(np.ravel(axs), images):
#         ax.imshow(np.array(img))
#         ax.axis('off')
    
#     plt.subplots_adjust(wspace=0, hspace=0)
#     plt.show()
