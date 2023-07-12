import seaborn as sns
import random
import numpy as np
import matplotlib.pyplot as plt

def generate_colors(num_colors, seed=None):
    """Generates a set of discrete colors."""
    if seed is not None:
        random.seed(seed)
    colors = ['#%06x' % random.randint(0, 0xFFFFFF) for _ in range(num_colors)]
    return colors

def channelDensity(I_list):
    # Convert the list to a numpy array and reshape it for easier access
    I_list_arr = np.array(I_list).reshape((-1, *I_list[0].shape))
    
    # Get a list of colors to use
    colors = generate_colors(num_colors=len(I_list), seed=12535)

    fig, axs = plt.subplots(3, 1, figsize=(10, 8))

    for idx, channel_name in enumerate(['R', 'G', 'B']):
        for i, color in enumerate(colors):
            sns.kdeplot(I_list_arr[i, :, :, idx].ravel(), color=color, alpha=0.7, label=f'Image {i+1}', ax=axs[idx])
        axs[idx].set_title(channel_name)
        axs[idx].legend()

    plt.tight_layout()
    plt.show()
