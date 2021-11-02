# Plotting functions
# P. McGIll 2021
from matplotlib.colors import PowerNorm
import matplotlib.pyplot as plt
import numpy as np
from utils import download_image_data, survey_list, construct_aperture
from astropy.coordinates import SkyCoord


def plot_image(image, axis):
    """Plots fits image"""
    axis.imshow(image[0].data, norm=PowerNorm(0.2), cmap='gray', origin='lower')



def plot_imaging_data(image_dict, aperture_dict, save_dir='scratch'):
    """Plots images with apertures over data"""

    survey_names = list(image_dict.keys())
    num_surveys = len(survey_names)
    n_rows = int(np.sqrt(num_surveys)) + 1
    fig, axs = plt.subplots(n_rows, n_rows)

    for survey ,ax in zip(survey_names, axs.reshape(-1)):
        plot_image(image_dict[survey], ax)
        ax.axis('off')
        ax.set_title(survey)

    plt.show()


supernova_position = SkyCoord(ra=188.5148408, dec=7.6991489, unit='deg')
survey_list = survey_list('survey_metadata.yml')

images = download_image_data(supernova_position, survey_list)
#apertures = {name: construct_aperture(image, supernova_position)
        #     for name, image in images.items()}

plot_imaging_data(images, {})







