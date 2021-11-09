# Plotting functions
# P. McGIll 2021
from matplotlib.colors import PowerNorm
import matplotlib.pyplot as plt
import numpy as np
from utils import download_image_data, survey_list, construct_aperture, cutout
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.units import Quantity

def plot_image(image, axis, color_scale=0.2):
    """Plots fits image"""
    axis.imshow(image[0].data, norm=PowerNorm(color_scale), cmap='gray', origin='lower')



def plot_imaging_data(image_dict, aperture_dict, save_dir='scratch'):
    """Plots images with apertures over data"""

    survey_names = list(image_dict.keys())
    num_surveys = len(survey_names)
    n_rows = int(np.sqrt(num_surveys)) + 1
    fig, axs = plt.subplots(n_rows, n_rows)

    for survey, ax in zip(survey_names, axs.reshape(-1)):
        plot_image(image_dict[survey], ax)
        ax.axis('off')
        ax.set_title(survey)



def plot_sn_host_position(sn_positon=None,
                          ghost_position=None,
                          true_host_position=None,
                          sn_name=None,
                          survey='PanSTARRS_y',
                          out_dir='sherlock_match_plots'):
    """Plot supernova, host and true host position over a survey image"""
    surveys = survey_list(survey_metadata_path='survey_metadata.yml')
    survey_to_plot = [survey_plot for survey_plot in surveys if survey_plot.name == survey][0]
    image = cutout(true_host_position, survey_to_plot, fov=Quantity(0.04, unit='deg'))
    wcs = WCS(image[0].header)

    fig, ax = plt.subplots(1,1)
    plot_image(image, ax, color_scale=1.0)
    ax.axis('off')

    sn_x, sn_y = sn_positon.to_pixel(wcs)
    ghost_x, ghost_y = ghost_position.to_pixel(wcs)
    true_x, true_y = true_host_position.to_pixel(wcs)

    ax.scatter(sn_x, sn_y, label='Supernova')
    ax.scatter(ghost_x, ghost_y, label='GHOST', s=200, facecolors='none', edgecolors='r')
    ax.scatter(true_x, true_y, label='Sherlock Survey Host')
    ax.set_title(sn_name.strip())
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_dir+ '/' + sn_name.strip() + '_host_match.png', dpi=150)
    plt.clf()
    plt.close()












#supernova_position = SkyCoord(ra=188.5148408, dec=7.6991489, unit='deg')
#survey_list = survey_list('survey_metadata.yml')

#images = download_image_data(supernova_position, survey_list)
#apertures = {name: construct_aperture(image, supernova_position)
        #     for name, image in images.items()}
#plot_imaging_data(images, {})







