# utils to extract galaxy photometry of the host from the cutout images
# P. McGill

from photutils.segmentation import detect_threshold, detect_sources,deblend_sources,SourceCatalog
from download_utils import cutout, find_host_data
from astropy.coordinates import SkyCoord
from photutils.background import Background2D, MedianBackground
from astropy.wcs import WCS
import numpy as np
import yaml
import astropy.units as u
from utils import survey_list, image_contains_data
from urllib3.exceptions import ReadTimeoutError
from photutils.aperture import EllipticalAperture
from astropy.wcs.utils import skycoord_to_pixel
from astropy.units import Quantity
#with open("survey_metadata.yml", "r") as stream:
#    survey_metadata = yaml.safe_load(stream)
#
##supernova_position = SkyCoord(ra=188.5126408, dec=7.6991489, unit='deg')
#image = cutout(position=supernova_position, survey=survey_metadata['PanSTARRS_z']['hips_id'],fov=0.3, width=300, height=300)
#data = image[0].data
#threshold = detect_threshold(data, nsigma=2.)


def build_source_catalog(image=None, background_estimator=MedianBackground(), box_size=50, threshhold_sigma=2.0, npixels=10):
    """
    takes images data and builds a source cataloge

    :param image_data: fits image data
    :param background_estimator:
    :param box_size:
    :param threshhold_sigma:
    :return source catalog:
    """

    image_data = image[0].data

    background = Background2D(image_data, box_size, bkg_estimator=background_estimator)
    background_subtracted_data = image_data - background.background
    threshold = threshhold_sigma * background.background_rms

    segmentation = detect_sources(background_subtracted_data, threshold, npixels=npixels)
    deblended_segmentation = deblend_sources(background_subtracted_data, segmentation, npixels=npixels)

    return SourceCatalog(background_subtracted_data, deblended_segmentation)

def match_sources_to_host(image, source_catalog, host_position):
    """
    Matches the host galaxy to a source in the source cataloge.

    :param image:
    :param source_catalog:
    :param host_position:
    :return Host galaxy information:
    """

    world_coordinate_system = WCS(image[0].header)
    host_x_pixel, host_y_pixel = world_coordinate_system.world_to_pixel(host_position)
    source_x_pixels, source_y_pixels = source_catalog.xcentroid, source_catalog.ycentroid
    closest_source_index = np.argmin(np.hypot(host_x_pixel - source_x_pixels, host_y_pixel - source_y_pixels))
    return source_catalog[closest_source_index]

def elliptical_sky_aperture(source_cat, world_coordinate_system, r=3.0):
    """Constructs elliptical sky aperture object
    """
    center = (source_cat.xcentroid, source_cat.ycentroid)
    semi_major_axis = source_cat.semimajor_sigma.value * r
    semi_minor_axis = source_cat.semiminor_sigma.value * r
    orientation_angle = source_cat.orientation.to(u.rad).value
    pixel_aperture = EllipticalAperture(center, semi_major_axis, semi_minor_axis, theta=orientation_angle)
    return pixel_aperture.to_sky(world_coordinate_system)

def host_aperture(image, host_position):
    source_catalog = build_source_catalog(image)
    host_catalog_data = match_sources_to_host(image, source_catalog, host_position)
    world_coordinate_system = WCS(image[0].header)
    return elliptical_sky_aperture(host_catalog_data, world_coordinate_system)


def find_largest_aperture(host_position, survey_images):
    """Find the largest fitted elliptical aperture from a list of images

    :param host_position:
    :param surveys:
    :return:
    """
    for image in survey_images:
        pass

import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm



def run_forced_host_photometry(supernova_position, survey_metadata_path=None):
    """
    main program for running host forced photometry

    :param supernova_position:
    :return:
    """

    host_data = find_host_data(supernova_position=supernova_position)
    host_position = host_data['position']

    all_surveys = survey_list(survey_metadata_path=survey_metadata_path)
    images = download_image_data(host_position, all_surveys)

    apertures = {name: host_aperture(image, host_position) for name, image in images.items()}
    plt.imshow(images['PanSTARRS_g'][0].data, cmap='gray', norm=PowerNorm(0.5))
    wcs = WCS(images['PanSTARRS_g'][0].header)

    for name, aperture in apertures.items():
        pixel_aperture = aperture.to_pixel(wcs)
        pixel_aperture.plot(color='black')
        sn_pixel = skycoord_to_pixel(supernova_position, wcs)
        host_pixel = skycoord_to_pixel(host_position, wcs)

    apertures['PanSTARRS_g'].to_pixel(wcs).plot(color='red', label='This image aperture')
    plt.scatter(sn_pixel[0], sn_pixel[1], label='SN position')
    plt.scatter(host_pixel[0], host_pixel[1], label='GOST host position')
    plt.title('PanSTARRS_g')
    plt.legend()
    plt.show()
supernova_position = SkyCoord(ra=188.5148408, dec=7.6991489, unit='deg')
print(run_forced_host_photometry(supernova_position, survey_metadata_path='survey_metadata.yml'))

#import astropy.units as u
#from photutils.aperture import EllipticalAperture
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

#cat = build_source_catalog(image)
#cat = match_sources_to_host(image, cat, supernova_position)
#position = (cat.xcentroid, cat.ycentroid)

#r = 1.0  # approximate isophotal extent
#a = cat.semimajor_sigma.value * r
#b = cat.semiminor_sigma.value * r
#theta = cat.orientation.to(u.rad).value
#apertures = EllipticalAperture(position, a, b, theta=theta)

#plt.imshow(data, origin='lower', cmap='gray', interpolation='nearest', norm=PowerNorm(0.5))
#apertures.plot(color='#d62728')
#plt.show()
#host = match_sources_to_host(image, cat, supernova_position)
#wcs = WCS(image[0].header)
#print(wcs.pixel_scale_matrix[0,0] / wcs.proj_plane_pixel_scales()[0].value)
#print(proj_plane_pixel_scales(wcs))
