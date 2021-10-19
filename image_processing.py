# utils to extract galaxy photometry of the host from the cutout images
# P. McGill

from photutils.segmentation import detect_threshold, detect_sources,deblend_sources,SourceCatalog
from download_utils import cutout
from astropy.coordinates import SkyCoord
from photutils.background import Background2D, MedianBackground
from astropy.wcs import WCS
import numpy as np

supernova_position = SkyCoord(ra=188.5126408, dec=7.6991489, unit='deg')
image = cutout(position=supernova_position, survey='CDS/P/2MASS/H')
data = image[0].data
threshold = detect_threshold(data, nsigma=2.)


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


#cat = build_source_catalog(image)
#host = match_sources_to_host(image, cat, supernova_position)

#print(host.kron_flux)