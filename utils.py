from collections import namedtuple
import yaml
import numpy as np
from urllib3.exceptions import ReadTimeoutError
from astropy.units import Quantity
import astropy.units as u
from astroquery.hips2fits import hips2fits


def survey_list(survey_metadata_path):
    """
    Build a list of survey objects from a metadata file.

    Parameters
    ----------
    :survey_metadata_path : str
        Path to a ymal data file containing survey metadata

    Returns
    -------
    :list of surveys: list[Survey]
        List of survey objects

    """
    with open(survey_metadata_path, "r") as stream:
        survey_metadata = yaml.safe_load(stream)

    # get first survey from the metadata in order to infer the data field names
    survey_name = list(survey_metadata.keys())[0]
    data_fields = list(survey_metadata[survey_name].keys())

    # create a named tuple class with all the survey data fields as attributes
    # including the survey name
    Survey = namedtuple('Survey', ['name'] + data_fields)

    survey_list = []
    for name in survey_metadata:
        field_dict = {field: survey_metadata[name][field] for field in data_fields}
        field_dict['name'] = name
        survey_list.append(Survey(**field_dict))

    return survey_list


def image_contains_data(image):
    """
    Checks if a fits image contains data and is not empty

    Parameters
    ----------
    :image : :class:`~astropy.io.fits.HDUList`
        Fits image to check.

    Returns
    -------
    : contains data: bool
        True if the image contains non-nan data, false otherwise.
    """
    return not np.all(np.isnan(image[0].data))


def cutout(position, survey, fov=Quantity(0.2, unit='deg')):
    """
    Download image cutout data from a survey.

    Parameters
    ----------
    :position : :class:`~astropy.coordinates.SkyCoord`
        Target centre position of the cutout image to be downloaded.
    :survey : :class: Survey
        Named tuple containing metadata for the survey the image is to be
        downloaded from.
    :fov : :class:`~astropy.units.Quantity`,
    default=Quantity(0.2,unit='deg')
        Field of view of the cutout image, angular length of one of the sides
        of the square cutout. Angular astropy quantity. Default is angular
        length of 0.2 degrees.

    Returns
    -------
    :cutout : :class:`~astropy.io.fits.HDUList` or None
        Image cutout in fits format or if the image cannot be download due to a
        `ReadTimeoutError` None will be returned.

    """
    num_pixels = fov.to(u.arcsec).value / survey.pixel_size_arcsec
    try:
        fits = hips2fits.query(hips=survey.hips_id, ra=position.ra,
                               dec=position.dec, width=int(num_pixels),
                               height=int(num_pixels), fov=fov,
                               projection='TAN', format='fits')
    except ReadTimeoutError:
        print(f'Conection timed out, could not download {survey.name} data')
        fits = None
    return fits


def download_image_data(position, survey_list, fov=Quantity(0.2, unit='deg')):
    """
    Download all available imaging from a list of surveys

    Parameters
    ----------

    :position : :class:`~astropy.coordinates.SkyCoord`
        Target centre position of the cutout image to be downloaded.
    :survey_list : list[Survey]
        List of surveys to download data from
    :fov : :class:`~astropy.units.Quantity`,
    default=Quantity(0.2,unit='deg')
        Field of view of the cutout image, angular length of one of the sides
        of the square cutout. Angular astropy quantity. Default is angular
        length of 0.2 degrees.

    Returns
    -------
    :images dictionary : dict[str: :class:`~astropy.io.fits.HDUList`]
        Dictionary of images with the survey names as keys and fits images
        as values.
    """
    images = [cutout(position, survey, fov=fov) for survey in survey_list]
    return {survey.name: image for survey, image in zip(survey_list, images)
              if image is not None and image_contains_data(image)}