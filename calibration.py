import numpy as np


def pixel_to_magnitude(pixel_value, survey, exposure_time=None):
    """
    Convert between image pixel units to magnitude

    Parameters
    ---------
    image_value: float: value of pixels in the image to be converted into magnitudes
    survey: :Survey: survey the image was taken from
    exposure_time: exposure time in seconds of the the image

    Returns
    -------
    magnitude: float: value of the magnitude
    """

    if survey.unit == 'counts' and exposure_time is None:
        raise ValueError(f'Image exposure time needs to be specified to calculate magnitude for {survey.name}')

    flux = pixel_value if survey.unit == 'counts' else image_value / exposure_time
    return survey.magnitude_zero_point - 2.5 * np.log10(flux)
