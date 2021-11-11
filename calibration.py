import numpy as np


def pixel_to_magnitude(image_value, survey):
    """
    Convert between image pixel units to magnitude

    Parmeters
    ---------

    Returns
    -------
    magnitude
    """
    zero_point = survey.magnitude_zero_point
    flux = image_value if survey.unit == 'counts' else image_value / survey.exposure_time
    magnitude = zero_point - 2.5 * np.log10(flux)

    return 0.0
