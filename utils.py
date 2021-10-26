from collections import namedtuple
import yaml
import numpy as np

Survey = namedtuple('Survey', 'name hips_id')


def survey_list(survey_metadata_path=None):
    """
    Build a list of survey objects from a metadata file

    :param survey_metadata_path: Path to the survey meta data
    :return:
    """

    with open(survey_metadata_path, "r") as stream:
        survey_metadata = yaml.safe_load(stream)

    return [Survey(name=name, hips_id=survey_metadata[name]['hips_id']) for name in survey_metadata]

def image_contains_data(image=None):
    """
    Checks if image contains data and is not empty

    :param image: fits image
    :return : true if image is contains data, false otherwise
    """
    return not np.all(np.isnan(image[0].data))





print(survey_list('survey_metadata.yml'))