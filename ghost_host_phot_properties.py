import utils as util
import pandas as pd
from astropy.coordinates import SkyCoord
import numpy as np

data = pd.read_csv('data/sherlock_ghost.csv').dropna()

surveys = util.survey_list('survey_metadata.yml')
survey = [survey for survey in surveys if survey.name == 'PanSTARRS_g'][0]
print(survey)

aperture_sma = np.zeros(len(data))

for row, _ in data.iterrows():
        ghost_ra, ghost_dec = data['ghost_ra'][row],data['ghost_dec'][row]
        host_pos = SkyCoord(ra=ghost_ra, dec=ghost_dec, unit='deg')
        image = util.cutout(host_pos,survey)
        aperture_sma[row] = util.construct_aperture(image, host_pos).a.value
        data['aperture_semi_major_axis'] = aperture_sma
        data.to_csv('data/sherklock_ghost_host_aperture.csv', index=False)


