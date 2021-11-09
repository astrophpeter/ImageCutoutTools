import pandas as pd
import matplotlib.pyplot as plt
from utils import find_host_data
from astropy.coordinates import SkyCoord



test_data = pd.read_csv('data/foley+18_table1_clean.csv')

ghost_host_ra = []
ghost_host_dec = []



for _, row in test_data.iterrows():

    if not pd.isna(row['sn_ra']):
        host_data = find_host_data(SkyCoord(row['sn_ra'],
                                            row['sn_dec'],
                                            unit='deg'))
    else:
        host_data = None

    if host_data is not None:
        ghost_position = host_data['position']
        ghost_host_ra.append(ghost_position.ra)
        ghost_host_dec.append(ghost_position.dec)
    else:
        ghost_host_ra.append(None)
        ghost_host_dec.append(None)


test_data['ghost_host_ra'], test_data['ghost_host_dec'] = ghost_host_ra, ghost_host_dec
print(test_data)
test_data.to_csv('data/foley+18_ghost.csv', index=False)














