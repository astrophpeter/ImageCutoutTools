import pandas as pd
import matplotlib.pyplot as plt
from utils import find_host_data
from astropy.coordinates import SkyCoord

test_data = pd.read_csv('data/foley+18_table1_clean.csv')

ghost_true_sep = []



for _, row in test_data.iterrows():
    host_data = find_host_data(SkyCoord(row['sn_ra'], row['sn_dec'], unit='deg'))
    if host_data is not None:
        ghost_position = host_data['position']
        true_host = SkyCoord(ra=row['host_ra'], dec=row['host_dec'], unit='deg')
        ghost_true_sep.append(true_host.separation(ghost_position).arcsecond)



plt.hist(ghost_true_sep, bins=20)
plt.xlabel('GHOST vs True Host Separation [arcsec]')
plt.ylabel('Count')
plt.title('Foley+18 Foundation Supernova Survey')
plt.savefig('GHOST_test.png', dpi=200)










