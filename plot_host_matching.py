import plotting_utils as utils
import pandas as pd
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np

matching_data = pd.read_csv('data/sherlock_ghost.csv').dropna(subset=['sn_ra', 'sn_dec', 'host_ra',
                                                                      'host_dec','ghost_ra','ghost_dec'])

seps = []

print(matching_data)
for _, row in matching_data.iterrows():
    sn_position = SkyCoord(ra=row['sn_ra'], dec=row['sn_dec'], unit='deg')
    ghost_position = SkyCoord(ra=row['ghost_ra'], dec=row['ghost_dec'], unit='deg')
    host_position = SkyCoord(ra=row['host_ra'], dec=row['host_dec'], unit='deg')
    seps.append(ghost_position.separation(host_position).arcsec)
    utils.plot_sn_host_position(sn_positon=sn_position,
                                ghost_position=ghost_position,
                                true_host_position=host_position,
                                sn_name=row['sn_name'])

hist, bins = np.histogram(seps, bins=30)
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
plt.hist(seps, bins=logbins)
plt.xscale('log')
plt.xlabel('Ghost vs True Host Separation [arcsec]')
plt.title('Sherlock match sample')
plt.ylabel('Count')
plt.savefig('sherlock_ghost_test.png')



#    print(ghost_position, host_position)
#    print(ghost_position.separation(host_position).arcsec)
#utils.plot_sn_host_position(sn_positon=sn_position,
#                                ghost_position=ghost_position,
#                                true_host_position=host_position,
#                                sn_name=row['sn name'])





