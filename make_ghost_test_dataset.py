import pandas as pd
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord



def simbad_to_skycoord(simbad_result):
    """
    Takes a simbad query objects and
    returns the sky coordinate of the object.
    """
    ra, dec = simbad_result['RA'].value[0], simbad_result['DEC'].value[0]
    return SkyCoord(ra + ' ' + dec, unit='deg')


raw_data = pd.read_fwf('data/foley+18_table1.txt', delim_whitespace='True', header=None)
sn_name_one, sn_name_two, host_names = raw_data[0].values, raw_data[1].values, raw_data[3].values



sn_coords = []
host_coords = []
sn_names = []

for sn_one, sn_two, name in zip(sn_name_one, sn_name_two, host_names):

    try:
        host_data = Simbad.query_object(name)
    except:
        host_data = None

    if host_data is not None and len(host_data['RA'].value[0]) > 0:
        host_coords.append(simbad_to_skycoord(host_data))
    else:
        host_coords.append(None)

    sn_name = 'sn' + sn_one
    sn_data = Simbad.query_object(sn_name)
    if sn_data is not None:
        sn_coords.append(simbad_to_skycoord(sn_data))
        sn_names.append(sn_one)
    else:

        sn_data = Simbad.query_object(sn_two)
        if sn_data is not None:
            sn_coords.append(simbad_to_skycoord(sn_data))
            sn_names.append(sn_two)
        else:
            sn_coords.append(None)
            sn_names.append(None)


sn_ra, sn_dec, host_ra, host_dec = [], [], [], []

for sn, host in zip(sn_coords, host_coords):

    if sn is not None and host is not None:
        sn_ra.append(sn.ra.degree)
        sn_dec.append(sn.dec.degree)
        host_ra.append(host.ra.degree)
        host_dec.append(host.dec.degree)
    else:
        sn_ra.append(None)
        sn_dec.append(None)
        host_ra.append(None)
        host_dec.append(None)


data = pd.DataFrame({'sn name': sn_names,
                     'sn_ra' : sn_ra,
                     'sn_dec' : sn_dec,
                     'host_ra' : host_ra,
                     'host_dec' : host_dec})

data.to_csv('data/foley+18_table1_clean.csv', index=False)