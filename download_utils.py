# P. McGill UCSC 2021
# utiliy functions to get image cutouts from various surveys

import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astroquery.hips2fits import hips2fits
from regions import CircleSkyRegion
from astroquery.cds import cds
import yaml
from astro_ghost.ghostHelperFunctions import getTransientHosts, getGHOST
from urllib3.exceptions import ReadTimeoutError
from astropy.units import Quantity

#supernova_position = SkyCoord(ra=188.5126408, dec=7.6991489, unit='deg')
#host = find_host_data(supernova_position=supernova_position)
#print(host)

#£with open("survey_metadata.yml", "r") as stream:
#    data = yaml.safe_load(stream)

#for image in data:
#    print(position_in_footprint(position=SkyCoord(ra=130, dec=30, unit='deg'),
#                                survey=data[image]['hips_id']))

#position = SkyCoord(ra=10, dec=30, unit='deg')
#print(find_host_position(position))

from astropy.io import fits
import io
import requests
from astropy.io.votable import parse
#def twomass_cutouts(position, image_size=0.0001):
#    """
#    Get the fits cutouts from the 2MASS survey

#    :param position: On Sky position of the centre of the cutout
#    :param image_size: Size of the cutout
#    :return: fits images of J, H, and K band images
#    """
#    base_url = 'http://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?'
#    query_url = base_url + f'POS={position.ra.degree},{position.dec.degree}&SIZE={image_size}'

#    xml = requests.get(url=query_url).content
#    votable = parse(io.BytesIO(xml))
#    table = votable.resources[0].tables[0]

#    # only take the first epoch if J H K fits files, (:3)
#     # TODO: grab the image with the longest exposure time
#    fits_mask = table.array['format'] == 'image/fits'
#    fits_url = table.array['download'][fits_mask][:3]
#    bands = table.array['band'][fits_mask][:3]
#    return {band: fits.open(url) for band, url in zip(bands, fits_url)}

#position = SkyCoord(ra=22, dec=11, unit='degree')
#all_cutouts(position)