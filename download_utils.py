# P. McGill UCSC 2021
# utiliy functions to get image cutouts from various surveys

from astropy.io import fits
import astropy.units as u
import requests
from astropy.io.votable import parse
from astropy.coordinates import SkyCoord
from astroquery.hips2fits import hips2fits
import io

def twomass_cutouts(position, image_size=0.0001):
    """
    Get the fits cutouts from the 2MASS survey

    :param position: On Sky position of the centre of the cutout
    :param image_size: Size of the cutout
    :return: fits images of J, H, and K band images
    """
    base_url = 'http://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?'
    query_url = base_url + f'POS={position.ra.degree},{position.dec.degree}&SIZE={image_size}'

    xml = requests.get(url=query_url).content
    votable = parse(io.BytesIO(xml))
    table = votable.resources[0].tables[0]

    # only take the first epoch if J H K fits files, (:3)
    # TODO: grab the image with the longest exposure time
    fits_mask = table.array['format'] == 'image/fits'
    fits_url = table.array['download'][fits_mask][:3]
    bands = table.array['band'][fits_mask][:3]
    return {band: fits.open(url) for band, url in zip(bands, fits_url)}

def cutout(position=None, survey=None, fov=0.001, width=1000, height=1000):
    """
    Get image cutout for a survey in a particular band

    :param position: On Sky position of the center of the cutout
    :param survey: Survey and band string, for allowed options see https://aladin.u-strasbg.fr/hips/list
    :param image_fov: Field of view in decimal degrees
    :param image_width: image width in pixels
    :param image_height: image height in pixels
    :return: fits image
    """
    fits = hips2fits.query(hips=survey,
                           ra=position.ra,
                           dec=position.dec,
                           width=width, height=height, fov=fov * (u.deg),
                           projection='TAN',
                           format='fits')
    return fits

def all_cutouts(position):

    names = {'2MASS/': ['H', 'J', 'K'],
             'WISE/':['W1', 'W2', 'W3', 'W4'],
             'GALEXGR6/AIS/': ['FUV','NUV'],
             'PanSTARRS/DR1/': ['g','r', 'i', 'z', 'y'],
             'SDSS9/' : ['i', 'r', 'u', 'z'],
             'DES-DR1/' : ['Y', 'g', 'r', 'i', 'z']             }

    for survey in names:
        filters = names[survey]
        for filter in filters:
            print(survey + filter)
            print(cutout(position=position, survey=survey + filter))

position = SkyCoord(ra=130, dec=-30, unit='degree')
all_cutouts(position)