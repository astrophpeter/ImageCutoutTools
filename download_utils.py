# P. McGill UCSC 2021
# utiliy functions to get image cutouts from various surveys

from astropy.io import fits
import requests
from astropy.io.votable import parse
from astropy.coordinates import SkyCoord
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

position = SkyCoord(ra=130.0, dec=-30.0, unit='degree')
print(twomass_cutouts(position))