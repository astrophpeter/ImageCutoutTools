# P. McGill UCSC 2021
# utiliy functions to get image cutouts from various surveys

from astropy.io import fits
import requests
import xml.etree.ElementTree as ET
from astropy.io.votable import parse
import io

ra = 150.0
dec = -30
size=0.001

base_url = 'http://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia?'
query_url = base_url + f'{ra},{dec}&SIZE={size}'


xml = requests.get(url=query_url).content
print(xml)
votable = parse(io.BytesIO(xml))
print(votable.resources)
table = votable.resources[0].tables[0]

# only take the first epoch if J H K fits files, (:3)
fits_mask = table.array['format'] == 'image/fits'
fits_url = table.array['download'][fits_mask][:3]
bands = table.array['band'][fits_mask][:3]
#print({band: url for band, url in zip(bands,fits_url)})


#print(irsa)

def twomass_cutout(ra, dec, size=0.001):
    """

    :param ra: Right acension [Deg]
    :param dec: declinatoin [Deg
    :param size: size of cutout
    :return: fits image
    """

    return 0.0