from astropy.io.fits import Header
from astropy.table import Table
from astropy.wcs import WCS

t_window = Table.read('./20190805T024026_f060_s00001.h5', path='window')
t_wcs = Table.read('20190805T024026_f060_s00001.h5', path='/wcs/cdmatrix')
t_c = Table.read('20190805T024026_f060_s00001.h5', path='/catalog')

tele = 0
h = Header()
h['NAXIS1'] = t_window['X1'].data[tele] - t_window['X0'].data[tele] + 1
h['NAXIS2'] = t_window['Y1'].data[tele] - t_window['Y0'].data[tele] + 1
h['CRVAL1'] = t_wcs['CRVAL1'].data[tele]
h['CRVAL2'] = t_wcs['CRVAL2'].data[tele]
# h['CRVAL1'] = t_c['RA'].data[0]
# h['CRVAL2'] = t_c['DEC'].data[0]
h['CRPIX1'] = t_wcs['CRPIX1'].data[tele]
h['CRPIX2'] = t_wcs['CRPIX2'].data[tele]
# h['CRPIX1'] = 1.0
# h['CRPIX2'] = 1.0
h['CDELT1'] = 4.0 / 3600.0
h['CDELT2'] = 4.0 / 3600.0
h['CD1_1'] = t_wcs['CD1_1'].data[tele]
h['CD1_2'] = t_wcs['CD1_2'].data[tele]
h['CD2_1'] = t_wcs['CD2_1'].data[tele]
h['CD2_2'] = t_wcs['CD2_2'].data[tele]

from caom2 import CoordPolygon2D, ValueCoord2D
import logging
w = WCS()
bounds = CoordPolygon2D()
result = w.calc_footprint(header=h)
for ii in result:
    logging.error(type(ii))
    vertex = ValueCoord2D(ii[0], ii[1])
    bounds.vertices.append(vertex)

