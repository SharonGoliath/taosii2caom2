Time:

# mjd start

>>> t = Time(1561057509.599084, format='unix')
>>> t.format = 'mjd'
>>> t.value
58654.795249989395
>>> 58654.795249989395
58654.795249989395

# mjd end

>>> t = Time(0.05000000074505806+1561057509.599084, format='unix')
>>> t.format = 'mjd'
>>> t.value
58654.7952505681

start = RefCoord(0.5, mjd_start)
end = RefCoord(1.5, mjd_end)


telescope location:

>>> from astropy.coordinates import EarthLocation
>>> x = EarthLocation.of_site('spm')
Downloading http://data.astropy.org/coordinates/sites.json
|=============================================================================|  23k/ 23k (100.00%)         0s
>>> x
<EarthLocation (-2354953.99637757, -4940160.36363812, 3270123.70695983) m>
>>>


energy:

min = 400 nm
max = 800 nm
central = min + max / 2.0 = 1200.0 / 2.0 == 600.0 nm
fwhm = 400 nm
ref_coord1 = RefCoord(0.5, central_wl - fwhm / 2.0) == 100.0 nm
ref_coord2 = RefCoord(1.5, central_wl + fwhm / 2.0) == 500.0 nm

Position:

# assume equatorial coordinates

>>> from astropy import wcs
>>> w = wcs.WCS(naxis=2)
>>> w.wcs.crpix = [-3850, 2310]
>>> w.wcs.crval = [72.0, 20.77222]
>>> w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
>>> print(w.to_header())
