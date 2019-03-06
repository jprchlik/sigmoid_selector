import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
import pandas as pd
from datetime import datetime

data = pd.read_csv('SigmoidCatalogAll_filament.csv')

#remove first 30 which were not reanalyzed
data = data.iloc[28:,:]
data = data.loc[[isinstance(i,str) for i in data.tobs],:]

dfmt = '%Y-%m-%dT%H:%M:%S.000'

for j in range(data.shape[0]):
    c = SkyCoord(data.iloc[j,:].X*u.arcsec, data.iloc[j,:].Y*u.arcsec, obstime=datetime.strptime(data.iloc[j,:].tobs,dfmt), frame=frames.Helioprojective)
    p = c.transform_to(frames.HeliographicCarrington(obstime=datetime.strptime(data.iloc[j,:].TBEST,dfmt)))
    print('ID = {0:3d} ==> SOL{1}L{3:03.0f}C{2:03.0f}'.format(data.iloc[j,:].ID,datetime.strptime(data.iloc[j,:].TBEST,dfmt).strftime('%Y-%m-%dT%H:%M:%S'),90.-p.lat.value,p.lon.value))

