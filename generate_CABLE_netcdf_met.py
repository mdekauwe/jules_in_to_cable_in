#!/usr/bin/env python

"""
Turn the JULES input file into a CABLE netcdf file.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (16.10.2019)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import pandas as pd
import numpy as np
import netCDF4 as nc
import datetime

def main(in_fname, out_fname, lat, lon, start_year, canopy_height,
         reference_height, co2_conc):

    DEG_2_KELVIN = 273.15
    SW_2_PAR = 2.3
    PAR_2_SW = 1.0 / SW_2_PAR
    HLFHR_2_SEC = 1.0 / 1800.

    header=["SWdown","LWdown","Rainf","Snowf","Tair","Wind","PSurf","Qair"]
    df = pd.read_csv(in_fname, comment='#', delimiter=" ", names=header)

    ndim = 1
    n_timesteps = len(df)
    times = []
    secs = 0.0
    for i in range(n_timesteps):
        times.append(secs)
        secs += 1800.

    # create file and write global attributes
    f = nc.Dataset(out_fname, 'w', format='NETCDF4')
    f.description = 'Loobos met data, created by Martin De Kauwe'
    f.history = "Created by: %s" % (os.path.basename(__file__))
    f.creation_date = "%s" % (datetime.datetime.now())
    f.contact = "mdekauwe@gmail.com"

    # set dimensions
    f.createDimension('time', None)
    f.createDimension('z', ndim)
    f.createDimension('y', ndim)
    f.createDimension('x', ndim)
    #f.Conventions = "CF-1.0"

    # create variables
    time = f.createVariable('time', 'f8', ('time',))
    time.units = "seconds since %s-01-01 00:00:00" % (start_year)
    time.long_name = "time"
    time.calendar = "standard"

    z = f.createVariable('z', 'f8', ('z',))
    z.long_name = "z"
    z.long_name = "z dimension"

    y = f.createVariable('y', 'f8', ('y',))
    y.long_name = "y"
    y.long_name = "y dimension"

    x = f.createVariable('x', 'f8', ('x',))
    x.long_name = "x"
    x.long_name = "x dimension"

    latitude = f.createVariable('latitude', 'f8', ('y', 'x',))
    latitude.units = "degrees_north"
    latitude.missing_value = -9999.
    latitude.long_name = "Latitude"

    longitude = f.createVariable('longitude', 'f8', ('y', 'x',))
    longitude.units = "degrees_east"
    longitude.missing_value = -9999.
    longitude.long_name = "Longitude"

    SWdown = f.createVariable('SWdown', 'f8', ('time', 'y', 'x',))
    SWdown.units = "W/m^2"
    SWdown.missing_value = -9999.
    SWdown.long_name = "Surface incident shortwave radiation"
    SWdown.CF_name = "surface_downwelling_shortwave_flux_in_air"

    Tair = f.createVariable('Tair', 'f8', ('time', 'z', 'y', 'x',))
    Tair.units = "K"
    Tair.missing_value = -9999.
    Tair.long_name = "Near surface air temperature"
    Tair.CF_name = "surface_temperature"

    Rainf = f.createVariable('Rainf', 'f8', ('time', 'y', 'x',))
    Rainf.units = "mm/s"
    Rainf.missing_value = -9999.
    Rainf.long_name = "Rainfall rate"
    Rainf.CF_name = "precipitation_flux"

    Snowf = f.createVariable('Snowf', 'f8', ('time', 'y', 'x',))
    Snowf.units = "mm/s"
    Snowf.missing_value = -9999.
    Snowf.long_name = "Snowfall rate"
    Snowf.CF_name = "snowfall_flux"

    Qair = f.createVariable('Qair', 'f8', ('time', 'z', 'y', 'x',))
    Qair.units = "kg/kg"
    Qair.missing_value = -9999.
    Qair.long_name = "Near surface specific humidity"
    Qair.CF_name = "surface_specific_humidity"

    Wind = f.createVariable('Wind', 'f8', ('time', 'z', 'y', 'x',))
    Wind.units = "m/s"
    Wind.missing_value = -9999.
    Wind.long_name = "Scalar windspeed" ;
    Wind.CF_name = "wind_speed"

    PSurf = f.createVariable('PSurf', 'f8', ('time', 'y', 'x',))
    PSurf.units = "Pa"
    PSurf.missing_value = -9999.
    PSurf.long_name = "Surface air pressure"
    PSurf.CF_name = "surface_air_pressure"

    LWdown = f.createVariable('LWdown', 'f8', ('time', 'y', 'x',))
    LWdown.units = "W/m^2"
    LWdown.missing_value = -9999.
    LWdown.long_name = "Surface incident longwave radiation"
    LWdown.CF_name = "surface_downwelling_longwave_flux_in_air"

    CO2 = f.createVariable('CO2air', 'f8', ('time', 'z', 'y', 'x',))
    CO2.units = "ppm"
    CO2.missing_value = -9999.
    CO2.long_name = ""
    CO2.CF_name = ""

    #LAI = f.createVariable('LAI', 'f8', ('time', 'y', 'x'))
    #LAI.setncatts({'long_name': u"Leaf Area Index",})

    hc = f.createVariable('hc', 'f8', ('y', 'x'))
    hc.units = "m"
    hc.missing_value = -9999.
    hc.long_name = "canopy height"

    za = f.createVariable('za', 'f8', ('y', 'x',))
    za.units = "m"
    za.missing_value = -9999.
    za.long_name = "reference height"

    # write data to file
    x[:] = ndim
    y[:] = ndim
    z[:] = ndim
    time[:] = times
    latitude[:] = lat
    longitude[:] = lon

    SWdown[:,0,0] = df.SWdown.values.reshape(n_timesteps, ndim, ndim)
    LWdown[:,0,0] = df.LWdown.values.reshape(n_timesteps, ndim, ndim)
    Tair[:,0,0,0] = df.Tair.values.reshape(n_timesteps, ndim, ndim, ndim)
    Rainf[:,0,0] = df.Rainf.values.reshape(n_timesteps, ndim, ndim)
    Snowf[:,0,0] = df.Snowf.values.reshape(n_timesteps, ndim, ndim)
    Qair[:,0,0,0] = df.Qair.values.reshape(n_timesteps, ndim, ndim, ndim)
    Wind[:,0,0,0] = df.Wind.values.reshape(n_timesteps, ndim, ndim, ndim)
    PSurf[:,0,0] = df.PSurf.values.reshape(n_timesteps, ndim, ndim)
    CO2[:,0,0] = co2_conc

    #LAI[:,0,0] = df.lai.values.reshape(n_timesteps, ndim, ndim)
    hc[:] = canopy_height
    za[:] = reference_height
    f.close()

def convert_rh_to_qair(rh, tair, press):
    """
    Converts relative humidity to specific humidity (kg/kg)

    Params:
    -------
    tair : float
        deg C
    press : float
        pa
    rh : float
        %
    """

    # Sat vapour pressure in Pa
    esat = calc_esat(tair)

    # Specific humidity at saturation:
    ws = 0.622 * esat / (press - esat)

    # specific humidity
    qair = (rh / 100.0) * ws

    return qair

def calc_esat(tair):
    """
    Calculates saturation vapour pressure

    Params:
    -------
    tair : float
        deg C

    Reference:
    ----------
    * Jones (1992) Plants and microclimate: A quantitative approach to
    environmental plant physiology, p110
    """

    esat = 613.75 * np.exp(17.502 * tair / (240.97 + tair))

    return esat


def estimate_lwdown(tairK, rh):
    """
    Synthesises downward longwave radiation based on Tair RH

    Reference:
    ----------
    * Abramowitz et al. (2012), Geophysical Research Letters, 39, L04808

    """
    zeroC = 273.15

    sat_vapress = 611.2 * np.exp(17.67 * ((tairK - zeroC) / (tairK - 29.65)))
    vapress = np.maximum(5.0, rh) / 100. * sat_vapress
    lw_down = 2.648 * tairK + 0.0346 * vapress - 474.0

    return lw_down


if __name__ == "__main__":


    start_year = 1997
    lat = 52.16658
    lon = 5.74356
    canopy_height = 15.5
    reference_height = 27.0
    in_fname = "raw/Loobos_1997.dat"
    out_fname = "Loobos_1997.nc"
    main(in_fname, out_fname, lat, lon, start_year, canopy_height,
         reference_height, co2_conc=400.0)
