#!/usr/bin/env python3
# Wenchang Yang (wenchang@princeton.edu)
# Sat Aug 24 15:02:58 EDT 2019
#readme: update from ERA5.py: potential_intensity -> potential_intensity_tcpypi; entropy_deficit also updated with new arg forGPI2010
#wy2025-03-06: add rename_time and rename_level functions to make new ERA5 data have compatible dimensions as the old data
# updated for DAILY data
if __name__ == '__main__':
    import sys
    wython = "/scratch/gpfs/GEOCLIM/wenchang/wython"
    if wython not in sys.path:
        sys.path.append(wython)
    from misc.timer import Timer
    tt = Timer('start ' + ' '.join(sys.argv))
import os.path#, os, sys
import xarray as xr, numpy as np#, pandas as pd
#import matplotlib.pyplot as plt
from numpy import absolute, exp, log
import calendar

from xtci.shared.entropy_deficit import entropy_deficit
#from xtci.shared.potential_intensity import potential_intensity
from xtci.shared.potential_intensity_tcpypi import potential_intensity 
from xtci.shared.wind_shear import wind_shear
from xtci.shared.absolute_vorticity import absolute_vorticity

if __name__ == '__main__':
    tt.check('end import')

def rename_time(da):
    """make the time dimension name in new data compatible to old data"""
    if 'valid_time' in da.dims:
        da = da.rename(valid_time='time')
        print('[renamed]: valid_time -> time')
    return da
def rename_level(da):
    """make the level dimension name in new data compatible to old data"""
    if 'pressure_level' in da.dims:
        da = da.rename(pressure_level='level')
        da = da.assign_coords(level=da.level.values.astype(int))
        print('[renamed]: pressure_level -> level')
    return da


def do_tci(year, month, odir=None):
    '''calculate TC indices (e.g. GPI, VI) and related variables given ERA5 DAILY reanalysis'''
    print('[year]:', year)
    if odir is None:
        odir = '.'
    ibasename = f'era5.daily.{year:04d}-{month:02d}.nc' 

    # sst and ocean mask
    ifile = f'/tigress/wenchang/data/era5/analysis_wy/sea_surface_temperature/daily/era5.sea_surface_temperature.daily.{year:04d}-{month:02d}.nc'
    sst = xr.open_dataarray(ifile)# units K
    sst = rename_time(sst)
    is_ocean = sst.isel(time=0).drop('time').pipe(lambda x: x*0==0)

    # slp
    ifile = f'/tigress/wenchang/data/era5/analysis_wy/mean_sea_level_pressure/daily/era5.mean_sea_level_pressure.daily.{year:04d}-{month:02d}.nc'
    slp = xr.open_dataarray(ifile) # units Pa
    slp = rename_time(slp)
    
    # t2m
    ifile = f'/tigress/wenchang/data/era5/analysis_wy/2m_temperature/daily/era5.2m_temperature.daily.{year:04d}-{month:02d}.nc'
    t2m = xr.open_dataarray(ifile) # units K
    t2m = rename_time(t2m)

    # RH2m
    ifile = f'/tigress/wenchang/data/era5/analysis_wy/2m_relative_humidity/daily/era5.2m_relative_humidity.daily.{year:04d}-{month:02d}.nc'
    RH2m = xr.open_dataarray(ifile) # units %
    RH2m = rename_time(RH2m)

    # RH600
    ifile = f'/projects/w/wenchang/data/era5/analysis_wy/plevels/rh600/daily/era5.rh600.daily.{year:04d}-{month:02d}.nc'
    RH600 = xr.open_dataarray(ifile) # in %
    RH600 = rename_time(RH600)
    RH600 = rename_level(RH600)

    # Ta
    ifile = f'/projects/w/wenchang/data/era5/analysis_wy/plevels/temperature/daily/era5.temperature.daily.{year:04d}-{month:02d}.nc'
    Ta = xr.open_dataarray(ifile) # units K
    Ta = rename_time(Ta)
    Ta = rename_level(Ta)

    # q 
    ifile = f'/projects/w/wenchang/data/era5/analysis_wy/plevels/specific_humidity/daily/era5.specific_humidity.daily.{year:04d}-{month:02d}.nc'
    q = xr.open_dataarray(ifile) # units K
    q = rename_time(q)
    q = rename_level(q)

    # u200
    ifile = f'/tigress/wenchang/data/era5/analysis_wy/plevels/u200/daily/era5.u200.daily.{year:04d}-{month:02d}.nc'
    u200 = xr.open_dataarray(ifile) # in m/s
    u200 = rename_time(u200)
    u200 = rename_level(u200)

    # u850
    ifile = f'/tigress/wenchang/data/era5/analysis_wy/plevels/u850/daily/era5.u850.daily.{year:04d}-{month:02d}.nc'
    u850 = xr.open_dataarray(ifile) # in m/s
    u850 = rename_time(u850)
    u850 = rename_level(u850)

    # v200
    ifile = f'/tigress/wenchang/data/era5/analysis_wy/plevels/v200/daily/era5.v200.daily.{year:04d}-{month:02d}.nc'
    v200 = xr.open_dataarray(ifile) # in m/s
    v200 = rename_time(v200)
    v200 = rename_level(v200)

    # v850
    ifile = f'/tigress/wenchang/data/era5/analysis_wy/plevels/v850/daily/era5.v850.daily.{year:04d}-{month:02d}.nc'
    v850 = xr.open_dataarray(ifile) # in m/s
    v850 = rename_time(v850)
    v850 = rename_level(v850)

    # vorticity
    ifile = f'/tigress/wenchang/data/era5/analysis_wy/plevels/vort850/daily/era5.vort850.daily.{year:04d}-{month:02d}.nc'
    vort850 = xr.open_dataarray(ifile) # in s**-1
    vort850 = rename_time(vort850)
    vort850 = rename_level(vort850)


    # entropy deficit: (s_m_star - s_m)/(s_sst_star - s_b)
    print('entropy deficit')
    dname = 'chi'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        chi = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        p_m = 600*100 # Pa
        chi = entropy_deficit(
            sst=sst,
            slp=slp,
            Tb=t2m,
            RHb=RH2m/100,
            p_m=p_m,
            Tm=Ta.sel(level=p_m/100).drop('level'),
            RHm=RH600/100 
            ).where(is_ocean)
        chi.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # entropy deficit for GPI2010: (s_b - s_m)/(s_sst_star - s_b)
    print('entropy deficit for GPI2010')
    dname = 'chi_sb'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        chi_sb = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        p_m = 600*100 # Pa
        chi_sb = entropy_deficit(
            sst=sst,
            slp=slp,
            Tb=t2m,
            RHb=RH2m/100,
            p_m=p_m,
            Tm=Ta.sel(level=p_m/100).drop('level'),
            RHm=RH600/100, 
            forGPI2010=True
            ).where(is_ocean)
        chi_sb.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # potential intensity
    print('potential intensity')
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.PI.nc') )
    if os.path.exists(ofile):
        PI = xr.open_dataset(ofile)
        print('[opened]:', ofile)
    else:
        #make sure pressure levels in decrease order, i.e. from bottom to top
        reverse_plevels = lambda x: x.isel(level=slice(-1, None, -1)) if x.level.values[-1] - x.level.values[0] > 0 else x # 
        """
        PI = potential_intensity(
            sst=sst,
            slp=slp.where(is_ocean),
            p=Ta.level.pipe(reverse_plevels),
            T=Ta.pipe(reverse_plevels).where(is_ocean),
            q=q.pipe(reverse_plevels).where(is_ocean),
            dim_x='longitude', dim_y='latitude', dim_z='level'
            )
        """
        PI = potential_intensity(
            sst=sst,
            slp=slp.where(is_ocean),
            p=Ta.level.pipe(reverse_plevels)*100,
            T=Ta.pipe(reverse_plevels).where(is_ocean),
            q=q.pipe(reverse_plevels).where(is_ocean),
            dim_z='level'
            )
        encoding = {dname:{'dtype': 'float32', 'zlib': True, 'complevel': 1} 
            for dname in ('pmin', 'vmax')}
        encoding['iflag'] = {'dtype': 'int32'}
        PI.to_netcdf(ofile, encoding=encoding, unlimited_dims='time')
        print('[saved]:', ofile)

    # wind shear: ( (u200-u850)**2 + (v200-v850)**2 )**0.5
    print('wind shear')
    dname = 'Vshear'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        Vshear = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        Vshear = wind_shear(
            u850,
            v850,
            u200,
            v200
            )
        Vshear.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # ventilation index: Vshear * chi_m /V_PI
    print('ventilation index')
    dname = 'VI'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        VI = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        VI = Vshear*chi/PI.vmax.pipe(lambda x: x.where(x>0))
        VI.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # absolute vorticity at 850hPa
    print('absolute vorticity')
    dname = 'eta'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        eta = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        eta = absolute_vorticity(
            vort850, 
            lat=vort850.latitude
            )
        eta.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)


    # relative humidity at 600hPa in %
    print('relative humidity in %')
    dname = 'H'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        H = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        H = RH600
        H.attrs['long_name'] = '600hPa relative humidity'
        H.attrs['units'] = '%'
        H.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)
    
    # GPI (Emanuel and Nolan 2004): |10**5\eta|**(3/2) * (H/50)**3 * (Vpot/70)**3 * (1+0.1*Vshear)**(-2)
    print('GPI')
    dname = 'GPI'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        GPI = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        GPI = (1e5 * absolute(eta) )**(3/2) \
            * (H/50)**3 \
            * (PI.vmax/70)**3 \
            * (1+0.1*Vshear)**(-2)
        GPI.attrs['long_name'] = 'Genesis Potential Index'
        GPI.attrs['history'] = '|10**5\eta|**(3/2) * (RH/50)**3 * (Vpot/70)**3 * (1+0.1*Vshear)**(-2)'
        GPI.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)

    # GPI2010 (Emanuel 2010): |\eta|**3 * chi**(-4/3) * max((Vpot-35),0)**2 * (25+Vshear)**(-4)
    print('GPI2010')
    dname = 'GPI2010'
    ofile = os.path.join(odir, ibasename.replace('.nc', f'.{dname}.nc') )
    if os.path.exists(ofile):
        GPI2010 = xr.open_dataset(ofile)[dname]
        print('[opened]:', ofile)
    else:
        GPI2010 = absolute(eta)**3 \
            * chi_sb.where(chi_sb>0)**(-4/3) \
            * (PI.vmax - 35).clip(min=0)**2 \
            * (25 + Vshear)**(-4)
        GPI2010.attrs['long_name'] = 'Genesis Potential Index of Emanuel2010'
        GPI2010.attrs['history'] = '|\eta|**3 * chi**(-4/3) * max((Vpot-35),0)**2 * (25+Vshear)**(-4)'
        GPI2010.to_dataset(name=dname) \
            .to_netcdf(ofile, 
                encoding={dname: {'dtype': 'float32', 'zlib': True, 'complevel': 1}},
                unlimited_dims='time')
        print('[saved]:', ofile)
    
if __name__ == '__main__':
    odir = '/home/el2358/GEOCLIM/el2358/projects/data/era5/analysis/daily/TCI_v202201/'
    if len(sys.argv) == 3:
        # year and month provided
        year = int(sys.argv[1])
        month = int(sys.argv[2])
        do_tci(year, month, odir=odir)
    elif len(sys.argv) == 2:
        # only year provided -> loop through all months
        year = int(sys.argv[1])
        for month in range(1, 13):
            do_tci(year, month, odir=odir)
    else:
        # default: run 2020 Jul–Nov
        year = 2020
        for month in range(7, 12):
            do_tci(year, month, odir=odir)

    tt.check('**done**')

