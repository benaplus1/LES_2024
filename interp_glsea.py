import xarray as xr
import numpy as np
import netCDF4 as nc
import xesmf
import h5py
import os
import shutil
import urllib

latlons = xr.open_dataset("glatglon.nc")
sstfile_dec = xr.open_dataset("/moonbow/ascheb/les/2010/hires_control/ssth-W-2009-12-27-000000-g1.h5", engine = "h5netcdf", phony_dims = "sort")["SEATF"]
sstfile_jan = xr.open_dataset("/moonbow/ascheb/les/2010/hires_control/ssth-W-2010-01-03-000000-g1.h5", engine = "h5netcdf", phony_dims = "sort")["SEATF"]

glsea_dec_url = "https://apps.glerl.noaa.gov/thredds/fileServer/glsea_nc/2009/12/2009_361_glsea_sst.nc"
glsea_dec_file = "/moonbow/ascheb/les/python/2009_361_glsea_sst.nc"
urllib.request.urlretrieve(glsea_dec_url, glsea_dec_file)
glsea_jan_url = "https://apps.glerl.noaa.gov/thredds/fileServer/glsea_nc/2010/01/2010_003_glsea_sst.nc"
glsea_jan_file = "/moonbow/ascheb/les/python/2010_003_glsea_sst.nc"
urllib.request.urlretrieve(glsea_jan_url, glsea_jan_file)

glsea_dec = nc.Dataset("/moonbow/ascheb/les/python/2009_361_glsea_sst.nc", "r")
glsea_jan = nc.Dataset("/moonbow/ascheb/les/python/2010_003_glsea_sst.nc", "r")

gllats = glsea_dec["lat"]
gllons = glsea_dec["lon"]
decsstarray = np.asarray(glsea_dec["sst"][0,:,:])
jansstarray = np.asarray(glsea_jan["sst"][0,:,:])

glsea_dec_xr = xr.Dataset({"sst": (["lat", "lon"], np.where(decsstarray>-99998, decsstarray+273.15, np.nan))}, coords = {"lat": gllats[:], "lon": gllons[:]})
glsea_jan_xr = xr.Dataset({"sst": (("lat", "lon"), np.where(jansstarray>-99998, jansstarray+273.15, np.nan))}, coords = {"lat": gllats[:], "lon": gllons[:]})

sstfile_dec_xr = xr.Dataset({"SEATF": (["y", "x"], sstfile_dec.values)}, coords = {"lat": (["y", "x"], latlons["GLAT"].values), "lon": (["y", "x"], latlons["GLON"].values)})
sstfile_jan_xr = xr.Dataset({"SEATF": (["y", "x"], sstfile_jan.values)}, coords = {"lat": (["y", "x"], latlons["GLAT"].values), "lon": (["y", "x"], latlons["GLON"].values)})

rgrid = xesmf.Regridder(glsea_dec_xr, sstfile_dec_xr, "bilinear", unmapped_to_nan = True)
rgrid.to_netcdf("glseaweights.nc")
interpglsea_dec = rgrid(glsea_dec_xr["sst"])
interpglsea_jan = rgrid(glsea_jan_xr["sst"])

sstfile_dec.close(); del sstfile_dec
sstfile_jan.close(); del sstfile_jan

if not os.path.exists("/moonbow/ascheb/les/2010/hires_vartemp"):
    os.mkdir("/moonbow/ascheb/les/2010/hires_vartemp")
shutil.copyfile("/moonbow/ascheb/les/2010/hires_control/ssth-W-2009-12-27-000000-g1.h5", "/moonbow/ascheb/les/2010/hires_vartemp/ssth-W-2009-12-27-000000-g1.h5")
shutil.copyfile("/moonbow/ascheb/les/2010/hires_control/ssth-W-2010-01-03-000000-g1.h5", "/moonbow/ascheb/les/2010/hires_vartemp/ssth-W-2010-01-03-000000-g1.h5")
newsst_dec = h5py.File("/moonbow/ascheb/les/2010/hires_vartemp/ssth-W-2009-12-27-000000-g1.h5", "r+")
newsst_jan = h5py.File("/moonbow/ascheb/les/2010/hires_vartemp/ssth-W-2010-01-03-000000-g1.h5", "r+")

nonansst_dec = np.where(np.isnan(interpglsea_dec.values), newsst_dec["SEATF"][:,:], interpglsea_dec.values)
nonansst_jan = np.where(np.isnan(interpglsea_jan.values), newsst_jan["SEATF"][:,:], interpglsea_dec.values)

newsst_dec["SEATF"][:,:] = np.copy(nonansst_dec)
newsst_jan["SEATF"][:,:] = np.copy(nonansst_jan)

newsst_dec.close()
newsst_jan.close()

