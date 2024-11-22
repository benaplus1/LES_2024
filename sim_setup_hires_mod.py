import xarray as xr
import numpy as np
# %matplotlib inline
from datetime import datetime, timedelta
import h5py
import os
import shutil
from glob import glob
import xesmf
import netCDF4 as nc

def mod_ndvi(controlsrfdir, nlhsrfdir):
    #This function fills in NDVI over Lake Huron and Georgian Bay in the NLH simulation, because by default NDVI is set to zero there.
    inndvpaths = sorted(glob(f"{controlsrfdir}/ndh-N*"))
    outndvpaths = sorted(glob(f"{nlhsrfdir}/ndh-N*"))
    srffile_control = xr.open_dataset(f"{controlsrfdir}/sfch-S-g1.h5", engine = "h5netcdf", phony_dims = "sort")
    for i, inndvpath in enumerate(inndvpaths):
        print(f"Modifying NDVI File {inndvpath}")
        ndvds = xr.open_dataset(inndvpath, engine = "h5netcdf", phony_dims = "sort")
        newndvi_lh = ndvds["VEG_NDVIF"][1,:,:].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 0.425)[lhymin:lhymax, lhxmin:lhxmax] #Set constant NDVI over water
        newndvi_gb = ndvds["VEG_NDVIF"][1,:,:].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 0.425)[gbymin:gbymax, gbxmin:gbxmax] #Set constant NDVI over water
        #Set a  NDVI value of 0.425 (mean of NDVI of mixed forest over Lake Huron and Georgian Bay boxes) over areas within the Lake Huron and Georgian Bay boxes which are covered by WATER (that's the (patch_area > 0.1 business). Leave NDVI over land as is)
        newndvi2_lh = ndvds["VEG_NDVIF"][2,:,:].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 0.05)[lhymin:lhymax, lhxmin:lhxmax] #Set constant NDVI over water
        newndvi2_gb = ndvds["VEG_NDVIF"][2,:,:].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 0.05)[gbymin:gbymax, gbxmin:gbxmax] #Set constant NDVI over water
        #This assigns an NDVI to patch 2 of 0.05. This is a placeholder value. THIS HAS NO EFFECT ON THE SIMULATION, SINCE PATCH 2 IS SET TO ZERO AREA OVER WATER
        ndvds.close(); del ndvds
        #Now, we're going to actually write the the new NDVI data to the ndh-N files
        ndvfile = h5py.File(outndvpaths[i], "r+")
        ndvfile["VEG_NDVIF"][1,gbymin:gbymax, gbxmin:gbxmax] = np.copy(newndvi_gb.values) #Apply constant NDVI value to Patch 1 over Georgian Bay only
        ndvfile["VEG_NDVIF"][1,lhymin:lhymax, lhxmin:lhxmax] = np.copy(newndvi_lh.values) #Apply constant NDVI value to Patch 1 over Lake Huron only
        ndvfile["VEG_NDVIF"][2,lhymin:lhymax, lhxmin:lhxmax] = np.copy(newndvi2_lh.values) #Apply constant NDVI value to Patch 2 over Lake Huron only
        ndvfile["VEG_NDVIF"][2,gbymin:gbymax, gbxmin:gbxmax] = np.copy(newndvi2_gb.values) #Apply constant NDVI value to Patch 2 over Georgian Bay only
        ndvfile.close()
        del ndvfile
    srffile_control.close(); del srffile_control
    return None

def mod_vfiles(controlsrfdir, controlvardir, nlhvardir):
    #This function modifies soil moisture, soil temperature, snow water equivalent, and snow depth in the varfiles over Lake Huron and Georgian Bay in the NLH simulation
    srffile_control = xr.open_dataset(f"{controlsrfdir}/sfch-S-g1.h5", engine = "h5netcdf", phony_dims = "sort")
    invarpaths = sorted(glob(f"{controlvardir}/var-V*.h5"))
    outvarpaths = sorted(glob(f"{nlhvardir}/var-V*.h5"))
    for i,invarpath in enumerate(invarpaths):
        print(f"Modifying varfile {invarpath}")
        varout_nolake = xr.open_dataset(invarpath, engine = "h5netcdf", phony_dims = "sort")
        newsoilmoist1_lh = varout_nolake["SOILMOIST1"].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 0.3)[lhymin:lhymax, lhxmin:lhxmax] #Set constant soil moisture of 0.3 over water over Lake Huron
        newsoilmoist1_gb = varout_nolake["SOILMOIST1"].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 0.3)[gbymin:gbymax, gbxmin:gbxmax] #Set constant soil moisture of 0.3 over water over Georgian Bay
        
        newsoilmoist2_lh = varout_nolake["SOILMOIST2"].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 0.3)[lhymin:lhymax, lhxmin:lhxmax] #Set constant soil moisture over water over Lake Huron
        newsoilmoist2_gb = varout_nolake["SOILMOIST2"].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 0.3)[gbymin:gbymax, gbxmin:gbxmax] #Set constant soil moisture over water over Georgian Bay

        newsoiltemp1_lh = varout_nolake["SOILTEMP1"].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 272.25)[lhymin:lhymax, lhxmin:lhxmax] #Set constant soil temperature over water over Lake Huron
        newsoiltemp1_lh = newsoiltemp1_lh.where(newsoiltemp1_lh<=272.24, other = 272.25)
        newsoiltemp1_gb = varout_nolake["SOILTEMP1"].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 272.25)[gbymin:gbymax, gbxmin:gbxmax] #Set constant soil temperature over water over Georgian Bay
        newsoiltemp1_gb = newsoiltemp1_gb.where(newsoiltemp1_gb<=272.24, other = 272.25)

        newsoiltemp2_lh = varout_nolake["SOILTEMP2"].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 271.5)[lhymin:lhymax, lhxmin:lhxmax] #Set constant soil temperature over water over Lake Huron
        newsoiltemp2_lh = newsoiltemp2_lh.where(newsoiltemp2_lh<272.25, other = 271.5)
        newsoiltemp2_gb = varout_nolake["SOILTEMP2"].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 271.5)[gbymin:gbymax, gbxmin:gbxmax] #Set constant soil temperature over water over Georgian Bay
        newsoiltemp2_gb = newsoiltemp2_gb.where(newsoiltemp2_gb<272.25, other = 271.5)

        newsnowmass_lh = varout_nolake["SNOWMASS"].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 25)[lhymin:lhymax, lhxmin:lhxmax] #Set constant swe over water
        newsnowmass_lh = newsnowmass_lh.where(newsnowmass_lh>=25, other = 25)
        newsnowmass_gb = varout_nolake["SNOWMASS"].where((srffile_control["PATCH_AREA"][1]+srffile_control["PATCH_AREA"][2]).data>0.1, other = 25)[gbymin:gbymax, gbxmin:gbxmax]
        newsnowmass_gb = newsnowmass_gb.where(newsnowmass_gb>=25, other = 25)
        newsnowdepth_lh = newsnowmass_lh/100 #Set constant snow depth over water
        newsnowdepth_gb = newsnowmass_gb/100
        
        varout_nolake.close(); del varout_nolake
        hfile = h5py.File(outvarpaths[i], "r+")
        hfile["SOILMOIST1"][lhymin:lhymax, lhxmin:lhxmax] = np.copy(newsoilmoist1_lh.values) #Apply constant soil moisture to Lake Huron only
        hfile["SOILMOIST1"][gbymin:gbymax, gbxmin:gbxmax] = np.copy(newsoilmoist1_gb.values) #Apply constant soil moisture to Georgian Bay only
        hfile["SOILMOIST2"][lhymin:lhymax, lhxmin:lhxmax] = np.copy(newsoilmoist2_lh.values) #Apply constant soil moisture to Lake Huron only
        hfile["SOILMOIST2"][gbymin:gbymax, gbxmin:gbxmax] = np.copy(newsoilmoist2_gb.values) #Apply constant soil moisture to Georgian Bay only
        hfile["SOILTEMP1"][lhymin:lhymax, lhxmin:lhxmax] = np.copy(newsoiltemp1_lh.values) #Apply constant soil temperature to Lake Huron only
        hfile["SOILTEMP1"][gbymin:gbymax, gbxmin:gbxmax] = np.copy(newsoiltemp1_gb.values) #Apply constant soil temperature to Georgian Bay only
        hfile["SOILTEMP2"][lhymin:lhymax, lhxmin:lhxmax] = np.copy(newsoiltemp2_lh.values) #Apply constant soil temperature to Lake Huron only
        hfile["SOILTEMP2"][gbymin:gbymax, gbxmin:gbxmax] = np.copy(newsoiltemp2_gb.values) #Apply constant soil temperature to Georgian Bay only
        hfile["SNOWMASS"][lhymin:lhymax, lhxmin:lhxmax] = np.copy(newsnowmass_lh.values) #Apply constant swe to Lake Huron only
        hfile["SNOWMASS"][gbymin:gbymax, gbxmin:gbxmax] = np.copy(newsnowmass_gb.values) #Apply constant swe to Georgian Bay only
        hfile["SNOWDEPTH"][lhymin:lhymax, lhxmin:lhxmax] = np.copy(newsnowdepth_lh.values) #Apply constant snow depth to Lake Huron only
        hfile["SNOWDEPTH"][gbymin:gbymax, gbxmin:gbxmax] = np.copy(newsnowdepth_gb.values) #Apply constant snow depth to Lake Huron only
        
        hfile.close()
        del hfile
    srffile_control.close(); del srffile_control
    return None

def mod_vegtype(controlsrfdir, nlhsrfdir):
    #Modify vegetation type, soil texture, and patch areas over Lake Huron and Georgian Bay in the NLH simulation
    sds = xr.open_dataset(f"{controlsrfdir}/sfch-S-g1.h5", engine = "h5netcdf", phony_dims = "sort")
    newleaf1_lh = sds["LEAF_CLASS"][1,:,:].where((sds["PATCH_AREA"][1]+sds["PATCH_AREA"][2]).data>0.1, other = 14)[lhymin:lhymax, lhxmin:lhxmax] #Set constant vegetation class over water over Lake Huron
    newleaf1_gb = sds["LEAF_CLASS"][1,:,:].where((sds["PATCH_AREA"][1]+sds["PATCH_AREA"][2]).data>0.1, other = 14)[gbymin:gbymax, gbxmin:gbxmax] #Set constant vegetation class over water over Georgian Bay
    newleaf2_lh = sds["LEAF_CLASS"][2,:,:].where((sds["PATCH_AREA"][1]+sds["PATCH_AREA"][2]).data>0.1, other = 3)[lhymin:lhymax, lhxmin:lhxmax] #Placeholder for Patch 2
    newleaf2_gb = sds["LEAF_CLASS"][2,:,:].where((sds["PATCH_AREA"][1]+sds["PATCH_AREA"][2]).data>0.1, other = 3)[gbymin:gbymax, gbxmin:gbxmax] #Placeholder for Patch 2
    newpatch0_lh = np.where(sds["PATCH_AREA"][0].data>0.10, 0, sds["PATCH_AREA"][0])[lhymin:lhymax, lhxmin:lhxmax] #Very important! Modifies the patch area over Lake Huron to set all water patches to zero area
    newpatch0_gb = np.where(sds["PATCH_AREA"][0].data>0.10, 0, sds["PATCH_AREA"][0])[gbymin:gbymax, gbxmin:gbxmax] #very important! Modifies the patch areas over Georgian Bay to set all water patches to zero area
    newpatch1_lh = np.where(sds["PATCH_AREA"][0].data>0.10, 1, sds["PATCH_AREA"][1])[lhymin:lhymax, lhxmin:lhxmax] #Sets all areas which were water patches to land
    newpatch1_gb = np.where(sds["PATCH_AREA"][0].data>0.10, 1, sds["PATCH_AREA"][1])[gbymin:gbymax, gbxmin:gbxmax] #Sets all areas which were water patches to land
    newpatch2_lh = np.where(sds["PATCH_AREA"][0].data>0.10, 0, sds["PATCH_AREA"][2])[lhymin:lhymax, lhxmin:lhxmax] #Placeholder for Patch 2
    newpatch2_gb = np.where(sds["PATCH_AREA"][0].data>0.10, 0, sds["PATCH_AREA"][2])[gbymin:gbymax, gbxmin:gbxmax] #Placeholder for Patch 2
    newsoil0_lh = 4*np.ones((11,lhymax-lhymin,lhxmax-lhxmin)) #Set constant soil type to use under Lake Huron
    newsoil0_gb = 4*np.ones((11,gbymax-gbymin,gbxmax-gbxmin)) #Set constant soil type to use under Georgian Bay
    newsoil1_lh = sds["SOIL_TEXT"][1,:,:,:].where((sds["PATCH_AREA"][1]+sds["PATCH_AREA"][2]).data>0.1, other = 4)[:, lhymin:lhymax, lhxmin:lhxmax] #Apply this soil type to areas which were previously covered by water under Lake Huron
    newsoil1_gb = sds["SOIL_TEXT"][1,:,:,:].where((sds["PATCH_AREA"][1]+sds["PATCH_AREA"][2]).data>0.1, other = 4)[:, gbymin:gbymax, gbxmin:gbxmax] #Apply this soil type to areas which were previously covered by water under Georgian Bay
    newsoil2_lh = sds["SOIL_TEXT"][2,:,:,:].where((sds["PATCH_AREA"][1]+sds["PATCH_AREA"][2]).data>0.1, other = 4)[:, lhymin:lhymax, lhxmin:lhxmax] #Placeholder for Patch 2
    newsoil2_gb = sds["SOIL_TEXT"][2,:,:,:].where((sds["PATCH_AREA"][1]+sds["PATCH_AREA"][2]).data>0.1, other = 4)[:, gbymin:gbymax, gbxmin:gbxmax] #Placeholder for Patch 2
    sds.close(); del sds

    sfile = h5py.File(f"{nlhsrfdir}/sfch-S-g1.h5", "r+")
    # print(sfile)
    sfile["PATCH_AREA"][0, lhymin:lhymax, lhxmin:lhxmax] = np.copy(newpatch0_lh)
    sfile["PATCH_AREA"][0, gbymin:gbymax, gbxmin:gbxmax] = np.copy(newpatch0_gb)
    sfile["PATCH_AREA"][1, lhymin:lhymax, lhxmin:lhxmax] = np.copy(newpatch1_lh)
    sfile["PATCH_AREA"][1, gbymin:gbymax, gbxmin:gbxmax] = np.copy(newpatch1_gb)
    sfile["PATCH_AREA"][2, lhymin:lhymax, lhxmin:lhxmax] = np.copy(newpatch2_lh)
    sfile["PATCH_AREA"][2, gbymin:gbymax, gbxmin:gbxmax] = np.copy(newpatch2_gb)
    sfile["LEAF_CLASS"][0, lhymin:lhymax, lhxmin:lhxmax] = np.zeros((lhymax-lhymin, lhxmax-lhxmin))
    sfile["LEAF_CLASS"][0, gbymin:gbymax, gbxmin:gbxmax] = np.zeros((gbymax-gbymin, gbxmax-gbxmin))
    sfile["LEAF_CLASS"][1, lhymin:lhymax, lhxmin:lhxmax] = np.copy(newleaf1_lh.values)
    sfile["LEAF_CLASS"][1, gbymin:gbymax, gbxmin:gbxmax] = np.copy(newleaf1_gb.values)
    sfile["LEAF_CLASS"][2, lhymin:lhymax, lhxmin:lhxmax] = np.copy(newleaf2_lh.values)
    sfile["LEAF_CLASS"][2, gbymin:gbymax, gbxmin:gbxmax] = np.copy(newleaf2_gb.values)
    sfile["SOIL_TEXT"][0, :, lhymin:lhymax, lhxmin:lhxmax] = np.copy(newsoil0_lh)
    sfile["SOIL_TEXT"][0, :, gbymin:gbymax, gbxmin:gbxmax] = np.copy(newsoil0_gb)
    sfile["SOIL_TEXT"][1, :, lhymin:lhymax, lhxmin:lhxmax] = np.copy(newsoil1_lh.values)
    sfile["SOIL_TEXT"][1, :, gbymin:gbymax, gbxmin:gbxmax] = np.copy(newsoil1_gb.values)
    sfile["SOIL_TEXT"][2, :, lhymin:lhymax, lhxmin:lhxmax] = np.copy(newsoil2_lh.values)
    sfile["SOIL_TEXT"][2, :, gbymin:gbymax, gbxmin:gbxmax] = np.copy(newsoil2_gb.values)
    sfile.close()
    del sfile
    return None

def mod_sst(controlsrfdir, vartempsrfdir, glseapath):
    latlons = xr.open_dataset("glatglon.nc", engine = "h5netcdf", phony_dims = "sort")[["GLAT", "GLON"]] 
    '''
    Pull lat and lons for the simulation from the 'glatglon' file. This file was made by the authors ahead of time
    because the version of RAMS used in this study only outputs latitude and longitude fields in the actual 'output' files,
    not the mksfc files. Therefore, this routine could not be run until one of the simulations has been run for a short
    time and produced some output. This was felt to be an unnecessary burden for others, so the latitudes and longitudes
    for all simulations have been saved so that this script can be run immediately after MKSFC and MKVFILE have been completed,
    after which you can run all your simulations with RUNTYPE=INITIAL.
    '''
    sstfile_dec = xr.open_dataset(f"{controlsrfdir}/ssth-W-2009-12-27-000000-g1.h5", engine = "h5netcdf", phony_dims = "sort")["SEATF"] #SST file 1
    sstfile_jan = xr.open_dataset(f"{controlsrfdir}/ssth-W-2010-01-03-000000-g1.h5", engine = "h5netcdf", phony_dims = "sort")["SEATF"] #SST file 2
    glsea_dec = nc.Dataset(f"{glseapath}/2009_361_glsea_sst.nc", "r") #SST data from GLERL for 27 December 2009
    glsea_jan = nc.Dataset(f"{glseapath}/2010_003_glsea_sst.nc", "r") #SST data from GLERL for 3 January 2010
    gllats = glsea_dec["lat"]
    gllons = glsea_dec["lon"]
    decsstarray = np.asarray(glsea_dec["sst"][0,:,:])
    jansstarray = np.asarray(glsea_jan["sst"][0,:,:])
    glsea_dec_xr = xr.Dataset({"sst": (["lat", "lon"], np.where(decsstarray>-99998, decsstarray+273.15, np.nan))}, coords = {"lat": gllats[:], "lon": gllons[:]}) #Mask missing data in GLSEA
    glsea_jan_xr = xr.Dataset({"sst": (("lat", "lon"), np.where(jansstarray>-99998, jansstarray+273.15, np.nan))}, coords = {"lat": gllats[:], "lon": gllons[:]}) #Mask missing data in GLSEA
    sstfile_dec_xr = xr.Dataset({"SEATF": (["y", "x"], sstfile_dec.values)}, coords = {"lat": (["y", "x"], latlons["GLAT"].values), "lon": (["y", "x"], latlons["GLON"].values)}) #Assign lat/lon coordinates to RAMS SST file for regridding
    rgrid = xesmf.Regridder(glsea_dec_xr, sstfile_dec_xr, "bilinear", unmapped_to_nan = True) #Use XESMF to make a regrid map from GLSEA grid to RAMS grid
    interpglsea_dec = rgrid(glsea_dec_xr["sst"]) #Regrid GLSEA data from 27 December 2009 to the RAMS grid
    interpglsea_jan = rgrid(glsea_jan_xr["sst"]) #Regrid GLSEA data from 3 January 2010 to the RAMS grid
    sstfile_dec.close(); del sstfile_dec
    sstfile_jan.close(); del sstfile_jan
    newsst_dec = h5py.File(f"{vartempsrfdir}/ssth-W-2009-12-27-000000-g1.h5", "r+") #Open the SST file in the VARTEMP directory
    newsst_jan = h5py.File(f"{vartempsrfdir}/ssth-W-2010-01-03-000000-g1.h5", "r+")
    nonansst_dec = np.where(np.isnan(interpglsea_dec.values), newsst_dec["SEATF"][:,:], interpglsea_dec.values)
    '''
    GLSEA data only covers the Great Lakes, with missing data elsewhere. Since RAMS expects the SST file to cover the entire grid,
    we'll fill in the missing data (everywhere but the Great Lakes) with the values which were already in the SST files
    '''
    nonansst_jan = np.where(np.isnan(interpglsea_jan.values), newsst_jan["SEATF"][:,:], interpglsea_dec.values)
    newsst_dec["SEATF"][:,:] = np.copy(nonansst_dec)
    newsst_jan["SEATF"][:,:] = np.copy(nonansst_jan)
    newsst_dec.close()
    newsst_jan.close()
    return None

def mod_snow(controlsrfdir, controlvardir, moresnowvardir):
    invarpaths = sorted(glob(f"{controlvardir}/var-V*.h5")) #Varfile paths in the CONTROl directory
    outvarpaths = sorted(glob(f"{moresnowvardir}/var-V*.h5")) #Varfile paths in the MORESNOW directory
    patcharea = xr.open_dataset(f"{controlsrfdir}/sfch-S-g1.h5")["PATCH_AREA"].values #Patch area in the CONTROL simulation
    for i, invarpath in enumerate(invarpaths):
        print(f"Modifying varfile {invarpath}")
        controlvfile = xr.open_dataset(invarpath, engine = "h5netcdf", phony_dims = "sort")
        snowmassfield = controlvfile["SNOWMASS"][:,:].values #CONTROL varfile snow mass
        snowdepthfield = controlvfile["SNOWDEPTH"][:,:].values #CONTROL varvfile snow depth
        snowmasscopy = np.copy(snowmassfield)
        snowdepthcopy = np.copy(snowdepthfield)
        snowmassmean = np.nanmean(np.where(patcharea[1,650:850,375:425]==1, 1, np.nan)*snowmasscopy[650:850,375:425])
        #Mean snow mass over the western Ontario Peninsula
        snowdepthmean = np.nanmean(np.where(patcharea[1,650:850,375:425]==1, 1, np.nan)*snowdepthcopy[650:850,375:425])
        #Ditto for snow depth
        snowmasscopy[650:850,375:600] = np.where(patcharea[1,650:850,375:600]==1, snowmassmean, snowmasscopy[650:850,375:600])
        '''
        Assign the value of mean snow mass over the western Ontario Peninsula to the ENTIRE Ontario Peninsula.
        The reason we do this is that the ERA5 snow data for some reason has an area of very little snow cover over
        the eastern Ontario Peninsula, which leads to the formation of Brown Lake-Effect bands. By filling in this
        area of little snow cover, we can eliminate the presence of these Brown Lake-Effect bands by reducing the
        sensible heat flux over the Ontario Peninsula.
        '''
        snowdepthcopy[650:850,375:600] = np.where(patcharea[1,650:850,375:600]==1, snowdepthmean, snowdepthcopy[650:850,375:600])
        #Ditto for snow depth
        controlvfile.close()
        del controlvfile
        moresnowvfile = h5py.File(outvarpaths[i], "r+") #Open the varfiles in the MORESNOW directory
        moresnowvfile["SNOWMASS"][:,:] = snowmasscopy #Replace the snow mass field with the one with more snow mass over the Ontario Peninsula
        moresnowvfile["SNOWDEPTH"][:,:] = snowdepthcopy #Ditto for snow depth
        moresnowvfile.close(); del moresnowvfile
    return None


# In[14]:
srfmade = input("Have you run RAMS with the RUNTYPE in the RAMSIN.les_control_hires set to 'MKSFC?'? yes or no: ")
if srfmade.strip().lower() == "no":
    raise Exception("MKSFC must be run before surface modifications can be performed!")
controlsrfdir = input("Enter the path of the folder containing mksfc files \n(sfch-S-g1.h5, toph-S-g1.h5, ndh-N files, ssth-W files): ")
srffolderexists = os.path.exists(controlsrfdir)
if not srffolderexists:
    raise Exception("Surface file folder not found! Aborting!")
else:
    controlsrfdir = controlsrfdir.rstrip("/")
    if not os.path.exists(f"{controlsrfdir}/sfch-S-g1.h5"):
        raise Exception("Surface file not found! Aborting!")
    elif not os.path.exists(f"{controlsrfdir}/toph-S-g1.h5"):
        raise Exception("Topography file not found! We don't need it here, but that means something went very wrong with makesfc! Aborting!")
    elif not all([os.path.exists(f"{controlsrfdir}/{entry}") for entry in ["ssth-W-2009-12-27-000000-g1.h5", "ssth-W-2010-01-03-000000-g1.h5", "ssth-W-2010-01-10-000000-g1.h5"]]):
        raise Exception("Sea Surface Temperature files not all found! Aborting!")
    elif not all([os.path.exists(f"{controlsrfdir}/ndh-N-0000-{ndvmonth}-g1.h5") for ndvmonth in ["01-16-120000", "02-15-000000", "03-16-120000", "04-16-000000", "05-16-120000", "06-16-000000", "07-16-120000", "08-16-120000", "09-16-000000", "10-16-120000", "11-16-000000", "12-16-120000"]]):
        raise Exception("NDVI files not all found! Aborting!")
varsmade = input("Have you run RAMS with the RUNTYPE in the RAMSIN.les_control_hires set to 'MKVFILE'? yes or no: ")
if varsmade.strip().lower() == "no":
    raise Exception("MKVFILE must be run before surface modifications can be performed!")
controlvardir = input("Enter the path of the folder containing all varfiles: ")
if not os.path.exists(controlvardir):
    raise Exception("Varfiles folder not found! Aborting!")
else:
    if len(glob(f"{controlvardir}/var-V*.h5")) == 0:
        raise Exception("No varfiles found within specified directory! Aborting!")

print("For this routine, we're going to modify the surface, ndvi, and ssth files produced by mksfc, as well as the varfiles produced by mkvfile in RAMSIN.les_control_hires,",
      "to use in the NLH, VARTEMP, and MORESNOW simulations. To ensure they're consistent with RAMS format, we'll modify them with h5py's 'r+' routine, which requires that the files to modify already",
      "exist. Therefore, we'll check if the mksfc files and varfiles already exist in the NLH, VARTEMP, and MORESNOW directories, and if not, copy them over.")
nlhsrfdir = input("Enter the directory for the mksfc files for the NLH simulation: ").rstrip("/")
nlhvardir = input("Enter the directory for the varfiles for the NLH simulation: ").rstrip("/")
vartempsrfdir = input("Enter the directory for the mksfc files for the VARTEMP simulation: ").rstrip("/")
vartempvardir = input("Enter the directory for the varfiles for the VARTEMP simulation: ").rstrip("")
moresnowsrfdir = input("Enter the directory for the mksfc files for the MORESNOW simulation: ").rstrip("/")
moresnowvardir = input("Enter the directory for the varfiles for the MORESNOW simulation: ").rstrip("/")
if not os.path.exists(nlhsrfdir):
    raise Exception("The provided mksfc file directory for the NLH simulation doesn't exist! Aborting!")
elif not os.path.exists(nlhvardir):
    raise Exception("The provided varfile directory for the NLH simulation doesn't exist! Aborting!")
if os.path.exists(f"{nlhsrfdir}/sfch-S-g1.h5"):
    print("Surface file already copied to NLH directory! Skipping copying!")
else:
    shutil.copy(f"{controlsrfdir}/sfch-S-g1.h5", f"{nlhsrfdir}/sfch-S-g1.h5")
    print("sfch-S-g1.h5 copied to NLH directory!")
if os.path.exists(f"{nlhsrfdir}/toph-S-g1.h5"):
    print("Topography file already copied to NLH directory! Skipping copying!")
else:
    shutil.copy(f"{controlsrfdir}/toph-S-g1.h5", f"{nlhsrfdir}/toph-S-g1.h5")
    print("toph-S-g1.h5 copied to NLH directory!")
if all([os.path.exists(f"{nlhsrfdir}/{entry}") for entry in ["ssth-W-2009-12-27-000000-g1.h5", "ssth-W-2010-01-03-000000-g1.h5", "ssth-W-2010-01-10-000000-g1.h5"]]):
    print("All SST files already copied to NLH directory! Skipping copying!")
else:
    for sstfile in sorted(glob(f"{controlsrfdir}/ssth-W*")):
        shutil.copy(sstfile, f"{nlhsrfdir}/")
    print("ssth-W files copied to NLH directory!")
if all([os.path.exists(f"{nlhsrfdir}/ndh-N-0000-{ndvmonth}-g1.h5") for ndvmonth in ["01-16-120000", "02-15-000000", "03-16-120000", "04-16-000000", "05-16-120000", "06-16-000000", "07-16-120000", "08-16-120000", "09-16-000000", "10-16-120000", "11-16-000000", "12-16-120000"]]):
    print("All NDVI files already copied to NLH directory! Skipping copying!")
else:
    for ndhfile in sorted(glob(f"{controlsrfdir}/ndh-N*")):
        shutil.copy(ndhfile, f"{nlhsrfdir}/")
    print("ndh-N files copied to NLH directory!")
if all([os.path.exists(f"{nlhvardir}/{vfile[vfile.index('var-V'):]}") for vfile in sorted(glob(f"{controlvardir}/var-V*"))]):
    print("All varfiles already copied to NLH directory! Skipping copying!")
else:
    print("Copying varfiles (var-V files). This may take some time.")
    for vfile in sorted(glob(f"{controlvardir}/var-V*")):
        print(f"Copying varfile {vfile} to NLH directory")
        shutil.copy(vfile, f"{nlhvardir}/")
    print("Varfiles copied to NLH directory!")
# srffile_control = xr.open_dataset(srfpath, engine = "h5netcdf", phony_dims = "sort")
lhymin = 640; lhymax = 1040; lhxmin = 140; lhxmax = 420 #lh - main body of Lake Huron
gbymin = 800; gbymax = 1000; gbxmin = 420; gbxmax = 700 #gb - Georgian Bay and downwind land
lakebounds = {"lhymin": lhymin, "lhymax": lhymax, "lhxmin": lhxmin, "lhxmax": lhxmax,
              "gbymin": gbymin, "gbymax": gbymax, "gbxmin": gbxmin, "gbxmax": gbxmax}
#NOTE THAT THIS ONLY WORKS FOR THE CURRENT HORIZONTAL GRID SPACING OF 1km - IF THAT CHANGES, NEED TO REWORK THE INDICES HERE
print("Modifying sfch-S-g1.h5 for the NLH simulation")
mod_vegtype(controlsrfdir, nlhsrfdir)
print("Modifying NDVI for the NLH simulation")
mod_ndvi(controlsrfdir, nlhsrfdir)
print("Modifying varfiles for the NLH simulation")
mod_vfiles(controlsrfdir, controlvardir, nlhvardir)

glseadir = input("Enter the directory containing the GLSEA water temperature files: ")
if not all([os.path.exists(f"{glseadir}/2009_361_glsea_sst.nc"), os.path.exists(f"{glseadir}/2010_003_glsea_sst.nc")]):
    raise Exception("Not all GLSEA files found! Aborting!")
if not os.path.exists(vartempsrfdir):
    raise Exception("The provided mksfc file directory for the VARTEMP simulation doesn't exist! Aborting!")
elif not os.path.exists(vartempvardir):
    raise Exception("The provided varfile directory for the VARTEMP simulation doesn't exist! Aborting!")
if os.path.exists(f"{vartempsrfdir}/sfch-S-g1.h5"):
    print("Surface file already copied to VARTEMP directory! Skipping copying!")
else:
    shutil.copy(f"{controlsrfdir}/sfch-S-g1.h5", f"{vartempsrfdir}/sfch-S-g1.h5")
    print("sfch-S-g1.h5 copied to VARTEMP directory!")
if os.path.exists(f"{vartempsrfdir}/toph-S-g1.h5"):
    print("Topography file already copied to VARTEMP directory! Skipping copying!")
else:
    shutil.copy(f"{controlsrfdir}/toph-S-g1.h5", f"{vartempsrfdir}/toph-S-g1.h5")
    print("toph-S-g1.h5 copied to VARTEMP directory!")
if all([os.path.exists(f"{vartempsrfdir}/{entry}") for entry in ["ssth-W-2009-12-27-000000-g1.h5", "ssth-W-2010-01-03-000000-g1.h5", "ssth-W-2010-01-10-000000-g1.h5"]]):
    print("All SST files already copied to VARTEMP directory! Skipping copying!")
else:
    for sstfile in sorted(glob(f"{controlsrfdir}/ssth-W*")):
        shutil.copy(sstfile, f"{vartempsrfdir}/")
    print("ssth-W files copied to VARTEMP directory!")
if all([os.path.exists(f"{vartempsrfdir}/ndh-N-0000-{ndvmonth}-g1.h5") for ndvmonth in ["01-16-120000", "02-15-000000", "03-16-120000", "04-16-000000", "05-16-120000", "06-16-000000", "07-16-120000", "08-16-120000", "09-16-000000", "10-16-120000", "11-16-000000", "12-16-120000"]]):
    print("All NDVI files already copied to VARTEMP directory! Skipping copying!")
else:
    for ndhfile in sorted(glob(f"{controlsrfdir}/ndh-N*")):
        shutil.copy(ndhfile, f"{vartempsrfdir}/")
    print("ndh-N files copied to VARTEMP directory!")
if all([os.path.exists(f"{vartempvardir}/{vfile[vfile.index('var-V'):]}") for vfile in sorted(glob(f"{controlvardir}/var-V*"))]):
    print("All varfiles already copied to VARTEMP directory! Skipping copying!")
else:
    print("Copying varfiles (var-V files). This may take some time.")
    for vfile in sorted(glob(f"{controlvardir}/var-V*")):
        print(f"Copying varfile {vfile} to VARTEMP directory")
        shutil.copy(vfile, f"{vartempvardir}/")
    print("Varfiles copied to VARTEMP directory!")

print("Modifying SST files for the VARTEMP simulation")
mod_sst(controlsrfdir, vartempsrfdir, glseadir)

if not os.path.exists(moresnowsrfdir):
    raise Exception("The provided mksfc file directory for the MORESNOW simulation doesn't exist! Aborting!")
elif not os.path.exists(moresnowvardir):
    raise Exception("The provided varfile directory for the MORESNOW simulation doesn't exist! Aborting!")
if os.path.exists(f"{moresnowsrfdir}/sfch-S-g1.h5"):
    print("Surface file already copied to MORESNOW directory! Skipping copying!")
else:
    shutil.copy(f"{controlsrfdir}/sfch-S-g1.h5", f"{moresnowsrfdir}/sfch-S-g1.h5")
    print("sfch-S-g1.h5 copied to MORESNOW directory!")
if os.path.exists(f"{moresnowsrfdir}/toph-S-g1.h5"):
    print("Topography file already copied to MORESNOW directory! Skipping copying!")
else:
    shutil.copy(f"{controlsrfdir}/toph-S-g1.h5", f"{moresnowsrfdir}/toph-S-g1.h5")
    print("toph-S-g1.h5 copied to MORESNOW directory!")
if all([os.path.exists(f"{moresnowsrfdir}/{entry}") for entry in ["ssth-W-2009-12-27-000000-g1.h5", "ssth-W-2010-01-03-000000-g1.h5", "ssth-W-2010-01-10-000000-g1.h5"]]):
    print("All SST files already copied to MORESNOW directory! Skipping copying!")
else:
    for sstfile in sorted(glob(f"{controlsrfdir}/ssth-W*")):
        shutil.copy(sstfile, f"{moresnowsrfdir}/")
    print("ssth-W files copied to MORESNOW directory!")
if all([os.path.exists(f"{moresnowsrfdir}/ndh-N-0000-{ndvmonth}-g1.h5") for ndvmonth in ["01-16-120000", "02-15-000000", "03-16-120000", "04-16-000000", "05-16-120000", "06-16-000000", "07-16-120000", "08-16-120000", "09-16-000000", "10-16-120000", "11-16-000000", "12-16-120000"]]):
    print("All NDVI files already copied to MORESNOW directory! Skipping copying!")
else:
    for ndhfile in sorted(glob(f"{controlsrfdir}/ndh-N*")):
        shutil.copy(ndhfile, f"{moresnowsrfdir}/")
    print("ndh-N files copied to MORESNOW directory!")
if all([os.path.exists(f"{moresnowvardir}/{vfile[vfile.index('var-V'):]}") for vfile in sorted(glob(f"{controlvardir}/var-V*"))]):
    print("All varfiles already copied to MORESNOW directory! Skipping copying!")
else:
    print("Copying varfiles (var-V files). This may take some time.")
    for vfile in sorted(glob(f"{controlvardir}/var-V*")):
        print(f"Copying varfile {vfile} to MORESNOW directory")
        shutil.copy(vfile, f"{moresnowvardir}/")
    print("Varfiles copied to MORESNOW directory!")

print("Modifying varfiles for the MORESNOW simulation.")
mod_snow(controlsrfdir, controlvardir, moresnowvardir)
print("Surface File and Varfile Corrections Complete!")




