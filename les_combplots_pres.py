import xarray as xr
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.lines import Line2D
# import matplotlib as mpl
# from metpy.units import units
# from metpy.plots import SkewT
# import metpy
from pandas import read_csv, date_range
from time import perf_counter
from glob import glob
from matplotlib.style import use as usestyle
from importlib import reload
import makevids as mkv
from concurrent.futures import ProcessPoolExecutor
from os import path, mkdir
from functools import partial

from matplotlib import font_manager as fm
fontdir = "/home/ascheb/libfonts/*.ttf"
for fpath in glob(fontdir):
    # print(fpath)
    fm.fontManager.addfont(fpath)
    
usestyle("presplots.mplstyle")
from matplotlib import rcParams
rcParams["figure.titlesize"] = 40
rcParams["axes.titlesize"] = 35
rcParams["axes.labelsize"] = 30
rcParams["xtick.labelsize"] = 27
rcParams["ytick.labelsize"] = 27

def make_combplots(t):
    folderprepath = "/moonbow/ascheb/les/"
    ftime = t.strftime("%Y-%m-%d-%H%M%S")
    print(ftime)
    afile_control = xr.open_dataset(f"{folderprepath}2010/hires_control/processed_data/mergedvars_{ftime}.nc")
    afile_nolake = xr.open_dataset(f"{folderprepath}2010/hires_nolake/processed_data/mergedvars_{ftime}.nc")
    # afile_pasturebroadforest = xr.open_dataset(f"{folderpath}pasturebroadforest-wind/processed_data/mergedvars_{ftime}.nc")
    if not path.exists(f"{folderprepath}combplots"):
        mkdir(f"{folderprepath}combplots")
    
    fields = ["SnowRate", "ShFlux", "LhFlux", "VapMix", "CldTop", "w", "Div"]
    # fields = ["SnowRate", "w", "div"]
    # fields = ["ShFlux", "LhFlux", "CldTop"]
    # fields = ["w"]
    for field in fields:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (20, 14), dpi = 200, layout = "constrained")
        for ax in (ax1, ax2):
            ax.set_xlabel("Longitude (Degrees E)")
            ax.set_ylabel("Latitude (Degrees N)")
            ax.set_yticks([41, 42, 43, 44, 45, 46])
            ax.set_xticks([-84, -82, -80])
            # ax.pcolormesh(afile_control["lon1d"], afile_control["lat1d"], afile_control["Patch"][0,:,:].where(afile_control["Patch"][0,:,:]==1), color = "Navy", zorder = 0)
        ax1.set_title("CONTROL", fontfamily = "Liberation Serif", y = 1.02)
        ax2.set_title("NLH", fontfamily = "Liberation Serif", y = 1.02)
        
        print(f"Plotting Field {field}")
        if field == "SnowRate":
            for ax in fig.get_axes():
                ax.set_yticks([38, 40, 42, 44, 46])
                ax.set_xticks([-84, -82, -80, -78, -76])
            fig.suptitle(f"Snowfall Rate at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
            from palettable.cmocean.sequential import Ice_15
            icecmap = Ice_15.get_mpl_colormap()
            fakecontour = Line2D([], [], linestyle = "-", linewidth = 0.7, color = "white", label = "Water Bodies")
            snowmp = ax1.pcolormesh(afile_control["lon1d"], afile_control["lat1d"], afile_control["SnowPrecipRate"]+afile_control["AggPrecipRate"]+afile_control["PrisPrecipRate"], shading = "nearest", cmap = icecmap, vmin = 0, vmax = 4)
            ax2.pcolormesh(afile_control["lon1d"], afile_control["lat1d"], afile_nolake["SnowPrecipRate"]+afile_nolake["AggPrecipRate"]+afile_nolake["PrisPrecipRate"], shading = "nearest", cmap = icecmap, vmin = 0, vmax = 4)
            ax1.contour(afile_control["lon1d"], afile_control["lat1d"], afile_control["Patch"][1,:,:], levels = [0.999], colors = "white", linestyles = "-", linewidths = 0.7)
            ax1.legend(loc = "lower left", fontsize = 20, handles = [fakecontour])
            ax2.contour(afile_control["lon1d"], afile_control["lat1d"], afile_nolake["Patch"][1,:,:], levels = [0.999], colors = "white", linestyles = "-", linewidths = 0.7)
            ax2.legend(loc = "lower left", fontsize = 20, handles = [fakecontour])
            cbar = fig.colorbar(snowmp, ax = fig.get_axes(), orientation = "horizontal", fraction = 0.05, extend = "max"); cbar.set_label("Snowfall Rate (mm/hr)")
            fig.savefig(f"{folderprepath}combplots/snowrate_{str(t.day).zfill(2)}{str(t.hour).zfill(2)}{str(t.minute).zfill(2)}z.png")

        elif field == "w":
            walt = 700 #m AMSL, altitude at which to evalute w
            fig.suptitle(f"700m AMSL Vertical Velocity at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
            fakecontour = Line2D([], [], color = "navy", linewidth = 0.7, label = "Water Bodies")
            wmp = ax1.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["w"].sel(z = walt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "RdBu_r", vmin = -3, vmax = 3)
            ax1.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "navy")
            ax2.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["w"].sel(z = walt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "RdBu_r", vmin = -3, vmax = 3)
            ax2.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "navy")
            ax1.legend(loc = "lower left", fontsize = 20, handles = [fakecontour], labelcolor = "black")
            ax2.legend(loc = "lower left", fontsize = 20, handles = [fakecontour], labelcolor = "black")
            cbar = fig.colorbar(wmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label("700m AMSL Vertical Velocity (m/s)")
            fig.savefig(f"{folderprepath}combplots/wcomp_{str(t.day).zfill(2)}{str(t.hour).zfill(2)}{str(t.minute).zfill(2)}z.png")
        elif field == "Div":
            divalt = 350 #m AMSL, altitude at which to evaluate horizontal divergence.
            fakecontour = Line2D([], [], color = "navy", linewidth = 0.7, label = "Water Bodies")
            fig.suptitle(f"{divalt}m AMSL Horizonal Divergence at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
            horizdiv_control = np.gradient(afile_control["u"].sel(z = divalt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)).values, 1000, axis = 1)+np.gradient(afile_control["v"].sel(z = divalt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)).values, 1000, axis = 0)
            horizdiv_nolake = np.gradient(afile_nolake["u"].sel(z = divalt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)).values, 1000, axis = 1)+np.gradient(afile_nolake["v"].sel(z = divalt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)).values, 1000, axis = 0)
            divmp = ax1.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), horizdiv_control, cmap = "BrBG", vmin = -5*10**(-3), vmax = 5*10**(-3), zorder = 1)
            ax1.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "navy")
            ax2.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), horizdiv_nolake, cmap = "BrBG", vmin = -5*10**(-3), vmax = 5*10**(-3), zorder = 1)
            ax2.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "navy")
            ax1.legend(loc = "lower left", fontsize = 20, handles = [fakecontour], labelcolor = "black")
            ax2.legend(loc = "lower left", fontsize = 20, handles = [fakecontour], labelcolor = "black")
            cbar = fig.colorbar(divmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"{divalt}m AMSL Horizontal Divergence ($s^{-1}$)")
            fig.savefig(f"{folderprepath}combplots/divcomp_{str(t.day).zfill(2)}{str(t.hour).zfill(2)}{str(t.minute).zfill(2)}z.png")
            
        elif field == "VapMix":
            vapalt = 700 #m AMSL, altitude at which to evaluate horizontal divergence.
            fakecontour = Line2D([], [], color = "navy", linewidth = 0.7, label = "Water Bodies")
            fig.suptitle(f"{vapalt}m AMSL Vapor Mixing Ratio at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
            vapmp = ax1.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["VaporMix"].sel(z = vapalt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = 0, vmax = 2, zorder = 1)
            ax1.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "navy")
            ax2.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["VaporMix"].sel(z = vapalt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = 0, vmax = 2, zorder = 1)
            ax2.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "navy")
            ax1.legend(loc = "lower left", fontsize = 20, handles = [fakecontour], labelcolor = "black")
            ax2.legend(loc = "lower left", fontsize = 20, handles = [fakecontour], labelcolor = "black")
            cbar = fig.colorbar(vapmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"{vapalt}m AMSL Vapor Mixing Ratio (g/kg)")
            fig.savefig(f"{folderprepath}combplots/vapcomp_{str(t.day).zfill(2)}{str(t.hour).zfill(2)}{str(t.minute).zfill(2)}z.png")
            
        elif field == "ShFlux":
            fakecontour = Line2D([], [], color = "black", linewidth = 0.7, label = "Water Bodies")
            fig.suptitle(f"Surface Sensible Heat Flux at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
            shmp = ax1.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["SensibleHeatFlux"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "bwr", vmin = -500, vmax = 500, zorder = 1)
            ax1.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black")
            ax2.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["SensibleHeatFlux"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "bwr", vmin = -500, vmax = 500, zorder = 1)
            ax2.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black")
            ax1.legend(loc = "lower left", fontsize = 20, handles = [fakecontour], labelcolor = "black")
            ax2.legend(loc = "lower left", fontsize = 20, handles = [fakecontour], labelcolor = "black")
            cbar = fig.colorbar(shmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"Sensible Heat Flux (w/m^2)")
            fig.savefig(f"{folderprepath}combplots/shcomp_{str(t.day).zfill(2)}{str(t.hour).zfill(2)}{str(t.minute).zfill(2)}z.png")
            
        elif field == "LhFlux":
            fakecontour = Line2D([], [], color = "black", linewidth = 0.7, label = "Water Bodies")
            fig.suptitle(f"Surface Latent Heat Flux at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
            lhmp = ax1.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["LatentHeatFlux"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = -500, vmax = 500, zorder = 1)
            ax1.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black")
            ax2.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["LatentHeatFlux"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = -500, vmax = 500, zorder = 1)
            ax2.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black")
            ax1.legend(loc = "lower left", fontsize = 20, handles = [fakecontour], labelcolor = "black")
            ax2.legend(loc = "lower left", fontsize = 20, handles = [fakecontour], labelcolor = "black")
            cbar = fig.colorbar(lhmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"Latent Heat Flux (w/m^2)")
            fig.savefig(f"{folderprepath}combplots/lhcomp_{str(t.day).zfill(2)}{str(t.hour).zfill(2)}{str(t.minute).zfill(2)}z.png")
            
        elif field == "CldTop":
            fakecontour = Line2D([], [], color = "white", linewidth = 0.7, label = "Water Bodies")
            fig.suptitle(f"Cloud Top Height at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
            cldtopmp = ax1.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["CloudTop"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "gist_earth", vmin = 0, vmax = 3000, zorder = 1)
            ax1.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "white")
            ax2.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["CloudTop"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "gist_earth", vmin = 0, vmax = 3000, zorder = 1)
            ax2.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "white")
            ax1.legend(loc = "lower left", fontsize = 20, handles = [fakecontour], labelcolor = "white")
            ax2.legend(loc = "lower left", fontsize = 20, handles = [fakecontour], labelcolor = "white")
            cbar = fig.colorbar(cldtopmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"Cloud Top Height (m AMSL)")
            fig.savefig(f"{folderprepath}combplots/cldtop_{str(t.day).zfill(2)}{str(t.hour).zfill(2)}{str(t.minute).zfill(2)}z.png")
            
        plt.close(); del fig; del ax1; del ax2; del cbar

    afile_nolake.close(); del afile_nolake
    afile_control.close(); del afile_control
    return ftime
runtype = input("Do you want to *plot*, *animate* existing plots, or do *both*? ")
if runtype.lower() == "plot":
    st = perf_counter()
    t0 = input("Enter the start time for plotting in yyyy-mm-dd-HHMMSS format: ")
    tf = input("Enter the end time for plotting in yyyy-mm-dd-HHMMSS format: ")
    tlist = date_range(t0, tf, freq = timedelta(minutes = 20)).to_pydatetime()
    # for t in tlist:
    #     make_combplots(t)
    ppool = ProcessPoolExecutor(max_workers = 4) 
    ppool.map(make_combplots, tlist)
    ppool.shutdown()
    et = perf_counter()
    print(f"Plotting took {et-st:.2f} seconds")
    
elif runtype.lower() == "animate":
    st = perf_counter()
    folderprepath = "/moonbow/ascheb/les/"
    reload(mkv)
    # fields = ["bowen"]
    # fields = ["srftemp", "rainrate", "srfpres", "vapplan", "cldtop", "uwind", "lhflux", "shflux", "bowen", "vertintcond", "vertintice", "vertintliq"]
    fields = ["snowrate", "shcomp", "lhcomp", "cldtop", "vapcomp", "divcomp", "wcomp"]
    # fields = ["w"]
    for field in fields:
        mkv.makevidcomb(folderprepath, "combplots", field)
    et = perf_counter()
    print(f"Animation took {et-st:.2f} seconds")
    
elif runtype.lower() == "both":
    st = perf_counter()
    t0 = input("Enter the start time for plotting in yyyy-mm-dd-HHMMSS format: ")
    tf = input("Enter the end time for plotting in yyyy-mm-dd-HHMMSS format: ")
    tlist = date_range(t0, tf, freq = timedelta(minutes = 20)).to_pydatetime()
    # for t in tlist:
    #     make_combplots(t)
    ppool = ProcessPoolExecutor(max_workers = 4) 
    ppool.map(make_combplots, tlist)
    ppool.shutdown()
    et = perf_counter()
    print(f"Plotting took {et-st:.2f} seconds")
    st = perf_counter()
    folderprepath = "/moonbow/ascheb/les/"
    reload(mkv)
    # fields = ["bowen"]
    # fields = ["srftemp", "rainrate", "srfpres", "vapplan", "cldtop", "uwind", "lhflux", "shflux", "bowen", "vertintcond", "vertintice", "vertintliq"]
    fields = ["snowrate", "shcomp", "lhcomp", "cldtop", "vapcomp", "divcomp", "wcomp"]
    for field in fields:
        mkv.makevidcomb(folderprepath, "combplots", field)
    et = perf_counter()
    print(f"Animation took {et-st:.2f} seconds")