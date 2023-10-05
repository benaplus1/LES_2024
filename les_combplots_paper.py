import xarray as xr
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib import colors
from matplotlib import colormaps as mcm
from matplotlib.lines import Line2D
# import matplotlib as mpl
# from metpy.units import units
# from metpy.plots import SkewT
# import metpy
from pandas import date_range
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
    
usestyle("paperplots.mplstyle")
from matplotlib import rcParams
rcParams["figure.titlesize"] = 10
rcParams["axes.titlesize"] = 8
rcParams["axes.labelsize"] = 8
rcParams["xtick.labelsize"] = 6.5
rcParams["ytick.labelsize"] = 6.5
rcParams["axes.linewidth"] = 0.2

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False):
        self.vcenter = vcenter
        super().__init__(vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        # Note also that we must extrapolate beyond vmin/vmax
        x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1.]
        return np.ma.masked_array(np.interp(value, x, y,
                                            left=-np.inf, right=np.inf))

    def inverse(self, value):
        y, x = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
        return np.interp(value, x, y, left=-np.inf, right=np.inf)

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def make_divplots(t):
    print("Plotting field Div")
    folderprepath = "/moonbow/ascheb/les/"
    figprepath = f"{folderprepath}python/PaperFigs/"
    ftime = t.strftime("%Y-%m-%d-%H%M%S")
    print(ftime)
    afile_control = xr.open_dataset(f"{folderprepath}2010/hires_control/processed_data/mergedvars_{ftime}.nc")
    afile_nolake = xr.open_dataset(f"{folderprepath}2010/hires_nolake/processed_data/mergedvars_{ftime}.nc")
    # afile_pasturebroadforest = xr.open_dataset(f"{folderpath}pasturebroadforest-wind/processed_data/mergedvars_{ftime}.nc")
    if not path.exists(f"{folderprepath}combplots"):
        mkdir(f"{folderprepath}combplots")
    modtopocmap = truncate_colormap(mcm.get_cmap("gist_earth"), minval = 0.5, maxval = 1, n = 128) #A truncated version of gist_earth to only show terrain height above sea level (gets rid of all the blue on the colorbar)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (5, 6), dpi = 200, layout = "constrained")
    axlabels = ["(a)", "(b)", "(c)", "(d)"]
    fakecontour = Line2D([], [], color = "black", linewidth = 0.375, label = "Water Bodies")
    for i, ax in enumerate(fig.get_axes()):
        ax.set_yticks([41, 42, 43, 44, 45, 46])
        ax.set_xticks([-84, -82, -80])
        ax.legend(loc = "lower left", fontsize = 5, handles = [fakecontour], labelcolor = "black")
        tercmap = ax.pcolormesh(afile_control["lon1d"][50:650], afile_control["lat1d"][450:1050], afile_control["Topo"][450:1050, 50:650], cmap = modtopocmap, vmin = 0, vmax = 1000, zorder = 0)
        ax.annotate(axlabels[i], (0.03, 0.965), xycoords = "axes fraction", color = "black", fontsize = 6, horizontalalignment = "center", bbox = {"edgecolor": "black", "linewidth": 0.2, "facecolor": "white", "pad": 1}, zorder = 3)
        # ax.pcolormesh(afile_control["lon1d"], afile_control["lat1d"], afile_control["Patch"][0,:,:].where(afile_control["Patch"][0,:,:]==1), color = "Navy", zorder = 0)
    ax3.set_xlabel("Longitude (Degrees E)")
    ax3.set_ylabel("Latitude (Degrees N)")
    divalt1 = 350 #m AMSL, altitude at which to evaluate horizontal divergence on bottom row of plots.
    divalt2 = 500 #m AMSL, altitude at which to evaluate horizontal divergence on top row of plot
    ax1.set_title(f"CONTROL {divalt1}m", fontfamily = "Liberation Serif", y = 1.02)
    ax2.set_title(f"NLH {divalt1}m", fontfamily = "Liberation Serif", y = 1.02)
    ax3.set_title(f"CONTROL {divalt2}m", fontfamily = "Liberation Serif", y = 1.02)
    ax4.set_title(f"NLH {divalt2}m", fontfamily = "Liberation Serif", y = 1.02)
    fig.suptitle(f"{divalt1}m and {divalt2}m AMSL Horizonal Divergence at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
    horizdiv_control_alt1 = np.gradient(afile_control["u"].sel(z = divalt1, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)).values, 1000, axis = 1)+np.gradient(afile_control["v"].sel(z = divalt1, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)).values, 1000, axis = 0)
    horizdiv_control_alt2 = np.gradient(afile_control["u"].sel(z = divalt2, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)).values, 1000, axis = 1)+np.gradient(afile_control["v"].sel(z = divalt2, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)).values, 1000, axis = 0)
    horizdiv_nolake_alt1 = np.gradient(afile_nolake["u"].sel(z = divalt1, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)).values, 1000, axis = 1)+np.gradient(afile_nolake["v"].sel(z = divalt1, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)).values, 1000, axis = 0)
    horizdiv_nolake_alt2 = np.gradient(afile_nolake["u"].sel(z = divalt2, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)).values, 1000, axis = 1)+np.gradient(afile_nolake["v"].sel(z = divalt2, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)).values, 1000, axis = 0)
    divmp = ax1.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), horizdiv_control_alt1, cmap = "BrBG", vmin = -5*10**(-3), vmax = 5*10**(-3), zorder = 1)
    ax1.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
    ax2.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), horizdiv_nolake_alt1, cmap = "BrBG", vmin = -5*10**(-3), vmax = 5*10**(-3), zorder = 1)
    ax2.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
    ax3.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), horizdiv_control_alt2, cmap = "BrBG", vmin = -5*10**(-3), vmax = 5*10**(-3), zorder = 1)
    ax3.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
    ax4.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), horizdiv_nolake_alt2, cmap = "BrBG", vmin = -5*10**(-3), vmax = 5*10**(-3), zorder = 1)
    ax4.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
    # ax1.legend(loc = "lower left", fontsize = 5, handles = [fakecontour], labelcolor = "black")
    # ax2.legend(loc = "lower left", fontsize = 5, handles = [fakecontour], labelcolor = "black")
    cbar1 = fig.colorbar(divmp, ax = [ax1, ax3], orientation = "horizontal", fraction = 0.05, extend = "both", pad = 0.02); cbar1.set_label("Horizontal Divergence ($\mathrm{{s^{{-1}}}}$)")
    cbar2 = fig.colorbar(tercmap, ax = [ax2, ax4], orientation = "horizontal", fraction = 0.05, extend = "max", pad = 0.02); cbar2.set_label("Terrain height (m)")
    fig.savefig(f"{figprepath}combplots/divcomp_4panel_{t.strftime('%d-%H%M')}z.png")
    plt.close(); del fig; del ax1; del ax2; del ax3; del ax4
    del cbar1; del cbar2

def make_combplots(fields, t):
    folderprepath = "/moonbow/ascheb/les/"
    figprepath = f"{folderprepath}python/PaperFigs/"
    ftime = t.strftime("%Y-%m-%d-%H%M%S")
    print(ftime)
    afile_control = xr.open_dataset(f"{folderprepath}2010/hires_control/processed_data/mergedvars_{ftime}.nc")
    afile_nolake = xr.open_dataset(f"{folderprepath}2010/hires_nolake/processed_data/mergedvars_{ftime}.nc")
    # afile_pasturebroadforest = xr.open_dataset(f"{folderpath}pasturebroadforest-wind/processed_data/mergedvars_{ftime}.nc")
    if not path.exists(f"{folderprepath}combplots"):
        mkdir(f"{folderprepath}combplots")
    
    modtopocmap = truncate_colormap(mcm.get_cmap("gist_earth"), minval = 0.5, maxval = 1, n = 128) #A truncated version of gist_earth to only show terrain height above sea level (gets rid of all the blue on the colorbar)
    for field in fields:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (5, 3.5), dpi = 200, layout = "constrained")
        axlabels = ["(a)", "(b)"]
        for i, ax in enumerate((ax1, ax2)):
            ax.set_xlabel("Longitude (Degrees E)")
            ax.set_ylabel("Latitude (Degrees N)")
            ax.set_yticks([41, 42, 43, 44, 45, 46])
            ax.set_xticks([-84, -82, -80])
            termp = ax.pcolormesh(afile_control["lon1d"][50:650], afile_control["lat1d"][450:1050], afile_control["Topo"][450:1050, 50:650], cmap = modtopocmap, vmin = 0, vmax = 1000, zorder = 0)
            ax.annotate(axlabels[i], (0.03, 0.965), xycoords = "axes fraction", color = "black", fontsize = 6, horizontalalignment = "center", bbox = {"edgecolor": "black", "linewidth": 0.2, "facecolor": "white", "pad": 1}, zorder = 3)
            # ax.pcolormesh(afile_control["lon1d"], afile_control["lat1d"], afile_control["Patch"][0,:,:].where(afile_control["Patch"][0,:,:]==1), color = "Navy", zorder = 0)
        ax1.set_title("CONTROL", fontfamily = "Liberation Serif", y = 1.02)
        ax2.set_title("NLH", fontfamily = "Liberation Serif", y = 1.02)
        
        print(f"Plotting Field {field}")
        if field == "SnowRate":
            for i, ax in enumerate(fig.get_axes()):
                ax.clear()
                ax.set_yticks([38, 40, 42, 44, 46])
                ax.set_xticks([-84, -82, -80, -78, -76])
                ax.set_xlabel("Longitude (Degrees E)")
                ax.set_ylabel("Latitude (Degrees N)")
                ax.annotate(axlabels[i], (0.03, 0.965), xycoords = "axes fraction", color = "black", fontsize = 6, horizontalalignment = "center", bbox = {"edgecolor": "black", "linewidth": 0.2, "facecolor": "white", "pad": 1}, zorder = 3)
                termp = ax.pcolormesh(afile_control["lon1d"], afile_control["lat1d"], afile_control["Topo"], cmap = modtopocmap, vmin = 0, vmax = 1000, zorder = 0)
            fig.suptitle(f"Snowfall Rate at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
            ax1.set_title("CONTROL", fontfamily = "Liberation Serif", y = 1.02)
            ax2.set_title("NLH", fontfamily = "Liberation Serif", y = 1.02)
            from palettable.lightbartlein.sequential import Blues10_10
            nohrsccmap = Blues10_10.get_mpl_colormap()
            snowbounds = [0.01, 0.1, 0.5, 1, 2, 3, 4, 5]
            snownorm = BoundaryNorm(snowbounds, ncolors = 256, extend = "max")
            fakecontour = Line2D([], [], linestyle = "-", linewidth = 0.375, color = "black", label = "Water Bodies")
            controlsnowrate = afile_control["SnowPrecipRate"]+afile_control["AggPrecipRate"]+afile_control["PrisPrecipRate"]
            nolakesnowrate = afile_nolake["SnowPrecipRate"]+afile_nolake["AggPrecipRate"]+afile_nolake["PrisPrecipRate"]
            snowmp = ax1.pcolormesh(afile_control["lon1d"], afile_control["lat1d"], controlsnowrate.where(controlsnowrate>1e-3), shading = "nearest", cmap = nohrsccmap, norm = snownorm, zorder = 1)
            ax2.pcolormesh(afile_control["lon1d"], afile_control["lat1d"], nolakesnowrate.where(nolakesnowrate>1e-3), shading = "nearest", cmap = nohrsccmap, norm = snownorm, zorder = 1)
            ax1.contour(afile_control["lon1d"], afile_control["lat1d"], afile_control["Patch"][1,:,:], levels = [0.999], colors = "black", linestyles = "-", linewidths = 0.375, zorder = 2)
            ax1.legend(loc = "lower left", fontsize = 5, handles = [fakecontour])
            ax2.contour(afile_control["lon1d"], afile_control["lat1d"], afile_nolake["Patch"][1,:,:], levels = [0.999], colors = "black", linestyles = "-", linewidths = 0.375, zorder = 2)
            ax2.legend(loc = "lower left", fontsize = 5, handles = [fakecontour])
            snowcbar = fig.colorbar(snowmp, ax = ax1, orientation = "horizontal", fraction = 0.05, extend = "max"); snowcbar.set_label("Snowfall Rate ($\mathrm{mm \ hr^{-1}}$)")
            tercbar = fig.colorbar(termp, ax = ax2, orientation = "horizontal", fraction = 0.05, extend = "max"); tercbar.set_label("Terrain height (m)")
            fig.savefig(f"{figprepath}combplots/snowrate_{str(t.day).zfill(2)}{str(t.hour).zfill(2)}{str(t.minute).zfill(2)}z.png")

        elif field == "w":
            walt = 700 #m AMSL, altitude at which to evalute w
            fig.suptitle(f"700m AMSL Vertical Velocity at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
            fakecontour = Line2D([], [], color = "black", linewidth = 0.375, label = "Water Bodies")
            wmp = ax1.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["w"].sel(z = walt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "RdBu_r", vmin = -3, vmax = 3, zorder = 1)
            ax1.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
            ax2.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["w"].sel(z = walt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "RdBu_r", vmin = -3, vmax = 3, zorder = 1)
            ax2.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
            ax1.legend(loc = "lower left", fontsize = 5, handles = [fakecontour], labelcolor = "black")
            ax2.legend(loc = "lower left", fontsize = 5, handles = [fakecontour], labelcolor = "black")
            cbar = fig.colorbar(wmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"{walt}m AMSL Vertical Velocity ($\mathrm{{m \ s^{{-1}}}}$)")
            fig.savefig(f"{figprepath}combplots/wcomp_{t.strftime('%d-%H%M')}z.png")
            
        elif field == "VapMix":
            vapalt = 700 #m AMSL, altitude at which to evaluate horizontal divergence.
            fakecontour = Line2D([], [], color = "black", linewidth = 0.375, label = "Water Bodies")
            fig.suptitle(f"{vapalt}m AMSL Vapor Mixing Ratio at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
            vapmp = ax1.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["VaporMix"].sel(z = vapalt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = 0, vmax = 2, zorder = 1)
            ax1.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
            ax2.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["VaporMix"].sel(z = vapalt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = 0, vmax = 2, zorder = 1)
            ax2.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
            ax1.legend(loc = "lower left", fontsize = 5, handles = [fakecontour], labelcolor = "black")
            ax2.legend(loc = "lower left", fontsize = 5, handles = [fakecontour], labelcolor = "black")
            cbar = fig.colorbar(vapmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"{vapalt}m AMSL Vapor Mixing Ratio ($\mathrm{{g \ kg^{{-1}}}}$)")
            fig.savefig(f"{figprepath}combplots/vapcomp_{t.strftime('%d-%H%M')}z.png")
            
        elif field == "ShFlux":
            fakecontour = Line2D([], [], color = "black", linewidth = 0.375, label = "Water Bodies")
            fig.suptitle(f"Surface Sensible Heat Flux at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
            shmp = ax1.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["SensibleHeatFlux"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "bwr", vmin = -500, vmax = 500, zorder = 1)
            ax1.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
            ax2.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["SensibleHeatFlux"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "bwr", vmin = -500, vmax = 500, zorder = 1)
            ax2.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
            ax1.legend(loc = "lower left", fontsize = 5, handles = [fakecontour], labelcolor = "black")
            ax2.legend(loc = "lower left", fontsize = 5, handles = [fakecontour], labelcolor = "black")
            cbar = fig.colorbar(shmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"Sensible Heat Flux ($\mathrm{{W \ m^{{-2}}}}$)")
            fig.savefig(f"{figprepath}combplots/shcomp_{t.strftime('%d-%H%M')}z.png")
            
        elif field == "LhFlux":
            fakecontour = Line2D([], [], color = "black", linewidth = 0.375, label = "Water Bodies")
            fig.suptitle(f"Surface Latent Heat Flux at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
            lhmp = ax1.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["LatentHeatFlux"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = -500, vmax = 500, zorder = 1)
            ax1.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
            ax2.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["LatentHeatFlux"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = -500, vmax = 500, zorder = 1)
            ax2.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
            ax1.legend(loc = "lower left", fontsize = 5, handles = [fakecontour], labelcolor = "black")
            ax2.legend(loc = "lower left", fontsize = 5, handles = [fakecontour], labelcolor = "black")
            cbar = fig.colorbar(lhmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"Latent Heat Flux ($\mathrm{{W \ m^{{-2}}}}$)")
            fig.savefig(f"{figprepath}combplots/lhcomp_{t.strftime('%d-%H%M')}z.png")
            
        elif field == "CldTop":
            from palettable.cmocean.sequential import Tempo_10
            cldcmap = Tempo_10.get_mpl_colormap().reversed()
            fakecontour = Line2D([], [], color = "black", linewidth = 0.375, label = "Water Bodies")
            fig.suptitle(f"Cloud Top Height at {(t-timedelta(hours=5)).strftime('%d')} Jan - {(t-timedelta(hours=5)).strftime('%H%M')} LT")
            cldtopmp = ax1.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["CloudTop"].where(afile_control["CloudTop"]>0).isel(x = slice(50, 650), y = slice(450, 1050)), cmap = cldcmap, vmin = 0, vmax = 3000, zorder = 1)
            ax1.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_control["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
            ax2.pcolormesh(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["CloudTop"].where(afile_nolake["CloudTop"]>0).isel(x = slice(50, 650), y = slice(450, 1050)), cmap = cldcmap, vmin = 0, vmax = 3000, zorder = 1)
            ax2.contour(afile_control["lon1d"].isel(x = slice(50, 650)), afile_control["lat1d"].isel(y = slice(450, 1050)), afile_nolake["Patch"].isel(patch = 1, x = slice(50, 650), y = slice(450, 1050)), levels = [0.99], colors = "black", linewidths = 0.375, zorder = 2)
            ax1.legend(loc = "lower left", fontsize = 5, handles = [fakecontour], labelcolor = "black")
            ax2.legend(loc = "lower left", fontsize = 5, handles = [fakecontour], labelcolor = "black")
            cldcbar = fig.colorbar(cldtopmp, ax = ax1, orientation = "horizontal", fraction = 0.05, extend = "max"); cldcbar.set_label(f"Cloud Top Height (m AMSL)")
            tercbar = fig.colorbar(termp, ax = ax2, orientation = "horizontal", fraction = 0.05, extend = "max"); tercbar.set_label("Terrain Height (m)")
            fig.savefig(f"{figprepath}combplots/cldtop_{t.strftime('%d-%H%M')}z.png")
            
        plt.close(); del fig; del ax1; del ax2;
        try:
            del cbar
        except:
            pass
        try:
            del cbar1
        except:
            pass
        try:
            del cbar2
        except:
            pass

    afile_nolake.close(); del afile_nolake
    afile_control.close(); del afile_control
    return ftime
runtype = input("Do you want to *plot*, *animate* existing plots, or do *both*? ")
if runtype.lower() == "plot":
    st = perf_counter()
    t0 = input("Enter the start time for plotting in yyyy-mm-dd-HHMMSS format: ")
    tf = input("Enter the end time for plotting in yyyy-mm-dd-HHMMSS format: ")
    fields = input("Enter the list of fields you want to plot, as comma-separated values: ")
    fields = fields.split(",")
    fields = [i.strip() for i in fields]
    if "Div" in fields:
        fields.remove("Div")
    plotdiv = True
    partial_combplots = partial(make_combplots, fields)
    tlist = date_range(t0, tf, freq = timedelta(minutes = 20)).to_pydatetime()
    seq = input("Is this a test run? Yes or No? ")
    if seq.lower() == "yes":
        for t in tlist:
            make_combplots(fields, t)
            if plotdiv:
                make_divplots(t)
    elif seq.lower() == "no":
        ppool = ProcessPoolExecutor(max_workers = 4) 
        ppool.map(partial_combplots, tlist)
        if plotdiv:
            ppool.map(make_divplots, tlist)
        ppool.shutdown()
    else:
        raise Exception("Must be *yes* or *no*!")
    et = perf_counter()
    print(f"Plotting took {et-st:.2f} seconds")
    
    
elif runtype.lower() == "animate":
    st = perf_counter()
    reload(mkv)
    # fields = ["bowen"]
    # fields = ["srftemp", "rainrate", "srfpres", "vapplan", "cldtop", "uwind", "lhflux", "shflux", "bowen", "vertintcond", "vertintice", "vertintliq"]
    fields = ["snowrate", "shcomp", "lhcomp", "cldtop", "vapcomp", "divcomp", "wcomp"]
    # fields = ["w"]
    for field in fields:
        mkv.makevidcomb(figprepath, "combplots", field)
    et = perf_counter()
    print(f"Animation took {et-st:.2f} seconds")
    
elif runtype.lower() == "both":
    st = perf_counter()
    t0 = input("Enter the start time for plotting in yyyy-mm-dd-HHMMSS format: ")
    tf = input("Enter the end time for plotting in yyyy-mm-dd-HHMMSS format: ")
    fields = fields.split(",")
    fields = [i.strip() for i in fields]
    if "Div" in fields:
        fields.remove("Div")
    plotdiv = True
    partial_combplots = partial(make_combplots, fields)
    tlist = date_range(t0, tf, freq = timedelta(minutes = 20)).to_pydatetime()
    seq = input("Is this a test run? Yes or No? ")
    if seq.lower() == "yes":
        for t in tlist:
            make_combplots(fields, t)
            if plotdiv:
                make_divplots(t)
    elif seq.lower() == "no":
        ppool = ProcessPoolExecutor(max_workers = 4) 
        ppool.map(partial_combplots, tlist)
        if plotdiv:
            ppool.map(make_divplots, tlist)
        ppool.shutdown()
    et = perf_counter()
    print(f"Plotting took {et-st:.2f} seconds")
    st = perf_counter()
    reload(mkv)
    # fields = ["bowen"]
    # fields = ["srftemp", "rainrate", "srfpres", "vapplan", "cldtop", "uwind", "lhflux", "shflux", "bowen", "vertintcond", "vertintice", "vertintliq"]
    fields = ["snowrate", "shcomp", "lhcomp", "cldtop", "vapcomp", "divcomp", "wcomp"]
    for field in fields:
        mkv.makevidcomb(figprepath, "combplots", field)
    et = perf_counter()
    print(f"Animation took {et-st:.2f} seconds")
else:
    raise Exception("Runtype must be either plot, animate, or both!")
