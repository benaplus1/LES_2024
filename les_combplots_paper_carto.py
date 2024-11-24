#test whether git pushing is working
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib import colors
from matplotlib import colormaps as mcm
from matplotlib.lines import Line2D
from pandas import date_range
from time import perf_counter
from glob import glob
from matplotlib.style import use as usestyle
from importlib import reload
import makevids as mkv
from concurrent.futures import ProcessPoolExecutor
from os import path, mkdir
from functools import partial
from cartopy import feature as cfeature
import cartopy.crs as ccrs
from matplotlib import ticker as mticker

from matplotlib import font_manager as fm
fontdir = "/home/ascheb/libfonts/*.ttf"
for fpath in glob(fontdir):
    # print(fpath)
    fm.fontManager.addfont(fpath)
    
usestyle("paperplots.mplstyle")
from matplotlib import rcParams
rcParams["figure.titlesize"] = 15
rcParams["axes.titlesize"] = 12
rcParams["axes.labelsize"] = 12
rcParams["xtick.labelsize"] = 9
rcParams["ytick.labelsize"] = 9
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

def add_cartofeatures_fulldomain(afile_control, geoax):
    ramscrs = ccrs.Stereographic(central_longitude = -80, central_latitude = 42)
    geoax.set_extent((afile_control["x"].min(), afile_control["x"].max(), afile_control["y"].min(), afile_control["y"].max()), crs = ramscrs)
    geoax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor = 'darkgrey', facecolor = "none", linewidth = 0.4, zorder = 3))
    geoax.add_feature(cfeature.NaturalEarthFeature(category="cultural", name = "admin_1_states_provinces_lakes", linewidth = 0.4, facecolor = "none", scale = "50m", edgecolor = "black", zorder = 3))
    return geoax

def add_cartofeatures_zoomdomain(afile_control, geoax):
    ramscrs = ccrs.Stereographic(central_longitude = -80, central_latitude = 42)
    geoax.set_extent((afile_control["x"][50], afile_control["x"][650], afile_control["y"][450], afile_control["y"][1050]), crs = ramscrs)
    geoax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor = 'darkgrey', facecolor = "none", linewidth = 0.4, zorder = 3))
    geoax.add_feature(cfeature.NaturalEarthFeature(category="cultural", name = "admin_1_states_provinces_lakes", linewidth = 0.4, facecolor = "none", scale = "50m", edgecolor = "black", zorder = 3))
    
    gl = geoax.gridlines(crs = ccrs.PlateCarree(), draw_labels = {"bottom": "x", "left": "y"}, x_inline = False, y_inline = False, dms = True)
    gl.xlocator = mticker.FixedLocator([-84, -82, -80, -78])
    gl.ylocator = mticker.FixedLocator([42, 43, 44, 45, 46])
    gl.bottom_labels = True
    gl.left_labels   = True
    gl.xlines = False
    gl.ylines = False
    gl.top_labels    = False
    gl.right_labels  = False
    gl.xlabel_style = {'size': 9, 'color': 'black', 'rotation': 0, "horizontalalignment": "right"}
    gl.ylabel_style = {'size': 9, 'color': 'black', 'rotation': 0, "horizontalalignment": "right"}
    return geoax



def make_divplots(figprepath, controlprepath, nlhprepath, t):
    ramscrs = ccrs.Stereographic(central_longitude = -80, central_latitude = 42)
    print("Plotting field Div")
    ftime = t.strftime("%Y-%m-%d-%H%M%S")
    print(ftime)
    afile_control = xr.open_dataset(f"{controlprepath}/mvars-cart-{ftime}-g1.nc")
    afile_nolake = xr.open_dataset(f"{nlhprepath}/mvars-cart-{ftime}-g1.nc")
    if not path.exists(f"{figprepath}combplots"):
        mkdir(f"{figprepath}combplots")
    modtopocmap = truncate_colormap(mcm.get_cmap("gist_earth"), minval = 0.5, maxval = 1, n = 128) #A truncated version of gist_earth to only show terrain height above sea level (gets rid of all the blue on the colorbar)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (5, 6), dpi = 200, layout = "compressed", subplot_kw = {"projection": ramscrs})
    axlabels = ["(a)", "(b)", "(c)", "(d)"]
    fakecontour = Line2D([], [], color = "black", linewidth = 0.375, label = "Water Bodies")
    for i, ax in enumerate(fig.get_axes()):
        ax = add_cartofeatures_zoomdomain(afile_control, ax)
        left, width = 0, 0.1
        bottom, height = 0.90, 0.1
        right = left + width
        top = bottom + height
        p = plt.Rectangle((left, bottom), width, height, fill=True, zorder = 3, edgecolor = "black", linewidth = 0.2, facecolor = "white")
        p.set_transform(ax.transAxes)
        p.set_clip_on(False)
        ax.add_patch(p)
        ax.text(left+0.5*width, bottom+0.5*height, axlabels[i], fontsize = 10, transform = ax.transAxes, horizontalalignment = "center", verticalalignment = "center", zorder = 5, color = "black")
        tercmap = ax.pcolormesh(afile_control["x"][50:650], afile_control["y"][450:1050], afile_control["Topo"][450:1050, 50:650], cmap = modtopocmap, vmin = 0, vmax = 1000, zorder = 0, transform = ramscrs)
    divalt1 = 350 #m AMSL, altitude at which to evaluate horizontal divergence on bottom row of plots.
    divalt2 = 1500 #m AMSL, altitude at which to evaluate horizontal divergence on top row of plot
    fig.suptitle(f"Horizontal Divergence at {t.strftime('%d')} Jan - {t.strftime('%H%M')} UTC")
    ax1.set_title(f"CONTROL {divalt1}m")
    ax2.set_title(f"NLH {divalt1}m")
    ax3.set_title(f"CONTROL {divalt2}m")
    ax4.set_title(f"NLH {divalt2}m")
    divmp = ax1.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), afile_control["HorizDiv"].sel(z = divalt1, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = -5*10**(-3), vmax = 5*10**(-3), zorder = 1, transform = ramscrs)
    ax2.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), afile_nolake["HorizDiv"].sel(z = divalt1, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = -5*10**(-3), vmax = 5*10**(-3), zorder = 1, transform = ramscrs)
    ax3.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), afile_control["HorizDiv"].sel(z = divalt2, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = -5*10**(-3), vmax = 5*10**(-3), zorder = 1, transform = ramscrs)
    ax4.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), afile_nolake["HorizDiv"].sel(z = divalt2, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = -5*10**(-3), vmax = 5*10**(-3), zorder = 1, transform = ramscrs)
    cbar1 = fig.colorbar(divmp, ax = [ax1, ax2, ax3, ax4], orientation = "horizontal", fraction = 0.05, aspect = 40, extend = "both", pad = 0.02); cbar1.set_label("Horizontal Divergence ($\mathrm{{s^{{-1}}}}$)", labelpad = 7)
    cbar2 = fig.colorbar(tercmap, ax = [ax1, ax2, ax3, ax4], orientation = "vertical", location = "left", fraction = 0.05, aspect = 40, extend = "max", pad = 0.02); cbar2.set_label("Terrain height (m)", labelpad = 7)
    for cbari in [cbar1, cbar2]:
        cbari.ax.tick_params(labelsize = 10)
    cbar1.formatter.set_powerlimits((-1, 1))
    fig.savefig(f"{figprepath}/combplots/div_carto_4panel_{t.strftime('%d-%H%M')}z.png")
    plt.close(); del fig; del ax1; del ax2; del ax3; del ax4
    del cbar1; del cbar2

def make_combplots(figprepath, controlprepath, nlhprepath, fields, t):
    ramscrs = ccrs.Stereographic(central_longitude = -80, central_latitude = 42)
    ftime = t.strftime("%Y-%m-%d-%H%M%S")
    print(ftime)
    afile_control = xr.open_dataset(f"{controlprepath}/mvars-cart-{ftime}-g1.nc")
    afile_nolake = xr.open_dataset(f"{nlhprepath}/mvars-cart-{ftime}-g1.nc")
    if not path.exists(f"{figprepath}/combplots"):
        mkdir(f"{figprepath}/combplots")
    
    modtopocmap = truncate_colormap(mcm.get_cmap("gist_earth"), minval = 0.5, maxval = 1, n = 128) #A truncated version of gist_earth to only show terrain height above sea level (gets rid of all the blue on the colorbar)
    for field in fields:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (5, 3.5), dpi = 200, layout = "compressed", subplot_kw = {"projection": ramscrs})
        axlabels = ["(a)", "(b)"]
        for i, ax in enumerate(fig.get_axes()):
            if field != "SnowRate":
                ax = add_cartofeatures_zoomdomain(afile_control, ax)
            left, width = 0, 0.1
            bottom, height = 0.90, 0.1
            right = left + width
            top = bottom + height
            p = plt.Rectangle((left, bottom), width, height, fill=True, zorder = 4, edgecolor = "black", linewidth = 0.2, facecolor = "white")
            p.set_transform(ax.transAxes)
            p.set_clip_on(False)
            ax.add_patch(p)
            ax.text(left+0.5*width, bottom+0.5*height, axlabels[i], fontsize = 10, transform = ax.transAxes, horizontalalignment = "center", verticalalignment = "center", zorder = 5, color = "black")
            tercmap = ax.pcolormesh(afile_control["x"][50:650], afile_control["y"][450:1050], afile_control["Topo"][450:1050, 50:650], cmap = modtopocmap, vmin = 0, vmax = 1000, zorder = 0, transform = ramscrs)
            
        ax1.set_title("CONTROL")
        ax2.set_title("NLH")
        
        print(f"Plotting Field {field}")
        if field == "snowrate":
            for i, ax in enumerate(fig.get_axes()):
                ax = add_cartofeatures_fulldomain(afile_control, ax)
                tercmap = ax.pcolormesh(afile_control["x"], afile_control["y"], afile_control["Topo"], cmap = modtopocmap, vmin = 0, vmax = 1000, zorder = 0, transform = ramscrs)
            ax1.set_title("CONTROL")
            ax2.set_title("NLH")
            gl = ax2.gridlines(crs = ccrs.PlateCarree(), draw_labels = {"bottom": "x", "left": "y"}, x_inline = False, y_inline = False, dms = True)
            gl.xlocator = mticker.FixedLocator([-84, -82, -80, -78, -76])
            gl.ylocator = mticker.FixedLocator([38, 40, 42, 44, 46])
            gl.bottom_labels = True
            gl.left_labels   = True
            gl.xlines = False
            gl.ylines = False
            gl.top_labels    = False
            gl.right_labels  = False
            gl.xlabel_style = {'size': 9, 'color': 'black', 'rotation': 0, "horizontalalignment": "right"}
            gl.ylabel_style = {'size': 9, 'color': 'black', 'rotation': 0, "horizontalalignment": "right"}
            from palettable.lightbartlein.sequential import Blues10_10
            nohrsccmap = Blues10_10.get_mpl_colormap()
            snowbounds = [0.01, 0.1, 0.5, 1, 2, 3, 4, 5]
            snownorm = BoundaryNorm(snowbounds, ncolors = 256, extend = "max")
            fakecontour = Line2D([], [], linestyle = "-", linewidth = 0.375, color = "black", label = "Water Bodies")
            fig.suptitle(f"Snowfall Rate at {t.strftime('%d')} Jan - {t.strftime('%H%M')} UTC")
            controlsnowrate = afile_control["SnowPrecipRate"]+afile_control["AggPrecipRate"]+afile_control["PrisPrecipRate"]
            nolakesnowrate = afile_nolake["SnowPrecipRate"]+afile_nolake["AggPrecipRate"]+afile_nolake["PrisPrecipRate"]
            ax1.pcolormesh(afile_control["x"], afile_control["y"], afile_control["PatchArea"][0,:,:].where(afile_control["PatchArea"][0,:,:]==1), cmap = "Blues", vmin = 0, vmax = 1, zorder = 1, transform = ramscrs)
            ax2.pcolormesh(afile_control["x"], afile_control["y"], afile_nolake["PatchArea"][0,:,:].where(afile_nolake["PatchArea"][0,:,:]==1), cmap = "Blues", vmin = 0, vmax = 1, zorder = 1, transform = ramscrs)
            snowmp = ax1.pcolormesh(afile_control["x"], afile_control["y"], controlsnowrate.where(controlsnowrate>1e-3), shading = "nearest", cmap = nohrsccmap, norm = snownorm, zorder = 2, transform = ramscrs)
            ax2.pcolormesh(afile_control["x"], afile_control["y"], nolakesnowrate.where(nolakesnowrate>1e-3), shading = "nearest", cmap = nohrsccmap, norm = snownorm, zorder = 2, transform = ramscrs)
            snowcbar = fig.colorbar(snowmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, aspect = 40, extend = "max"); snowcbar.set_label("Snowfall Rate ($\mathrm{mm \ hr^{-1}}$)")
            snowcbar.set_ticks(ticks = snowbounds, labels = ["0.01", "0.1", "0.5", "1", "2", "3", "4", "5"])
            tercbar = fig.colorbar(tercmap, ax = [ax1, ax2], orientation = "vertical", location = "left", fraction = 0.05, extend = "max", pad = 0.02); tercbar.set_label("Terrain height (m)");
            fig.savefig(f"{figprepath}/combplots/snowrate_carto_{t.strftime('%d-%H%M')}z.png")

        elif field == "w":
            walt = 700 #m AMSL, altitude at which to evalute w
            fig.suptitle(f"{walt} m AMSL Vertical Velocity at {t.strftime('%d')} Jan - {t.strftime('%H%M')} UTC")
            wmp = ax1.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), afile_control["w"].sel(z = walt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "RdBu_r", vmin = -3, vmax = 3, zorder = 1, transform = ramscrs)
            ax2.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), afile_nolake["w"].sel(z = walt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "RdBu_r", vmin = -3, vmax = 3, zorder = 1, transform = ramscrs)
            cbar = fig.colorbar(wmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"{walt}m AMSL Vertical Velocity ($\mathrm{{m \ s^{{-1}}}}$)")
            fig.savefig(f"{figprepath}/combplots/w_carto_{t.strftime('%d-%H%M')}z.png")
            
        elif field == "vapmix":
            vapalt = 700 #m AMSL, altitude at which to evaluate vapor mixing ratio.
            fig.suptitle(f"{vapalt} m AMSL Water Vapor at {t.strftime('%d')} Jan - {t.strftime('%H%M')} UTC")
            vapmp = ax1.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), afile_control["VaporMix"].sel(z = vapalt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = 0, vmax = 2, zorder = 1, transform = ramscrs)
            ax2.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), afile_nolake["VaporMix"].sel(z = vapalt, method = "nearest").isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = 0, vmax = 2, zorder = 1, transform = ramscrs)
            cbar = fig.colorbar(vapmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"{vapalt}m AMSL Vapor Mixing Ratio ($\mathrm{{g \ kg^{{-1}}}}$)")
            fig.savefig(f"{figprepath}/combplots/vapmix_carto_{t.strftime('%d-%H%M')}z.png")
            
        elif field == "shflux":
            fig.suptitle(f"Surface Sensible Heat Flux at {t.strftime('%d')} Jan - {t.strftime('%H%M')} UTC")
            shmp = ax1.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), afile_control["SensibleHeatFlux"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "bwr", vmin = -500, vmax = 500, zorder = 1, transform = ramscrs)
            ax2.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), afile_nolake["SensibleHeatFlux"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "bwr", vmin = -500, vmax = 500, zorder = 1, transform = ramscrs)
            cbar = fig.colorbar(shmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"Sensible Heat Flux ($\mathrm{{W \ m^{{-2}}}}$)")
            fig.savefig(f"{figprepath}/combplots/shflux_carto_{t.strftime('%d-%H%M')}z.png")
            
        elif field == "lhflux":
            fig.suptitle(f"Surface Latent Heat Flux at {t.strftime('%d')} Jan - {t.strftime('%H%M')} UTC")
            lhmp = ax1.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), afile_control["LatentHeatFlux"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = -500, vmax = 500, zorder = 1, transform = ramscrs)
            ax2.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), afile_nolake["LatentHeatFlux"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "BrBG", vmin = -500, vmax = 500, zorder = 1, transform = ramscrs)
            cbar = fig.colorbar(lhmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"Latent Heat Flux ($\mathrm{{W \ m^{{-2}}}}$)")
            fig.savefig(f"{figprepath}/combplots/lhflux_carto_{t.strftime('%d-%H%M')}z.png")

        elif field == "srftemp":
            fig.suptitle(f"15 m Air Temperature at {t.strftime('%d')} Jan - {t.strftime('%H%M')} UTC")
            tempmp = ax1.pcolormesh(afile_control["x"].isel(x = slice(50,650)), afile_control["y"].isel(y = slice(450, 1050)), afile_control["SrfTemp"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "Spectral_r", vmin = 245, vmax = 275, zorder = 1, transform = ramscrs)
            ax2.pcolormesh(afile_nolake["x"].isel(x = slice(50, 650)), afile_nolake["y"].isel(y = slice(450, 1050)), afile_nolake["SrfTemp"].isel(x = slice(50, 650), y = slice(450, 1050)), cmap = "Spectral_r", vmin = 245, vmax = 275, zorder = 1, transform = ramscrs)
            cbar = fig.colorbar(tempmp, ax = [ax1, ax2], orientation = "horizontal", fraction = 0.05, extend = "both"); cbar.set_label(f"15 m AGL Air Temperature (K)")
            fig.savefig(f"{figprepath}/combplots/srftemp_carto_{t.strftime('%d-%H%M')}z.png")
        elif field == "cldtop":
            from palettable.cmocean.sequential import Tempo_10
            cldcmap = Tempo_10.get_mpl_colormap().reversed()
            import warnings
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', r'All-NaN slice encountered')
            #Numpy will freak out when we try to find the cloud height on a column where there are no clouds, so this quiets that down.
            #Need to limit our cloud top filter to 2km to screen out mid-level overlying clouds
            cldtop_control = np.nanmax(np.where(afile_control["CloudMassMix"].isel(x = slice(50, 650), y = slice(450, 1050)).sel(z = slice(0, 2000)).values+afile_control["PrisMassMix"].isel(x = slice(50, 650), y = slice(450, 1050)).sel(z = slice(0, 2000)).values > 1e-5, afile_control["z"].sel(z = slice(0,2000)).values[:,None,None], np.nan), axis = 0)
            cldtop_nolake = np.nanmax(np.where(afile_nolake["CloudMassMix"].isel(x = slice(50, 650), y = slice(450, 1050)).sel(z = slice(0, 2000)).values+afile_nolake["PrisMassMix"].isel(x = slice(50, 650), y = slice(450, 1050)).sel(z = slice(0, 2000)).values > 1e-5, afile_nolake["z"].sel(z = slice(0,2000)).values[:,None,None], np.nan), axis = 0)
            fig.suptitle(f"Cloud Top Height at {t.strftime('%d')} Jan - {t.strftime('%H%M')} UTC")
            cldtopmp = ax1.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), cldtop_control, cmap = cldcmap, vmin = 0, vmax = 2000, zorder = 1, transform = ramscrs)
            ax2.pcolormesh(afile_control["x"].isel(x = slice(50, 650)), afile_control["y"].isel(y = slice(450, 1050)), cldtop_nolake, cmap = cldcmap, vmin = 0, vmax = 2000, zorder = 1, transform = ramscrs)
            cldcbar = fig.colorbar(cldtopmp, ax = [ax1], orientation = "horizontal", fraction = 0.05, extend = "max"); cldcbar.set_label(f"Cloud Top Height (m AMSL)")
            tercbar = fig.colorbar(tercmap, ax = [ax2], orientation = "horizontal", fraction = 0.05, extend = "max"); tercbar.set_label("Terrain Height (m)")
            fig.savefig(f"{figprepath}/combplots/cldtop_carto_{t.strftime('%d-%H%M')}z.png")
            
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
    figprepath = input("Enter the directory where you want the figures to go: ").rstrip("/ ")
    print("Note that whatever folder path you just picked, a new subfolder will be created inside called "
          "'combplots', where these plots will be placed. This is done to avoid cluttering up the directory if you decide"
          " to put all figures for recreation in the same folder.")
    controlprepath = input("Enter the directory containing the post-processed data for the CONTROL simulation: ").rstrip("/")
    if not path.exists(controlprepath):
        raise Exception("CONTROL post-processed files not found!")
    nlhprepath = input("Enter the directory containing the post-processed data for the NLH simulation: ").rstrip("/ ")
    st = perf_counter()
    if not path.exists(nlhprepath):
        raise Exception("NLH post-processed files not found!")
    t0 = input("Enter the start time for plotting in yyyy-mm-dd-HHMMSS format: ")
    tf = input("Enter the end time for plotting in yyyy-mm-dd-HHMMSS format: ")
    if (datetime.strptime(tf, "%Y-%m-%d-%H%M%S") - datetime.strptime(t0, "%Y-%m-%d-%H%M%S") < timedelta(seconds = 1)):
        raise Exception("End time must be later than start time!")
    print(f"{controlprepath}/mvars-cart-{t0}-g1.nc") 
    print(f"{nlhprepath}/mvars-cart-{t0}-g1.nc")
    if (not path.exists(f"{controlprepath}/mvars-cart-{t0}-g1.nc")) or (not path.exists(f"{nlhprepath}/mvars-cart-{t0}-g1.nc")):
        raise Exception("Selected start time not available in post-processed data!")
    if (not path.exists(f"{controlprepath}/mvars-cart-{tf}-g1.nc")) or (not path.exists(f"{nlhprepath}/mvars-cart-{tf}-g1.nc")):
        raise Exception("Selected end time not available in post-processed data!")
    fields = input("Enter the list of fields you want to plot, as comma-separated values: ")
    '''
    Available fields are:
    snowrate: Plots the instantaneous snowfall rate across the entire domain
    cldtop: Plots the instantaneous cloud top height over a subdomain downwind of Lake Erie
    shflux: Plots the instantaneous surface sensible heat flux over a subdomain downwind of Lake Erie
    lhflux: Plots the instantaneous latent heat flux over a subdomain downwind of Lake Erie
    vapmix: Plots the instantaneous water vapor mixing ratio at a user-chosen altitude over a subdomain downwind of Lake Erie
    w: Plots the instantaneous vertical velocity at a user-chosen altitude over a subdomain downwind of Lake Erie
    div: Plots the instantaneous horizontal divergence at two user-chosen altitudes over a subdomain downwind of Lake Erie
    '''
    fields = fields.split(",")
    fields = [i.strip().lower() for i in fields]
    print(fields)
    plotdiv = False
    if "div" in fields:
        fields.remove("div")
        plotdiv = True
    partial_combplots = partial(make_combplots, figprepath, controlprepath, nlhprepath, fields)
    partial_divplots = partial(make_divplots, figprepath, controlprepath, nlhprepath)
    tlist = date_range(t0, tf, freq = timedelta(minutes = 20)).to_pydatetime()
    seq = input("Is this a test run? Yes or No? ")
    if seq.lower() == "yes":
        for t in tlist:
            if len(fields) > 0:
                make_combplots(figprepath, controlprepath, nlhprepath, fields, t)
            if plotdiv:
                make_divplots(figprepath, controlprepath, nlhprepath, t)
    elif seq.lower() == "no":
        ppool = ProcessPoolExecutor(max_workers = 4) 
        if len(fields)> 0:
            ppool.map(partial_combplots, tlist)
        if plotdiv:
            ppool.map(partial_divplots, tlist)
        ppool.shutdown()
    else:
        raise Exception("Must be *yes* or *no*!")
    et = perf_counter()
    print(f"Plotting took {et-st:.2f} seconds")
    
    
elif runtype.lower() == "animate":
    figprepath = input("Enter the directory containing the 'combplots' folder"
                       " where the plan view frames are saved: ").rstrip("/ ")
    figprepath = f"{figprepath}/combplots"
    if not path.exists(figprepath):
        raise Exception("combplots folder not found in the specified subdirectory!")
    st = perf_counter()
    reload(mkv)
    fields = input("Enter the list of fields you want to animate, as comma-separated values: ")
    '''
    Available fields are:
    snowrate: Instantaneous snowfall rate across the entire domain
    cldtop: Instantaneous cloud top height over the northwest part of the model domain
    shflux: Instantaneous surface sensible heat flux over the northwest part of the model domain
    lhflux: Instantaneous latent heat flux over the northwest part of the model domain
    srftemp: Instantaneous 15 m AGL temperature over the northwest part of the model domain
    vapmix: Instantaneous water vapor mixing ratio at a user-chosen altitude over the northwest part of the model domain
    w: Instantaneous vertical velocity at a user-chosen altitude over the northwest part of the model domain
    div: Instantaneous horizontal divergence at two user-chosen altitudes over the northwest part of the model domain
    '''
    fields = fields.split(",")
    fields = [i.strip().lower() for i in fields]
    # fields = ["w"]
    for field in fields:
        mkv.makevidcomb(figprepath, f"{field}_carto")
    et = perf_counter()
    print(f"Animation took {et-st:.2f} seconds")
    
elif runtype.lower() == "both":
    st = perf_counter()
    figprepath = input("Enter the directory where you want the figures to go: ").rstrip("/ ")
    print("Note that whatever folder path you just picked, a new subfolder will be created inside called"
          "'combplots', where these plots will be placed. This is done to avoid cluttering up the directory if you decide"
          "to put all figures for recreation in the same folder.")
    controlprepath = input("Enter the directory containing the post-processed data for the CONTROL simulation: ").rstrip("/")
    if not path.exists(controlprepath):
        raise Exception("CONTROL post-processed files not found!")
    nlhprepath = input("Enter the directory containing the post-processed data for the NLH simulation: ").rstrip("/ ")
    if not path.exists(nlhprepath):
        raise Exception("NLH post-processed files not found!")
    t0 = input("Enter the start time for plotting in yyyy-mm-dd-HHMMSS format: ")
    tf = input("Enter the end time for plotting in yyyy-mm-dd-HHMMSS format: ")
    if (datetime.strptime(tf, "%Y-%m-%d-%H%M%S") - datetime.strptime(t0, "%Y-%m-%d-%H%M%S") < timedelta(seconds = 1)):
        raise Exception("End time must be later than start time!")
    if (not path.exists(f"{controlprepath}/mvars-cart-{t0}-g1.nc")) or (not path.exists(f"{nlhprepath}/mvars-cart-{t0}-g1.nc")):
        raise Exception("Selected start time not available in post-processed data!")
    if (not path.exists(f"{controlprepath}/mvars-cart-{tf}-g1.nc")) or (not path.exists(f"{nlhprepath}/mvars-cart-{tf}-g1.nc")):
        raise Exception("Selected end time not available in post-processed data!")
    fields = input("Enter the list of fields you want to plot, as comma-separated values: ")
    '''
    Available fields are:
    snowrate: Instantaneous snowfall rate across the entire domain
    cldtop: Instantaneous cloud top height over the northwest part of the model domain
    shflux: Instantaneous surface sensible heat flux over the northwest part of the model domain
    lhflux: Instantaneous latent heat flux over the northwest part of the model domain
    srftemp: Instantaneous 15 m AGL temperature over the northwest part of the model domain
    vapmix: Instantaneous water vapor mixing ratio at a user-chosen altitude over the northwest part of the model domain
    w: Instantaneous vertical velocity at a user-chosen altitude over the northwest part of the model domain
    div: Instantaneous horizontal divergence at two user-chosen altitudes over the northwest part of the model domain
    '''
    fields = fields.split(",")
    fields = [i.strip().lower() for i in fields]
    print(fields)
    plotdiv = False
    if "div" in fields:
        fields.remove("div")
        plotdiv = True
    partial_combplots = partial(make_combplots, figprepath, controlprepath, nlhprepath, fields)
    partial_divplots = partial(make_divplots, figprepath, controlprepath, nlhprepath)
    tlist = date_range(t0, tf, freq = timedelta(minutes = 20)).to_pydatetime()
    seq = input("Is this a test run? Yes or No? ")
    if seq.lower() == "yes":
        for t in tlist:
            if len(fields) > 0:
                make_combplots(figprepath, controlprepath, nlhprepath, fields, t)
            if plotdiv:
                make_divplots(figprepath, controlprepath, nlhprepath, t)
    elif seq.lower() == "no":
        ppool = ProcessPoolExecutor(max_workers = 4) 
        if len(fields) > 0:
            ppool.map(partial_combplots, tlist)
        if plotdiv:
            ppool.map(partial_divplots, tlist)
        ppool.shutdown()
    et = perf_counter()
    print(f"Plotting took {et-st:.2f} seconds")
    figpath = f"{figprepath}/combplots"
    if not path.exists(figpath):
        raise Exception("combplots folder not found in the specified subdirectory!")
    st = perf_counter()
    reload(mkv)
    if plotdiv == True:
        fields.append("div")
    for field in fields:
        mkv.makevidcomb(figpath, field)
    et = perf_counter()
    print(f"Animation took {et-st:.2f} seconds")
elif runtype.lower() == "quicktest":
    print("quicktest")
    t0 = "2010-01-02-180000"
    t1 = "2010-01-02-181000"
    tlist = date_range(t0, t1, freq = timedelta(minutes = 20)).to_pydatetime()
    fields = ["w","SnowRate"]
    for t in tlist:
        make_divplots(t)
        make_combplots(fields, t)
    print("done")
else:
    raise Exception("Runtype must be either plot, animate, or both!")
