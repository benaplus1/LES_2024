# LES_2024

Instructions for recreating the figures seen in the manuscript "Lake Huron Enhances Snowfall Downwind of Lake Erie: A Modeling Study of the 2010 New Year's Lake Effect Snowfall Event."

For the "sim_setup_hires_mod.py" file, which you will use to alter the mksfc and varfiles to get the NLH, MORESNOW, and VARTEMP simulations ready, prompts when executing the python script will guide you through the process.

For "PaperPlots.ipynb", just run all the cells. The only thing you might need to change is in the third cell, when the font paths for "Liberation Serif" are installed. You can also delete that if you don't care about fonts.

For "les_combplots_paper_carto.py", executing the script will first lead to a prompt of whether you want to plot, animate, or do both. Just select "plot" for now, which will make the individual frames of plots. It will then ask for the start and end times of plotting. You need to enter these in YYYY-mm-dd-HHMMSS format. If you want to plot the whole simulation (ignoring spin-up), then put "2010-01-02-120000" for the start time, and "2010-01-03-120000" for the end time. It will then ask you for which fields you want to plot. The field names are currently hardcoded with the following options: "Div" for 350 and 500 m horizontal divergence, "SnowRate" for snowfall rate, "w" for 700 m vertical velocity, "VapMix" for 700 m water vapor mixing ratio, "ShFlux" for surface sensible heat flux, "LhFlux" for surface latent heat flux, and "CldTop" for cloud top heights. Note that the altitudes for horizontal divergence are customizable as "divalt1" and "divalt2" variables in the code, as are the altitudes for w ("walT") and vapor mixing ratio ("vapalt").

If instead of just plotting frames you want to animate them, select "animate" and then enter the path to where the frames are stored. By default this is hardcoded to attempt to plot all fields possible ("snowrate_carto" for SnowRate, "lhcomp_carto" for LhFlux, "shcomp_carto" for ShFlux, "wcomp_carto" for w, "vapcomp_carto" for VapMix, and "cldtop_carto" for CldTop). It will print some error messages if some of the fields it is looking for haven't been plotted yet, but it shouldn't crash. The "animate" function of this script uses a bash script originally developed by Sean Freeman while he was a member of the van den Heever group.

Finally, run "aensitivity_analysis.ipynb". Make sure to change file paths, and then it should work fine.
