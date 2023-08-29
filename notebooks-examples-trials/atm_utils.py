# -*- coding: utf-8 -*-
# isort: off

import cfgrib
from datetime import datetime
from ecmwf.opendata import Client
from ipyleaflet import Choropleth, Map, basemap_to_tiles, LayersControl, LegendControl
import ipywidgets as widgets
from IPython.display import display
from localtileserver import get_leaflet_tile_layer, TileClient
# from Magics import macro as magics
# from Magics.macro import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import rasterio
# import requests
import rioxarray as rxr
import xarray as xr

#%% Download

def dwnl_atmdata_step(variables, stepsdict):
    """
    Function which downloads 1, 2, 5 and 10 days ahead forecasts from today
    """
    fnames = []
    today = datetime.today().strftime('%Y%m%d')
    for var in variables:
        if var == "10fgg15": steps = stepsdict["10fgg15"]
        else: steps = stepsdict["base"]
        for s in steps:
            if var == "10fgg15":
                rqt = {
                    "date": 0,      #date start of the forecast: today
                    "time": 0,      #time start of the forecast, can be 0 or 12
                    "step": s,      #step of the forecast: 1, 2, 5, 10 days
                    "stream": "enfo",
                    "type": "ep",
                    "param": var,
                }
            else:
                rqt = {
                    "date": 0,      #date start of the forecast: today
                    "time": 0,      #time start of the forecast, can be 0 or 12
                    "step": s,      #step of the forecast: 1, 2, 5, 10 days
                    "stream": "oper",
                    "type": "fc",
                    "levtype": "sfc",
                    "param": var,
                }
            filename = f"data/atm/{var}_{rqt['date']}_{today}_time{rqt['time']}_step{rqt['step']}_{rqt['stream']}_{rqt['type']}.grib"
            if not os.path.exists(filename):
                client = Client(source = "ecmwf", beta = True)
                client.retrieve(
                    request = rqt,
                    target = filename
                )
            fnames.append(filename)
    return(fnames)

#%% Load

def gen_raster(var, filename):
    """
    Generates .tiff raster files from .grib files downloaded from ECMWF's Open Data through dwnl_atmdata

    var: str
        Code of the variable
    d: str
        Date of the downloaded variable file

    Returns:
    tiffpath: str
        Path to the created .tiff file
    """
    tiffpath = ".".join([filename.split(".")[0], "tiff"])
    if not os.path.exists(tiffpath):
        print(var, " conversion")
        f = xr.load_dataset(filename, engine = "cfgrib")
        f = f.rio.write_crs("epsg:4326")
        f.rio.to_raster(tiffpath)
        f.close()
    return(tiffpath)

def load_atmdata(varlst, fnames):
    """
    varlst: str, list of str
        List containing strings of the variables to load
    fnames: str, list of str
        List containing the paths to downloaded data
    
    Returns:
    vardict: dict
        A dictionary of rasterio DatasetReader objects, each assigned to the corresponding variable
    """
    vardict = {}
    for i, var in enumerate(varlst):
        lst = []
        tool = [x for x in fnames if var in x] #filenames containing the variable
        for filename in tool:
            lst.append(rasterio.open(gen_raster(var, filename)))
        vardict[f"{var}"] = lst
        print(var, " loaded")
    return(vardict)

#%% Plot

def get_colordict(array, cmap):
    """
    array: numpy.array
        Array containing the values of which the quantiles will be obtained
    cmap: str
        Name of the matplotlib colormap associated to the array

    Returns:
    dict
        Dictionary containing Hex codes as keys and the array quantiles as values
    """
    n_col = 8
    cmap = plt.get_cmap(cmap, n_col)
    custom_palette = [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
    quantiles = [round(np.quantile(array, q)) for q in np.linspace(0,1,n_col)]
    return(dict(zip(quantiles, custom_palette)))

def plot_atmdata_step(vardict, step, coord, stepsdict):
    """
    vardict: dict
        Dictionary containing the rasters uploaded trhough load_atmdata
    step: int
        Index of the desired step inside the list defined in stepsdict
    coord: list or tuple
        Coordinates of the map central point. Provide them as lat, lon

    Returns:
    m: ipyleaflet.Map
    """
    namedict = {
        "msl": "Mean sea level pressure [Pa]",
        "2t": "2 meter temperature [K]",
        "tp": "Total Precipitation [m]",
        "10fgg15": "10 metre wind gust of at least 15 m/s [%]",
    }
    palettedict = {
        "msl": "cool",
        "2t": "RdBu_r",
        "tp": "PuBu",
        "10fgg15": "winter",
    }
    m = Map(center=(coord[0], coord[1]), zoom = 3)
    for var in vardict.keys():
        if var == "10fgg15": steps = stepsdict["10fgg15"]
        else: steps = stepsdict["base"]
        r = [x for x in vardict[var] if f"step{steps[step]}" in x.name][0]
        print("Plotting ", namedict[var])
        client = TileClient(r)
        t = get_leaflet_tile_layer(client, name = namedict[var], opacity = 0.7, palette = palettedict[var], n_colors = 8)
        m.add_layer(t)

        # find a way to add a colormap
        # in leaflet you can do it as:
        # m.add_colorbar(colors=custom_palette, vmin=r.read(1).ravel().min(), vmax=r.read(1).ravel().max())
        #
        # doing as below works but the legends created are too big:
        # colordict = get_colordict(r.read(1).ravel(), palettedict[var])
        # legend = LegendControl(legend = colordict, name = namedict[var], position="topright")
        # # legend = LegendControl({"low":"#24dbff", "medium":"#A55", "High":"#500"}, name = namedict[var], position="bottomleft")
        # m.add_control(legend)
    m.add_control(LayersControl())
    m.layout.height = "700px"
    return(m)

#%% Kept for compatibility

def dwnl_atmdata(variables, dates):
    """
    Function which replicates the code written by Laura
    It downloads 10 day ahead forecasts from the dates specified in dates
    It need the dates from which the forecast starts
    """
    fnames = []
    today = datetime.today().strftime('%Y%m%d')
    for var in variables:            
        for d in dates:
            if var == "10fgg15":
                rqt = {
                    "date": str(d), #date start of the forecast
                    "time": 0,      #time start of the forecast, can be 0 or 12
                    "step": "216-240",    #step of the forecast: 10 days
                    "stream": "enfo",
                    "type": "ep",
                    "param": var,
                }
            else:
                rqt = {
                    "date": str(d), #date start of the forecast
                    "time": 0,      #time start of the forecast, can be 0 or 12
                    "step": 240,    #step of the forecast: 10 days
                    "stream": "oper",
                    "type": "fc",
                    "levtype": "sfc",
                    "param": var,
                }
            filename = f"data/atm/{var}_{rqt['date']}_{today}_time{rqt['time']}_step{rqt['step']}_{rqt['stream']}_{rqt['type']}.grib"
            if not os.path.exists(filename):
                client = Client(source = "ecmwf", beta = True)
                client.retrieve(
                    request = rqt,
                    target = filename
                )
            fnames.append(filename)
    return(fnames)

def plot_atmdata_date(vardict, date, coord):
    """
    vardict: dict
        Dictionary containing the rasters uploaded trhough load_atmdata
    date: str, list of str
        Dates considered
    coord: list or tuple
        Coordinates of the map central point. Provide them as lat, lon

    Returns:
    m: ipyleaflet.Map
    """
    namedict = {
        "msl": "Mean sea level pressure [Pa]", #meansealevelpressure
        "2t": "2 meter temperature [K]",
        "tp": "Total Precipitation [m]",
        "10fgg15": "10 metre wind gust of at least 15 m/s [%]",
    }
    palettedict = {
        "msl": "cool", #meansealevelpressure
        "2t": "RdBu_r",
        "tp": "PuBu",
        "10fgg15": "winter",
    }
    m = Map(center=(coord[0], coord[1]), zoom = 3)
    for var in vardict.keys():
        r = [x for x in vardict[var] if date in x.name][0]
        print("Plotting ", namedict[var])
        client = TileClient(r)
        t = get_leaflet_tile_layer(client, name = namedict[var], opacity = 0.7, palette = palettedict[var])
        m.add_layer(t)
    m.add_control(LayersControl())
    m.layout.height = "700px"
    return(m)