# -*- coding: utf-8 -*-
# isort: off

import branca.colormap as bc
import cfgrib
from datetime import datetime, timedelta
from ecmwf.opendata import Client
from ipyleaflet import Map, ColormapControl, LayersControl, basemaps, FullScreenControl
from ipyleaflet.velocity import Velocity
from IPython.display import display
from localtileserver import get_leaflet_tile_layer, TileClient
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import rasterio
import rioxarray as rxr
import xarray as xr

#%% Download

def dwnl_atmdata(variables, stepsdict, stdate = 0, source = "azure", pr = True):
    """
    Function which downloads forecasts from today for the provided steps

    variables: str, list of str
        List containing the codes of the variables to be downloaded.
        Codes can be found here: https://github.com/ecmwf/ecmwf-opendata/tree/main#parameters-and-levels
    stepsdict: dict
        Dictionary containing the steps codes needed for each variable. The standard step format is under "base".
    stdate: int (0/1/-1) or str (YYYYMMDD), optional
        Integer or string defining the starting date to download the data from. Default is 0, meaning
        the function will download the forecasts from the current day.
    source: str, optional
        Parameter needed for ecmwf.opendata Client class. Can be "azure" or "ecmwf". Default: azure
    pr: bool, optional
        If True, prints a message. Default is True

    Returns:
    fnames: list of str
        List containing the paths to the downloaded data
    """
    fnames = []
    if stdate == 0:
        stdate = datetime.today().strftime('%Y%m%d')
    else:
        stdate = stdate.strftime('%Y%m%d')
    if not os.path.exists(f'data/atm/{stdate}'):
        os.mkdir(f'data/atm/{stdate}')
    for var in variables:
        if var == "10fgg25":
            steps = stepsdict["10fgg25"]
        else:
            steps = stepsdict["base"]
        if var == "10fgg25":
                    rqt = {
                        "date": stdate,     #date start of the forecast
                        "time": 0,          #time start of the forecast, can be 0 or 12
                        "step": steps,      #step of the forecast: 1, 2, 5, 10 days
                        "stream": "enfo",
                        "type": "ep",
                        "param": var,
                    }
        elif var == "wind":
            rqt = {
                "date": stdate,     #date start of the forecast
                "time": 0,          #time start of the forecast, can be 0 or 12
                "step": steps,      #step of the forecast: 1, 2, 5, 10 days
                "stream": "oper",
                "type": "fc",
                "levtype": "sfc",
                "param": ["10u", "10v"],
            }
        else:
            rqt = {
                "date": stdate,     #date start of the forecast
                "time": 0,          #time start of the forecast, can be 0 or 12
                "step": steps,      #step of the forecast: 1, 2, 5, 10 days
                "stream": "oper",
                "type": "fc",
                "levtype": "sfc",
                "param": var,
            }
        client = Client(source = source, preserve_request_order = True)
        try:
            filename = f"data/atm/{rqt['date']}/{var}_{rqt['date']}_time{rqt['time']}_steps{rqt['step'][0]}-{rqt['step'][-1]}_{rqt['stream']}_{rqt['type']}.grib"
            if not os.path.exists(filename):
                client.retrieve(
                    request = rqt,
                    target = filename
                );
        except:
            # Usually early in the morning the forecast of the current day is not available
            # > get the forecast of the day before
            rqt['date'] = datetime.strptime(rqt['date'], '%Y%m%d') - timedelta(days = 1)
            tooldate = rqt['date'].strftime('%Y%m%d')
            filename = f"data/atm/{tooldate}/{var}_{tooldate}_steps{rqt['step'][0]}-{rqt['step'][-1]}_{rqt['stream']}_{rqt['type']}.grib"
            if not os.path.exists(filename):
                if not os.path.exists(f"data/atm/{tooldate}"):
                    os.mkdir(f"data/atm/{tooldate}")
                client.retrieve(
                    request = rqt,
                    target = filename
                );
            if pr: print(f"Today's forecast not available, downloaded yesterday's forecast: {rqt['date'].strftime('%d/%b/%Y')}")       
        fnames.append(filename)
    return(fnames)

#%% Load

def load_atmdata(varlst, fnames, pr = False):
    """
    Loads the downloaded files

    varlst: str, list of str
        List containing strings of the variables to load
    fnames: str, list of str
        List containing the paths to downloaded data
    pr: bool, optional
        If True, prints the variables which are being loaded. Default is False
    
    Returns:
    vardict: dict
        A dictionary of rasterio DatasetReader objects and xarray.Dataset objects,
        each assigned to the corresponding variable, or a dictionary containing the
        geenrated files from gen_raster
    """
    vardict = {}
    for var in varlst:
        tool = [x for x in fnames if var in x] #filename containing the variable
        vardict[f"{var}"] = gen_raster(var, tool[0])
        if pr: print(var, " loaded")
    return(vardict)

def gen_raster(var, filename, delete = False, pr = False):
    """
    Generates .tiff raster files from .grib files downloaded from ECMWF's Open Data through dwnl_atmdata

    var: str
        Code of the variable
    filename: str
        Path of the .grib file to be converted to .tiff
    delete: bool, optional
        If True, it will delete the original .grib file once the .tiff file is created. Default: False
    pr: bool, optional
        If True, prints the variables which are being converted. Default is False

    Returns:
    tiffpath: str
        Path to the created .tiff file
    """
    fix = "_".join(filename.split("_")[0:2] + filename.split("_")[3:])
    if var != "wind":
        fix = ".".join([fix.split(".")[0], "tiff"])
    else:
        fix = ".".join([fix.split(".")[0], "nc"])
    if pr: print(var, " conversion")
    f = xr.load_dataset(filename, engine = "cfgrib")
    pathlst = []
    for step in range(0, len(f.step)):
        outpath = ".".join(["_".join([fix.split(".")[0], f"step{step}"]), fix.split(".")[1]])
        if not os.path.exists(outpath):
            outf = f.isel(step = step)
            if var == "msl":
                outf["msl"] = outf.msl/100 #hPa
            if var == "skt":
                outf["skt"] = outf["skt"] - 273.15 #Celsius degrees
            outf = outf.rio.write_crs("epsg:4326")
            if var != "wind":
                outf = outf.rio.reproject("epsg:3857") #to avoid errors when creating layers with get_leaflet_tile_layer
                outf.rio.to_raster(outpath)
            else:
                outf.to_netcdf(outpath)
        pathlst.append(outpath)
    f.close()
    os.remove(f"{filename}.923a8.idx")
    if delete:
        os.remove(filename)
    return(pathlst)

#%% Plot

def plot_atmdata_step(vardict, step, coord, stepsdict):
    """
    Plots the atmospheric variables for a specified forecast step

    vardict: dict
        Dictionary containing the path to the rasters, obtained through load_atmdata.
        Returned by load_atmdata
    step: str
        String returned from the widget used to select the forecast period
    coord: tuple
        Coordinates of the map central point. Provide them as lat, lon
    stepsdict: dict
        Dictionary containing the steps codes needed for each variable. The standard step format is under "base".
        
    Returns:
    None
    """
    namedict = {
        "msl": "Mean sea level pressure [hPa]",
        "skt": "Skin temperature [°C]",
        "tp": "Total Precipitation [m]",
        "10fgg25": "Probability of 10 metre wind gust of at least 25 m/s [%]",
        "wind": "10 metre wind component [m/s]",
    }
    m = Map(basemap = basemaps.Esri.WorldTopoMap, center = coord, zoom = 3)
    for var in vardict.keys():
        palette = get_palette(var)
        s = sel_forecast(var, step, stepsdict)
        r = [x for x in vardict[var] if f"step{s}" in x][0] #extract the correct raster path
        if var != "wind":
            client = TileClient(r)
            t = get_leaflet_tile_layer(client, name = namedict[var], opacity = 0.7, palette = palette, max_zoom = 30)
            m.add_layer(t)
            with rasterio.open(r) as raster:
                minv = "%.2f" % round(raster.read(1).ravel().min(), 1)
                maxv = "%.2f" % round(raster.read(1).ravel().max(), 1)
            # add colorbar
            cmap_control = ColormapControl(
                                            caption = namedict[var],
                                            colormap = bc.StepColormap(palette),
                                            value_min = float(minv),
                                            value_max = float(maxv),
                                            position = 'topright',
                                            transparent_bg = True
                                            )
            m.add(cmap_control)
        else:
            pass
            display_options = {
                'velocityType': 'Global Wind',
                'displayPosition': 'bottomleft',
                'displayEmptyString': 'No wind data'
            }
            with xr.load_dataset(r) as netcdf:
                wind_layer = Velocity(  name = namedict[var],
                                        data = netcdf,
                                        zonal_speed = 'u10',
                                        meridional_speed = 'v10',
                                        latitude_dimension = 'latitude',
                                        longitude_dimension = 'longitude',
                                        velocity_scale = 0.01,
                                        max_velocity = float(max(netcdf["v10"].values.ravel().max(), netcdf["u10"].values.ravel().max())),
                                        display_options = display_options,
                                        color_scale = palette
                                     )
            m.add_layer(wind_layer)
    m.add_control(LayersControl())
    m.add_control(FullScreenControl())
    m.layout.height = "700px"
    display(m)

def sel_forecast(var, str, stepsdict):
    """
    Returns the step used to compose the filename of the specified variable

    str: str
        String returned from the widget used to select the forecast period
    stepsdict: dict
        Dictionary containing the steps codes needed for each variable. The standard step format is under "base".
    
    Returns:
    out: int or str
        Forecast step, from 0 to 3. Used to extract the selected timeframe from the downloaded file
    """
    toollist = [f"{i}h from selected date" for i in stepsdict["base"]]
    out = toollist.index(str)
    if var == "wind" and out != 0:
        out = out - 1
    return(out)

def get_palette(var):
    """
    Returns hex codes for colors in the palette corresponding to the variable requested.
    Palettes' RGB values are retrieved from the IPCC guide:
    - https://www.ipcc.ch/site/assets/uploads/2022/09/IPCC_AR6_WGI_VisualStyleGuide_2022.pdf
    - https://github.com/IPCC-WG1/colormaps

    var: str
        Code of the variable for which to return the palette
    
    Returns:
    hex: list of str
        List of 10 hex codes corresponding to 
    """
    def rgb_to_hex(rgb):
        return '#{:02x}{:02x}{:02x}'.format(rgb[0], rgb[1], rgb[2])
    palettedict = {
        "msl": [(230, 240, 240),(182, 217, 228),(142, 192, 226),(118, 163, 228),(116, 130, 222),(121, 97, 199),(118, 66, 164),(107, 40, 121),(86, 22, 75),(54,14, 36)],
        "skt": [(254, 254, 203),(251, 235, 153),(244, 204, 104),(235, 167, 84),(228, 134, 80),(209, 98, 76),(164, 70, 66),(114, 55, 46),(66, 40, 24),(25, 25, 0)],
        "tp": [(255, 255, 229),(217, 235, 213),(180, 216, 197),(142, 197, 181),(105, 177, 165),(67, 158, 149),(44, 135, 127),(29, 110, 100),(14, 85, 74),(0, 60, 48)],
        "10fgg25": [(254, 252, 205),(235, 219, 144),(210, 192, 83),(170, 171, 32),(121, 154, 5),(70, 136, 22),(24, 114, 39),(13, 88, 44),(24, 61, 36),(23, 35, 18)],
        "wind": [(254, 252, 205),(235, 219, 144),(210, 192, 83),(170, 171, 32),(121, 154, 5),(70, 136, 22),(24, 114, 39),(13, 88, 44),(24, 61, 36),(23, 35, 18)]
    }
    hex = [rgb_to_hex(x) for x in palettedict[var]]
    return(hex)

# %% Additional functions

def reverse_branca(brancacmap):
    """
    Reverses branca.colormap.linear instances
    """
    out = bc.LinearColormap(colors = list(reversed(brancacmap.colors)))
    return(out)

def get_colordict(array, cmap, var = None, n_col = 8):
    """
    Creates a dictionary containing values and cmap colors associated to an array's
    values quantiles

    array: numpy.array
        Array containing the values of which the quantiles will be obtained
    cmap: str
        Name of the matplotlib colormap associated to the array
    var: str, optional
        Name of the variable. Needed for "tp"
    n_col: int, optional
        Number of colors to be extracted from cmap

    Returns:
    dict
        Dictionary containing Hex codes as keys and the array quantiles as values
    """
    cmap = plt.get_cmap(cmap, n_col)
    custom_palette = [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
    if var == "tp":
        quantiles = [round(np.quantile(array, q), 3) for q in np.linspace(0,1,n_col)]
    else:
        quantiles = [round(np.quantile(array, q), 1) for q in np.linspace(0,1,n_col)]
    return(dict(zip(quantiles, custom_palette)))

def get_colorbar(colordict, title, var):
    """
    Create a plot with the colormap to be added to the ipyleaflet map

    colordict: dict
        Output of get_colordict
    title: str
        Title to be assigned to the colormap
    var: str
        Name of the variable for which the colorbar is being computed. Needed to
        create the output path
    
    Returns:
    imgpath: str
        Path to the colorbar plot, saved as .png
    """
    squares = [(i, 1) for i in range(len(colordict.keys()))]
    colors = [colordict[i] for i in colordict.keys()]
    fig, ax = plt.subplots(1,1, )
    plt.broken_barh(squares, (0,1), facecolors = colors)
    # set x labels
    ax.set_xticks([x+0.5 for x in range(len(colordict.keys()))])
    ax.set_xticklabels(colordict.keys())
    ax.tick_params('x', length = 0)
    for tick in ax.get_xticklabels():
        tick.set_fontname("tahoma")
    # remove bounding box and y ticks
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_ticks([])
    # set aspect ration
    ax.set_aspect(0.5)
    # set title
    ax.set_title(title, fontname = "tahoma")
    # save and return path
    imgpath = f"data/colorbars/colorbar_{var}.png"
    fig.savefig(imgpath, bbox_inches = "tight")
    plt.close(fig)
    return(imgpath)
