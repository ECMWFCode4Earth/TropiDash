
# -*- coding: utf-8 -*-
# isort: off

#Necessary packages
import branca.colormap as bc
from ipyleaflet import Choropleth, Map, LayersControl, ColormapControl, basemaps, FullScreenControl
from IPython.display import display
import ipywidgets as widgets
from localtileserver import get_leaflet_tile_layer, TileClient
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import rasterio
from rasterio.enums import Resampling
import requests
import json

# %% General functions

def rgb_to_hex(rgb):
    """
    Returns the HEX code of the RGB tuple provided

    rgb: tuple of int
        rgb[0]: red, rgb[1]: green, rgb[2]: blue
    """
    return '#{:02x}{:02x}{:02x}'.format(rgb[0], rgb[1], rgb[2])

def get_colorlist(cmap, n_col = 10):
    """
    Creates a list containing cmap colors in hex format

    cmap: str
        Name of the matplotlib colormap associated to the array
    n_col: int, optional
        Number of colors to be extracted from cmap. Default is 10

    Returns:
    list
        List containing Hex codes (str)
    """
    cmap = plt.get_cmap(cmap, n_col)
    return([mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)])

# %% Coastal hazard functions

def dwnl_coastalhaz(rp):
    """
    Downloads coastal hazard raster maps from The World Bank Data Catalog service
    by providing the hazard return period.
    Data source: https://datacatalog.worldbank.org/search/dataset/0038579/Global-coastal-flood-hazard
    
    rp: str
        Return period.
        Can be one of the following: 5yr, 10yr, 50yr, 100yr, 250yr, 500yr, 1000yr.

    Returns:
    None    
    """
    links = {
            "5yr": "https://www.geonode-gfdrrlab.org/uploaded/layers/2019/10/16/ss_muis_rp0005m.tif",
            "10yr": "https://www.geonode-gfdrrlab.org/uploaded/layers/2019/10/16/ss_muis_rp0010m.tif",
            "50yr": "https://www.geonode-gfdrrlab.org/uploaded/layers/2019/10/16/ss_muis_rp0050m.tif",
            "100yr": "https://www.geonode-gfdrrlab.org/uploaded/layers/2019/10/16/ss_muis_rp0100m.tif",
            "250yr": "https://www.geonode-gfdrrlab.org/uploaded/layers/2019/10/16/ss_muis_rp0250m.tif",
            "500yr": "https://www.geonode-gfdrrlab.org/uploaded/layers/2019/10/16/ss_muis_rp0500m.tif",
            "1000yr": "https://www.geonode-gfdrrlab.org/uploaded/layers/2019/10/16/ss_muis_rp1000m.tif",
            }
    filename = f"data/impacts/coastalhaz_{rp}.tif"
    if not os.path.exists(filename):
        tif = requests.get(links[rp])
        with open(filename, "wb") as tiffile:
                tiffile.write(tif.content)
                # print("Coastal hazard data - Download complete")

def load_coastalhaz(rp, resample = False, open = True):
    """
    Loads coastal hazard raster maps downloaded through dwnl_coastalhaz
    by providing the return period. Returns a path
    
    rp: str
        Return period.
        Can be one of the following: 5yr, 10yr, 50yr, 100yr, 250yr, 500yr, 1000yr.
    resample: bool, optional
        If True, the raster is resampled at a lower spatial resolution. Default is False
    open: bool, optional
        If True, it will open the file and return a rasterio.DatasetReader object
        otherwise, it will return the path to the file. Default is True.
    
    Returns:
    coh: str or rasterio.DatasetReader object
        Path to the coastal hazard .tif file
    """
    path = f"data/impacts/coastalhaz_{rp}.tif"
    if resample:
        print("Resampling coastal hazard layer to ease its plotting")
        coh = resample_raster(path) #returns the path to resampled file
        if open:
            coh = rasterio.open(coh)
    elif open:
        coh = rasterio.open(path)
    else:
        coh = path
    return(coh)

def plot_coastalhaz(coh, rp, addlayer = True, coord = None, m = None):
    """
    Plot coastal hazard raster provided as a layer in an ipyleaflet map

    coh: rasterio.DatasetReader object
        rasterio.DatasetReader object of the coastal hazard .tiff file to be plotted
    rp: str
        Return period associated to coh. Can be one of the following: 5yr, 10yr, 50yr, 100yr, 250yr, 500yr, 1000yr.
    addlayer: bool, optional
        If True, the function will only add a layer to the ipyleaflet.Map provided
        through m. If False, it will create a new ipyleaflet.Map. Default is True
    coord: tuple of float, optional
        Tuple containing latitude and longitude where to center the map. Only needed if 
        addlayer = False. Default is None
    m: ipyleaflet.Map, optional
        Only if addlayer is True, m needs to be provided as an ipyleaflet.Map object. Default is None
    
    Returns:
    m: ipyleaflet.Map
    """
    palette = [(8, 29, 88),(34, 84, 163),(49, 165, 194),(151, 214, 184),(239, 248, 182),(255, 239, 165),(254, 191, 90),(252, 113, 52),(218, 20, 30),(128, 0, 38)]
    palettehex = [rgb_to_hex(x) for x in palette]
    client = TileClient(coh)
    t = get_leaflet_tile_layer(client, name = f"Coastal hazard - RP: {rp}", opacity = 0.7, palette = palettehex)
    minv = "%.2f" % round(coh.read(1).ravel().min(), 1)
    maxv = "%.2f" % round(coh.read(1).ravel().max(), 1)
    if addlayer:
        cmap_control = ColormapControl(
                                        caption = f"Coastal hazard - RP: {rp}",
                                        colormap = bc.StepColormap(palettehex),
                                        value_min = float(minv),
                                        value_max = float(maxv),
                                        position = 'topright',
                                        transparent_bg = True
                                        )
        m.add(cmap_control)
        m.add_layer(t)
    else:
        m = Map(center = coord, zoom = 3)
        m.add_layer(t)
        m.add_control(LayersControl())
        m.layout.height="700px"
    return(m)

# %% Cyclone Hazard functions

def dwnl_cyclonehaz(rp):
    """
    Downloads cyclone hazard raster maps from The World Bank Data Catalog service
    by providing the hazard return period
    Data source: https://datacatalog.worldbank.org/search/dataset/0038577/Global-cyclone-hazard

    rp: str
        Return period.
        Can be one of the following: 50yr, 100yr, 250yr, 500yr, 1000yr.
    
    Returns:
    None
    """
    links = {
            "50yr": "https://www.geonode-gfdrrlab.org/uploaded/layers/2019/02/26/viento_mundo_tr50_int1_Ai5zJup.tif",
            "100yr": "https://www.geonode-gfdrrlab.org/uploaded/layers/2019/02/26/viento_mundo_tr100_int1_afg5XsL.tif",
            "250yr": "https://www.geonode-gfdrrlab.org/uploaded/layers/2019/02/26/viento_mundo_tr250_int1.tif",
            "500yr": "https://www.geonode-gfdrrlab.org/uploaded/layers/2019/02/26/viento_mundo_tr500_int1.tif",
            "1000yr": "https://www.geonode-gfdrrlab.org/uploaded/layers/2019/02/26/viento_mundo_tr1000_int1.tif",
            }
    filename = f"data/impacts/cyclonehaz_{rp}.tif"
    if not os.path.exists(filename):
        tif = requests.get(links[rp])
        with open(filename, "wb") as tiffile:
                tiffile.write(tif.content)
                # print("Cyclone hazard download complete")

def load_cyclonehaz(rp, open = True):
    """
    Loads coastal hazard raster maps downloaded through dwnl_coastalhaz
    by providing the return period.

    rp: str
        Return period.
        Can be one of the following: 50yr, 100yr, 250yr, 500yr, 1000yr.
    open: bool, optional
        If True, it will open the file and return a rasterio.DatasetReader object
        otherwise, it will return the path to the file. Default is True.
    
    Returns:
    cyh: str or rasterio.DatasetReader object
    """
    path = f"data/impacts/cyclonehaz_{rp}.tif"
    if open:
        cyh = rasterio.open(path)
    else:
        cyh = path
    return(cyh)

def plot_cyclonehaz(cyh, rp, addlayer = True, coord = None, m = None):
    """
    Plots the cyclone hazard layer provided as input in an ipyleaflet map

    cyh: rasterio.DatasetReader object
        rasterio.DatasetReader object of the cyclone hazard .tiff file to be plotted
    rp: str
        Return period associated to cyh. Can be one of the following: 50yr, 100yr, 250yr, 500yr, 1000yr.
    addlayer: bool, optional
        If True, the function will only add a layer to the ipyleaflet.Map provided
        through m. If False, it will create a new ipyleaflet.Map. Default is True
    coord: tuple of float, optional
        Tuple containing latitude and longitude where to center the map. Only needed if 
        addlayer = False. Default is None
    m: ipyleaflet.Map, optional
        Only if addlayer is True, m needs to be provided as an ipyleaflet.Map object. Default is None
    
    Returns:
    m: ipyleaflet.Map
    """
    palettehex = get_colorlist("magma_r")
    client = TileClient(cyh)
    t = get_leaflet_tile_layer(client, name = f"Cyclone hazard - RP: {rp}", opacity = 0.7, palette = palettehex)
    minv = "%.2f" % round(cyh.read(1).ravel().min(), 1)
    maxv = "%.2f" % round(cyh.read(1).ravel().max(), 1)
    if addlayer:
        cmap_control = ColormapControl(
                                        caption = f"Cyclone hazard - RP: {rp}",
                                        colormap = bc.StepColormap(palettehex),
                                        value_min = float(minv),
                                        value_max = float(maxv),
                                        position = 'topright',
                                        transparent_bg = True
                                        )
        m.add(cmap_control)
        m.add_layer(t)
    else:
        m = Map(center = coord, zoom = 3)
        cmap_control = ColormapControl(
                                        caption = f"Cyclone hazard - RP: {rp}",
                                        colormap = bc.StepColormap(palettehex),
                                        value_min = float(minv),
                                        value_max = float(maxv),
                                        position = 'topright',
                                        transparent_bg = True
                                        )
        m.add(cmap_control)
        m.add_layer(t)
        m.add_control(LayersControl())
        m.layout.height="700px"
    return(m)

#%% Population functions

def load_poplayer(path = None, returnpath = True):
    """
    Loads the population raster with 10 km spatial resolution saved in data/impacts folder
    The raster file was downloaded from https://hub.worldpop.org/geodata/summary?id=24777 with 1 km resolution
    Then, it was resampled to 10 km resolution through QGIS
    
    path: str, optional
        Path to the population raster file. Default is None
    returnpath: bool, optional
        If True, returns the path to the file. If False, a rasterio.DatasetReader object. Deafult is True

    Returns:
    r: rasterio.DatasetReader object or str
    """
    if path is None:
        path = "data/impacts/ppp_2020_resampled_10km_log.tif"
    if returnpath:
        r = path
    else:
        r = rasterio.open(path)
    return(r)

def plot_poplayer(addlayer = True, coord = None, m = None):
    """
    Plots the population layer with 10km spatial resolution saved in data/impacts folder

    addlayer: bool, optional
        If True, the function will only add a layer to the ipyleaflet.Map provided
        through m. If False, it will create a new ipyleaflet.Map. Default is True
    coord: tuple of float, optional
        Tuple containing latitude and longitude where to center the map. Only needed if 
        addlayer = False. Default is None
    m: ipyleaflet.Map, optional
        Only if addlayer is True, m needs to be provided as an ipyleaflet.Map object. Default is None
    
    Returns:
    m: ipyleaflet.Map
    """
    palette = [(255, 255, 229),(217, 235, 213),(180, 216, 197),(142, 197, 181),(105, 177, 165),(67, 158, 149),(44, 135, 127),(29, 110, 100),(14, 85, 74),(0, 60, 48)]
    palettehex = [rgb_to_hex(x) for x in palette]
    with rasterio.open(load_poplayer(returnpath = True)) as r:
        minv = "%.2f" % 0
        maxv = "%.2f" % round(np.exp(r.read(1).ravel().max())/100, 1) #approximate density (pop/km2)
        if addlayer:
            client = TileClient(r)
            t = get_leaflet_tile_layer(client, name = "Population density - n° of people/sq. km (2020)", opacity = 0.7, palette = palettehex, nodata = r.nodata)
            cmap_control = ColormapControl(
                                            caption = "Population density - n° of people/sq. km (2020)",
                                            colormap = bc.StepColormap(palettehex),
                                            value_min = float(minv),
                                            value_max = float(maxv),
                                            position = 'topright',
                                            transparent_bg = True
                                            )
            m.add(cmap_control)
            m.add_layer(t)
        else:
            m = Map(center = coord, zoom = 3)
            client = TileClient(r)
            t = get_leaflet_tile_layer(client, name = "Population density - n° of people/sq. km (2020)", opacity = 0.7, palette = palettehex, nodata = r.nodata)
            cmap_control = ColormapControl(
                                            caption = "Population density - n° of people/sq. km (2020)",
                                            colormap = bc.StepColormap(palettehex),
                                            value_min = float(minv),
                                            value_max = float(maxv),
                                            position = 'topright',
                                            transparent_bg = True
                                            )
            m.add(cmap_control)
            m.add_layer(t)
            m.add_control(LayersControl())
            m.layout.height = "700px"
    return(m)

# %% World Risk Index functions

def dwnl_riskidx():
    """
    Downloads the World Risk Index data for 2022

    Data source: https://data.humdata.org/dataset/1efb6ee7-051a-440f-a2cf-e652fecccf73
    """
    links = {
     #     "meta": "https://data.humdata.org/dataset/1efb6ee7-051a-440f-a2cf-e652fecccf73/resource/9ab4cdc6-9682-45a6-8314-db5e0d49e7a7/download/worldriskindex-meta.xlsx",
         "data": "https://data.humdata.org/dataset/1efb6ee7-051a-440f-a2cf-e652fecccf73/resource/1b47b40c-f746-427c-bbf8-8caa157e03da/download/worldriskindex-2022.csv",
     #     "trend": "https://data.humdata.org/dataset/1efb6ee7-051a-440f-a2cf-e652fecccf73/resource/9ae68502-6011-4276-a99d-118bb2826323/download/worldriskindex-trend.csv"
    }
    for id in links:
        f = requests.get(links[id])
        if id == "meta":
             fformat = "xlsx"
        else:
             fformat = "csv"
        filename = f"data/impacts/worldriskindex22_{id}.{fformat}"
        if not os.path.exists(filename):
            with open(filename, "wb") as file:
                    file.write(f.content)
            print("World Risk Index ", id, " downloaded")

def load_riskidx():
    """
    Load data in a DataFrame and keep only columns with the annual percentage of people of that
    country exposed to severe tsunamis, severe coastal floods and sea level rise.

    Returns:
    csv: pandas.DataFrame
        DataFrame containing countries names, ISO codes and indexes of exposition
    """
    csv = pd.read_csv("data/impacts/worldriskindex22_data.csv")
    # cols = ["Country", "ISO3", "EI_02d_Base", "EI_03d_Base", "EI_07b_Base"]
    cols = ["Country", "ISO3", "EI_02d_Norm", "EI_03d_Norm", "EI_07b_Norm"]
    csv = csv.loc[:, cols]
    cols = ["Country", "ISO3", "Tsunamis", "Coastal_floods", "Sea_level_rise"]
    csv.columns = cols
    return(csv)

def plot_riskidx(var, csv = None, addlayer = True, coord = None, m = None):
    """
    Plots variables loaded from load_riskidx
    World countries boundaries loaded from: https://public.opendatasoft.com/explore/dataset/world-administrative-boundaries/export/

    var: str or list of str
        Tags to be plotted selected from the World Risk Index possible indexes
    csv: pandas.DataFrame, optional
        load_riskidx output. Deafult is None, meaning csv will be obtained by calling load_riskidx()
    addlayer: bool, optional
        If True, the function will only add a layer to the ipyleaflet.Map provided
        through m. If False, it will create a new ipyleaflet.Map. Default is True
    coord: tuple of float, optional
        Tuple containing latitude and longitude where to center the map. Only needed if 
        addlayer = False. Default is None
    m: ipyleaflet.Map, optional
        Only if addlayer is True, m needs to be provided as an ipyleaflet.Map object. Default is None
    
    Returns:
    m: ipyleaflet.Map
    """
    #Load risk data
    if csv is None: csv = load_riskidx()
    #Load world countries boundaries
    f = requests.get("https://public.opendatasoft.com/api/explore/v2.1/catalog/datasets/world-administrative-boundaries/exports/geojson?lang=en&timezone=Europe%2FBerlin")
    geo_json_data = json.loads(f.content)
    for d in geo_json_data["features"]:
        d["iso"] = d["properties"]["iso3"]
    #Tool function
    def createandadd(geo_json_data, csv, v, m):
        #Layer selection dictionary
        namedict = {
            "Tsunamis": "Tsunamis Exposition Index - 2022",
            "Coastal_floods": "Coastal floods Exposition Index - 2022",
            "Sea_level_rise": "Sea level rise Exposition Index - 2022"
        }
        #Palettes to be used
        palettedict = {
            "Tsunamis": [(254, 254, 203),(251, 235, 153),(244, 204, 104),(235, 167, 84),(228, 134, 80),(209, 98, 76),(164, 70, 66),(114, 55, 46),(66, 40, 24),(25, 25, 0)],
            "Coastal_floods": [(254, 254, 203),(251, 235, 153),(244, 204, 104),(235, 167, 84),(228, 134, 80),(209, 98, 76),(164, 70, 66),(114, 55, 46),(66, 40, 24),(25, 25, 0)],
            "Sea_level_rise": [(254, 254, 203),(251, 235, 153),(244, 204, 104),(235, 167, 84),(228, 134, 80),(209, 98, 76),(164, 70, 66),(114, 55, 46),(66, 40, 24),(25, 25, 0)]
        }
        palettehex = [rgb_to_hex(x) for x in palettedict[v]]
        mapping = dict(zip(csv["ISO3"].str.strip(), csv[v]))
        for d in geo_json_data["features"]:
            if d["iso"] not in mapping:
                mapping[d["iso"]] = 0
        layer = Choropleth(
                geo_data = geo_json_data,
                choro_data = mapping,
                name = namedict[v],
                style = {'fillOpacity': 0.75, "color":"black"},
                key_on = "iso",
                colormap = bc.StepColormap(palettehex)
                )
        m.add_layer(layer)
        return(m)
    #Plot
    if addlayer:
        if type(var) is list:
            for v in var:
                m = createandadd(geo_json_data, csv, v, m)        
        else:
            m = createandadd(geo_json_data, csv, var, m)
        palette = [(254, 254, 203),(251, 235, 153),(244, 204, 104),(235, 167, 84),(228, 134, 80),(209, 98, 76),(164, 70, 66),(114, 55, 46),(66, 40, 24),(25, 25, 0)]
        palettehex = [rgb_to_hex(x) for x in palette]
        cmap_control = ColormapControl(
                                        caption = "Exposition Indexes",
                                        colormap = bc.StepColormap(palettehex),
                                        value_min = 0,
                                        value_max = 100,
                                        position = 'topright',
                                        transparent_bg = True
                                        )
        m.add(cmap_control)
    else:
        m = Map(center = coord, zoom = 3)
        if type(var) is list:
            for v in var:
                m = createandadd(geo_json_data, csv, v, m)  
        else:
            m = createandadd(geo_json_data, csv, var, m)
        m.add_control(LayersControl())
        m.layout.height="700px"
    return(m)

#%% Resample raster function

def resample_raster(path, fact = 0.5, rio = True, nodata = -3.40282e+38):
    """
    Resample a raster from its file path applying a factor
    ## DOESN'T WORK PROPERLY ##
    Issues:
        - This resampling method will create artifacts in areas with nodata
        - Does not work with Resampling.sum, which was the algorithm needed to 
        - The output file weigth is way higher than the original file 
    Sources used to write the code:
        1. https://rasterio.readthedocs.io/en/stable/topics/resampling.html
        2. https://pygis.io/docs/e_raster_resample.html
    
    path: str
        Path to raster file
    fact: float
        Upscaling or downscaling factor. Default: 0.5 (i.e. halving the raster spatial resolution)
    rio: bool, optional
        If True, a method using rasterio is used. Not yet implemented a method using other 
        packages (e.g. xarray)
    nodata: float, optional
        Value associated to nodata
    
    Returns:
    outpath: str
        Path to the resampled .tiff file
    """
    outpath = path.split(".")[0] + "_resampled." + path.split(".")[1]
    if rio:
        if not os.path.exists(outpath):
            with rasterio.open(path, "r") as dataset:
                # resample data to target shape
                data = dataset.read(
                    out_shape=(
                        dataset.count,
                        int(dataset.height * fact),
                        int(dataset.width * fact)
                    ),
                    resampling = Resampling.bilinear,
                    masked = True
                )
                # scale image transform
                transform = dataset.transform * dataset.transform.scale(
                    (dataset.width / data.shape[-1]),
                    (dataset.height / data.shape[-2])
                )

                # Write outputs
                # set properties for output
                dst_kwargs = dataset.meta.copy()
                dst_kwargs.update(
                    {
                        "crs": dataset.crs,
                        "transform": transform,
                        "width": data.shape[-1],
                        "height": data.shape[-2],
                        "nodata": 0,  
                    }
                )
                with rasterio.open(outpath, "w", **dst_kwargs) as dst:
                    # iterate through bands
                    for i in range(data.shape[0]):
                            dst.write(data[i].astype(rasterio.float64), i+1)
        return(outpath)
    else:
        # method with xarray if needed
        # https://docs.xarray.dev/en/stable/generated/xarray.DataArray.coarsen.html
        pass

#%% Joint plot

def impacts_plot(rp_coh, rp_cyh, coord):
    """
    Displays Section 3 interactive plot by creating an ipyleaflet.Map object 
    and adding a layer for each impact variable

    rp_coh: str
        Return period associated to Coastal hazard. Can be one of the following: 5yr, 10yr, 50yr, 100yr, 250yr, 500yr, 1000yr.
        Provided by a widget in TropiDash_backcone
    rp_cyh: str
        Return period associated to Cyclone hazard. Can be one of the following: 50yr, 100yr, 250yr, 500yr, 1000yr.
        Provided by a widget in TropiDash_backcone
    coord: tuple of float
        Coordinates where to center the ipyleaflet.Map

    Returns:
    None
    """
    #Download
    dwnl_coastalhaz(rp_coh)
    dwnl_cyclonehaz(rp_cyh)

    #Load
    coh = load_coastalhaz(rp_coh, open = True)
    cyh = load_cyclonehaz(rp_cyh, open = True)

    #Create the plot
    m = Map(basemap = basemaps.Esri.WorldTopoMap, center = coord, zoom = 3)
    m = plot_poplayer(m = m)
    m = plot_coastalhaz(coh, rp_coh, m = m)
    m = plot_cyclonehaz(cyh, rp_cyh, m = m)
    m.add_control(LayersControl())
    m.add_control(FullScreenControl())
    m.layout.height = "700px"

    #Show the plot
    display(m)