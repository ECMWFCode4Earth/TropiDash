# -*- coding: utf-8 -*-
# isort: off

import branca.colormap as bc
from ipyleaflet import Map, ColormapControl, LayersControl, basemaps, FullScreenControl
from ipyleaflet.velocity import Velocity
from IPython.display import display
from localtileserver import get_leaflet_tile_layer, TileClient
import rasterio
import rioxarray as rxr
import xarray as xr

from utils_atm import sel_forecast, get_palette
from utils_impacts import *

def plot_section4(vardict, step, coord, stepsdict, rp_coh, rp_cyh, cyclonelayers):
    """
    Plots everything together

    vardict: dict
        Dictionary containing the rasters uploaded trhough load_atmdata.
        Returned by load_atmdata
    step: str
        String returned from the widget used to select the forecast period
    coord: tuple
        Coordinates of the map central point. Provide them as lat, lon
    stepsdict: dict
        Dictionary containing the steps codes needed for each variable. The standard step format is under "base".
    rp_coh:
        Return period of the coastal hazard layer to be added
    rp_cyh:
        Return period of the cyclone hazard layer to be added
    cyclonelayers: 
        Layers of cyclone related variables
        
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
    
    # Add impact layers
    #Download
    dwnl_coastalhaz(rp_coh)
    dwnl_cyclonehaz(rp_cyh)
    #Load
    coh = load_coastalhaz(rp_coh, open = True)
    cyh = load_cyclonehaz(rp_cyh, open = True)
    #Add
    m = plot_coastalhaz(coh, rp_coh, m = m)
    m = plot_cyclonehaz(cyh, rp_cyh, m = m)
    m = plot_poplayer(m = m)

    # Add cyclone layers
    for layer in cyclonelayers:
        m.add_layer(layer)

    m.add_control(LayersControl())
    m.add_control(FullScreenControl())
    m.layout.height = "700px"
    display(m)