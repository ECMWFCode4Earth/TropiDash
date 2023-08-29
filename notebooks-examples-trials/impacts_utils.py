
# -*- coding: utf-8 -*-
# isort: off

#Necessary packages
from ipyleaflet import Choropleth, Map, basemap_to_tiles, LayersControl
import ipywidgets as widgets
from IPython.display import display
from localtileserver import get_leaflet_tile_layer, TileClient
import matplotlib
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
import requests
import json

#%% Resample function

def resample_raster(path, fact = 0.5, rio = True, nodata = -3.40282e+38):
    """
    Resample a raster from its file path applying a factor
    Issues:
        - This resampling method will create artifacts in areas with nodata
        - Does not work with Resampling.sum, which was the algorithm needed to 
    
    path: str
        Path to raster file
    fact: float
        Upscaling or downscaling factor. Default: 0.5 (i.e. halving the raster spatial resolution)
    """
    outpath = path.split(".")[0] + "_resampled." + path.split(".")[1]
    if rio:
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
