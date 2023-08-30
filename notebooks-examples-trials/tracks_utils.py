# -*- coding: utf-8 -*-
# isort: off

import os
import birdy
import requests
import pdbufr
import sys
import traceback
import re
import math
import ipyleaflet
import rasterio
import rioxarray
import pandas as pd
import numpy as np
import xarray as xr
import haversine as hs
import ipywidgets as widgets

from math import isnan
from eccodes import *
from ecmwf.opendata import Client
from datetime import datetime, timedelta
from itertools import chain, tee
from ipywidgets import interact
from IPython.display import display
from scipy.optimize import root_scalar
from scipy.spatial import KDTree

def download_tracks_forecast(start_date):
    """
    Downloads the forecast of the tropical cyclone tracks from ECMWF's open data dataset azure
    and saves them in a bufr file "data/tc_test_track_data.bufr".
    If today's forecast is not available yet, it downloads yesterday's forecast.
    
    start_date: datetime object
        Day of the forecast of which the data want to be downloaded
    """
    client = Client(source="azure")
    try:
        client.retrieve(
            date=int(start_date.strftime("%Y%m%d")),
            time=0,
            stream="enfo",
            type="tf",
            step=240,
            target="data/tc_test_track_data.bufr",
        );
    # Usually early in the morning the forecast of the current day is not available
    except:
        start_date = start_date - timedelta(days=1)
        client.retrieve(
            date=int(start_date.strftime("%Y%m%d")),
            time=0,
            stream="enfo",
            type="tf",
            step=240,
            target="data/tc_test_track_data.bufr",
        );
        print(f'Today\'s forecast not available, downloaded yesterday\'s forecast: {start_date.strftime("%d %b %Y")}')

def create_storms_df():
    """
    Creates a dataframe containing data of the cyclone tracks starting from the forecast file
    "data/tc_test_track_data.bufr". 
    
    Returns:
    df_storms: pandas DataFrame
        Dataframe containing the forecast data for the active cyclones.
    """
    # Load cyclone dataframe with Mean sea level pressure value
    df_storms = pdbufr.read_bufr('data/tc_test_track_data.bufr',
        columns=("stormIdentifier", "longStormName", "ensembleMemberNumber", "year", "month", "day", "hour", "latitude", "longitude",
                 "pressureReducedToMeanSeaLevel"))
    # Load cyclone dataframe with Wind speed at 10m value
    df1 = pdbufr.read_bufr('data/tc_test_track_data.bufr',
        columns=("stormIdentifier", "longStormName", "ensembleMemberNumber", "latitude", "longitude",
                 "windSpeedAt10M"))
    # Load cyclone dataframe with timeperiod column
    df2 = pdbufr.read_bufr('data/tc_test_track_data.bufr',
        columns=("stormIdentifier", "longStormName", "ensembleMemberNumber", "latitude", "longitude",
                 "timePeriod"))
    # Add the Wind speed at 10m column to the storms dataframe 
    df_storms["windSpeedAt10M"] = df1.windSpeedAt10M
    df_storms["timePeriod"] = df2.timePeriod
    # Storms with number higher than 10 are not real storms (according to what Fernando said)
    drop_condition = df_storms.stormIdentifier < '11'
    df_storms = df_storms[drop_condition]
    return df_storms

def meanposit(knpf, rlatpf, rlonpf):
    """
    Computes the average coordinates from the ensemble tracks at a specific timestep.
    
    knpf: int
        Number of ensemble tracks 
    rlatpf: numpy array
        Array containing the latitudes of the ensemble tracks
    rlonpf: numpy array
        Array containing the longitudes of the ensemble tracks
    
    Returns:
    rlatmean: float
        Mean latitude at the timestep
    rlonmean: float
        Mean longitude at the timestep
    """
    rpi = math.acos(0.0)
    rnomin = 0.0
    rdenom = 0.0
    rphisum = 0.0
    
    for k in range(knpf):
        rlat = rlatpf[k] * rpi / 180.0
        rlon = rlonpf[k] * rpi / 180.0
        rcosphi = math.cos(rlat)
        rnomin += rcosphi * rlon
        rdenom += rcosphi
        rphisum += rlat
    
    rlabda = rnomin / rdenom
    rlonmean = rlabda * 180.0 / rpi
    
    rnomin = 0.0
    repsilon = 0.0
    
    for k in range(knpf):
        rlat = rlatpf[k] * rpi / 180.0
        rlon = rlonpf[k] * rpi / 180.0
        rcosphi = math.cos(rlat)
        rnomin += rcosphi * (rlabda - rlon) ** 2
    
    repsilon = rnomin / (2.0 * knpf)
    rphimean = rphisum / float(knpf)
    rphi = rphimean + repsilon * math.sin(rphimean)
    rlatmean = rphi * 180.0 / rpi
    
    return rlatmean, rlonmean

def mean_forecast_track(df_storm):
    """
    Computes the average forecast track of a tropical cyclone given its ensemble members.
    
    df_storm: pandas DataFrame
        Dataframe containing the cyclone forecast data
        
    Returns:
    mean_track_coord: list
        List containing (lat,lon) coordinates of the average track
    timesteps: list
        List containing the cyclone's timesteps
    """
    # Create 2 dataframe with latitude and logitude coordinates for each ensemble as columns
    members = df_storm.ensembleMemberNumber.unique()
    df_lat_tracks = pd.DataFrame()
    df_lon_tracks = pd.DataFrame()
    for member in members:
        df_track = df_storm[df_storm.ensembleMemberNumber == member]
        df_track.reset_index(drop=True, inplace=True)
        df_lat_tracks[f'latitude{member}'] = df_track.latitude
        df_lon_tracks[f'longitude{member}'] = df_track.longitude
    # Add date information for the average track
    df_track.reset_index(drop=True, inplace=True)
    start = datetime(df_track.iloc[0]['year'], df_track.iloc[0]['month'], df_track.iloc[0]['day'], df_track.iloc[0]['hour'])
    dates = start + timedelta(hours=6) * (df_track.index+1)
    
    # Cycle through the rows of df_lat_track and df_lon_tracks to compute the average track lat,lon
    mean_track_coord = []
    timesteps = []
    # radii = []
    for t in range(len(df_lat_tracks)):
        lat = df_lat_tracks.iloc[t].dropna().to_numpy()
        lon = df_lon_tracks.iloc[t].dropna().to_numpy()
        date = dates[t].strftime("%d-%m-%Y %H:%M")
        if len(lat) > 0:
            mean_lat_lon = meanposit(len(lat), lat, lon)
            mean_track_coord.append(mean_lat_lon)
            timesteps.append(date)
        
    return mean_track_coord, timesteps

def forecast_tracks_locations(df_storm_forecast):
    """
    Returns the list containing all the cyclone's forecast ensemble tracks for plotting in ipyleaflet Map
    
    df_storm_forecast: pandas.DataFrame
        Dataframe containing the cyclone forecast data
        
    Returns:
    locations: list
        List containing the lists of the ensemble tracks coordinates
    timesteps: list
        List containing the lists of the ensemble tracks timesteps
    """
    members = df_storm_forecast.ensembleMemberNumber.unique()
    locations = []
    timesteps = []
    radii = []
    for member in members:
        df_track = df_storm_forecast[df_storm_forecast.ensembleMemberNumber == member]
        df_track.reset_index(drop=True, inplace=True)
        # create timestep column in dataframe
        start = datetime(df_track.iloc[0]['year'], df_track.iloc[0]['month'], df_track.iloc[0]['day'], df_track.iloc[0]['hour'])
        df_track['date'] = start + timedelta(hours=6) * (df_track.index+1)
        df_track.dropna(subset = ['latitude', 'longitude'], inplace=True)
        latitude = df_track.latitude.tolist()
        longitude = df_track.longitude.tolist()
        dates = df_track.date.tolist()
        locs = []
        tmtstps = []
        for i in range(len(latitude)):
            loc = (latitude[i], longitude[i])
            locs.append(loc)
            tmtstps.append(dates[i].strftime("%d-%m-%Y %H:%M"))
        locations.append(locs)
        timesteps.append(tmtstps)
    return locations, timesteps

def observed_track_locations(df_storm_observed):
    """
    Returns the list containing all the coordinates of the observed track for plotting in ipyleaflet Map
    
    df_storm_observed: pandas.DataFrame
        Dataframe containing the cyclone observed data
        
    Returns:
    locations: list
        List containing the lists of the observed track coordinates
    timesteps: list
        List containing the lists of the observed track timesteps
    """
    latitude = df_storm_observed.LAT.squeeze().tolist()
    longitude = df_storm_observed.LON.squeeze().tolist()
    dates = df_storm_observed.ISO_TIME.squeeze().tolist()
    timesteps = []
    locations = []
    for i in range(len(latitude)):
        loc = (latitude[i], longitude[i])
        date = dates[i][:-3]
        locations.append(loc)
        timesteps.append(date)
    return locations, timesteps

# Series of function needed to compute the strike probability map
def previous_and_current(some_iterable):
    prevs, currs = tee(some_iterable, 2)
    prevs = chain([None], prevs)
    return zip(prevs, currs)

def ll_to_ecef(lat, lon, height=0.0, radius=6371229.0):
    lonr = np.radians(lon)
    latr = np.radians(lat)

    x = (radius + height) * np.cos(latr) * np.cos(lonr)
    y = (radius + height) * np.cos(latr) * np.sin(lonr)
    z = (radius + height) * np.sin(latr)
    return x, y, z

def distance_from_overlap(radius, overlap):
    assert 0.0 < radius
    assert 0.0 <= overlap < 1.0
    if overlap <= 0.0:
        return np.inf

    def overlap_unit_circles(d_over_r):
        assert 0.0 <= d_over_r
        if 2.0 <= d_over_r:
            return 0.0

        hd = d_over_r / 2.0
        ha_inter = np.arccos(hd) - hd * np.sqrt(1.0 - hd * hd)
        ha_union = np.pi - ha_inter
        return ha_inter / ha_union

    d = root_scalar(lambda d: overlap_unit_circles(d) - overlap, bracket=[0, 2], x0=1)
    return radius * d.root

def delta_hours(a: datetime, b: datetime) -> int:
    delta = a - b
    return delta // timedelta(hours=1)

def storm_df_reorganization(df):
    df.rename(columns={"ensembleMemberNumber":"number", "latitude":"lat", "longitude":"lon", "pressureReducedToMeanSeaLevel":"msl", "windSpeedAt10M":"wind"}, inplace=True)
    df.drop(columns=["stormIdentifier"])
    df['msl'] = df['msl'] / 100
    dates = datetime(df.year[0], df.month[0], df.day[0], df.hour[0]) + timedelta(hours=1) * df.timePeriod
    df["date"] = dates.dt.strftime("%Y%m%d")
    df["step"] = dates.dt.strftime("%H").astype("int") * 100
    df.drop(columns=['year','month','day','hour'])
    df = df.reindex(['lat','lon','number','date','step','wind','msl'], axis=1)
    return df

def strike_probability_map(df_storm_forecast):
    """
    Computes the strike probability map of a cyclone and saves it to raster file
    "data/pts_raster.tif"
    
    df_storm_forecast: pandas.DataFrame
        Dataframe containing the cyclone forecast tracks data
        
    Returns:
    strike_map: xarray.DataArray
        Dataarray containing the coordinates and the values of the strike probability map
    tif_path: str
        Filename of the raster file containing the strike probability map
    """
    # Reorganization of the dataframe to adjust it to pts algorithm
    df = storm_df_reorganization(df_storm_forecast)
    df.dropna(subset = ['lat', 'lon'], inplace=True)
    df["id"] = df.number
    
    # Set of general parameters for the pts algorithm
    distance = 300.0e3
    overlap = 0.7
    dist_circle = distance_from_overlap(distance, overlap)
    basetime = datetime.strptime(df.date[0], '%Y%m%d')
    
    # K-d tree building
    lat_max = df.lat.max() + 5
    lat_min = df.lat.min() - 5
    lon_max = df.lon.max() + 5
    lon_min = df.lon.min() - 5
    
    lats = np.flip(np.arange(lat_min, lat_max, 0.25))
    lons = np.arange(lon_min, lon_max, 0.25)
    
    N = len(lats) * len(lons)
    Ni = len(lons)
    Nj = len(lats)
    
    P = np.empty([N,3])
    
    i = 0
    for lat in lats:
        for lon in lons:
            P[i, :] = ll_to_ecef(lat, lon)
            i+=1
    
    tree = KDTree(P)
    
    # pre-process (calculate/drop columns)
    if df.empty:
        df[["lat", "lon", "number", "t", "wind", "msl"]] = None
        print("Warning:", df)
    else:
        datestep = [
            datetime.strptime(k, "%Y%m%d%H%M")
            for k in (
                df.date.astype(str) + df.step.apply(lambda s: str(s).zfill(4))
            )
        ]
        if not basetime:
            basetime = min(datestep)
        df["t"] = [delta_hours(ds, basetime) for ds in datestep]

        df.drop(["date", "step"], axis=1, inplace=True)
        
    # probability field (apply filter_number)
    val = np.zeros(N)
    numbers = sorted(set(df.number.tolist()))

    for number in numbers:
        pts = set()

        tracks = df[df.number == number]
        for id in set(tracks.id.tolist()):
            track = tracks[tracks.id == id].sort_values("t")

            # special cases
            if track.shape[0] == 1:
                if args.verbosity >= 1:
                    print(f"number={number} segments=0 len=1")
                p = ll_to_ecef(track.lat.iat[0], track.lon.iat[0])
                pts.update(tree.query_ball_point(p, r=args.distance))
                continue

            ti = np.array([])
            tend = None
            npoints = 1
            for a, b in previous_and_current(track.itertuples()):
                if a is not None:
                    tend = b.t
                    npoints += 1

                    # approximate distance(a, b) with Cartesian distance
                    ax, ay, az = ll_to_ecef(a.lat, a.lon)
                    bx, by, bz = ll_to_ecef(b.lat, b.lon)
                    dist_ab = np.linalg.norm(np.array([bx - ax, by - ay, bz - az]))

                    num = max(1, int(np.ceil(dist_ab / dist_circle)))
                    ti = np.append(
                        ti, np.linspace(a.t, b.t, num=num, endpoint=False)
                    )
            if not tend:
                continue
            ti = np.append(ti, tend)

            assert 0 < npoints == track.shape[0]

            lati = np.interp(ti, track.t, track.lat)
            loni = np.interp(ti, track.t, track.lon)

            # track points
            x, y, z = ll_to_ecef(lati, loni)
            for p in zip(x, y, z):
                pts.update(tree.query_ball_point(p, r=distance))

        for i in pts:
            assert i < N
            val[i] = val[i] + 1.0

    if numbers:
        val = (val / len(numbers)) * 100.0  # %

    if len(numbers) > 1:

        def ranges(i):
            from itertools import groupby

            for _, b in groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
                b = list(b)
                yield b[0][1], b[-1][1]

        mx = max(val)

    assert 0 <= min(val) and max(val) <= 100.0
    
    # Format the algorithm result to a xarray.DataArray
    strike_map = val.reshape((Nj, Ni))
    strike_map_xr = xr.DataArray(strike_map, dims=('latitude', 'longitude'), coords={'latitude': lats, 'longitude': lons})
    
    tif_path = "data/pts_raster.tif"
    strike_map_xr = strike_map_xr.rio.write_crs("epsg:4326")
    strike_map_xr.rio.to_raster(tif_path)
    
    return strike_map_xr, tif_path

def plot_cyclone_tracks_ipyleaflet(ens_members, df_storm_forecast, df_storm_observed):
    """
    Creates the ipyleaflet Map containing the cyclone track data
    
    ens_members: list
        List containing the ensemble members to consider
    df_storm_forecast: pandas.DataFrame
        Dataframe containing the cyclone forecast data
    df_storm_observed: pandas.DataFrame
        Dataframe containing the cyclone observed data
        
    Returns:
    tc_track_map: ipyleaflet.Map
        Map with cyclone track forecast and observation
    """
    # storm data preparation for plotting
    df_f = df_storm_forecast[df_storm_forecast.ensembleMemberNumber.isin(ens_members)]
    df_f.reset_index(drop=True, inplace=True)
    initial_lat_lon = (df_storm_forecast.latitude.iloc[0], df_storm_forecast.longitude.iloc[0])
    
    # create lists of locations
    locations_f, timesteps_f = forecast_tracks_locations(df_f)
    locations_o, timesteps_o = observed_track_locations(df_storm_observed)
    locations_avg, timesteps_avg = mean_forecast_track(df_storm_forecast)
    
    # Create the basemap for plotting
    tc_track_map = ipyleaflet.Map(
        center=initial_lat_lon,
        basemap=ipyleaflet.basemaps.OpenStreetMap.France,
        zoom = 3.0,
        # scroll_wheel_zoom=True,
    )

    # Define forecasted tracks polyline element for the map
    colours = ["red", "blue", "green", "yellow", "purple", "orange", "cyan", "brown"]
    tracks_layer_group = ipyleaflet.LayerGroup()
    markers_layer_group = ipyleaflet.LayerGroup()
    colour = 0
    i = 0
    # Cycle on the ensembles of the forecast track
    for locs in locations_f:

        tmtstps = timesteps_f[i]
        
        # Define the ensemble track polyline for the map
        track = ipyleaflet.Polyline(
            locations=locs,
            color=colours[colour],
            fill=False,
            weight=2,
        )
        # Define the markers element for each position of the cyclone ensemble forecast
        markers = []
        for j in range(len(locs)):
            marker = ipyleaflet.CircleMarker(
                location=locs[j],
                radius=1,
                color=colours[colour],
                popup=widgets.HTML(value=f'<b> {tmtstps[j]} </b>')
            )
            markers.append(marker)

        colour += 1
        if colour == len(colours):
            colour = 0

        marker_layer = ipyleaflet.LayerGroup(layers=markers)

        tracks_layer_group.add_layer(track)
        markers_layer_group.add_layer(marker_layer)
        
        i += 1
        
    # Define average forecast polyline for the map
    track_avg = ipyleaflet.Polyline(
            locations=locations_avg,
            color="black",
            fill=False,
            weight=3,
            # name='Average Forecast Track',
        )
    
    # Define the markers element and the circles element for each position of the average track
    marker_avg = []
    for avg in range(len(locations_avg)):
        marker = ipyleaflet.CircleMarker(
            location = locations_avg[avg],
            radius=1,
            color="black",
            popup=widgets.HTML(value=f'<b> {timesteps_avg[avg]} </b>')
        )
        marker_avg.append(marker)
    
    # Define observed tracks polyline element for the map
    track_o = ipyleaflet.Polyline(
            locations=locations_o,
            color= "black",
            fill=False,
            weight=2,
        )

    # Define the markers element for each position of the observed track
    marker_o = []
    for o in range(len(locations_o)):
        marker = ipyleaflet.CircleMarker(
            location = locations_o[o],
            radius=1,
            color="black",
            popup=widgets.HTML(value=f'<b> {timesteps_o[o]} </b>')
        )
        marker_o.append(marker)
    
    # Add observed track to the map
    markers_layer_o = ipyleaflet.LayerGroup(layers=marker_o)
    layer_group_o = ipyleaflet.LayerGroup(layers=[track_o, markers_layer_o], name='Observed Track')

    tc_track_map.add_layer(layer_group_o)

    # Add average forecast track to the map
    markers_layer_avg = ipyleaflet.LayerGroup(layers=marker_avg)
    layer_group_avg = ipyleaflet.LayerGroup(layers=[track_avg, markers_layer_avg], name='Average Forecast Track')
    
    tc_track_map.add_layer(layer_group_avg)

    # Add forecasted ensemble tracks to the map
    layer_group_f = ipyleaflet.LayerGroup(layers=[tracks_layer_group, markers_layer_group], name='Forecasted Ensemble Tracks')

    tc_track_map.add_layer(layer_group_f)

    # Add layers widget to the map
    layers_control = ipyleaflet.LayersControl()
    tc_track_map.add_control(layers_control);

    return tc_track_map