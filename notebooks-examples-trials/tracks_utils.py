# -*- coding: utf-8 -*-
# isort: off

import pdbufr
import sys
import traceback
 
from math import isnan
import math

from eccodes import *
from ecmwf.opendata import Client
from datetime import datetime, timedelta

import os
import birdy
import pandas as pd
import numpy as np
# import xarray as xr
import haversine as hs
import requests

from ipywidgets import interact
import ipyleaflet
import ipywidgets as widgets

from Magics import macro as magics
from IPython.display import display
from Magics.macro import *

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

# Function to plot the interactive map with ipyleaflet given a cyclone
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