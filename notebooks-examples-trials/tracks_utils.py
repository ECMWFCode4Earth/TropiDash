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
import geopandas as gpd
import pandas as pd
import numpy as np
import xarray as xr
import haversine as hs
import requests

from ipywidgets import interact
import ipyleaflet
import ipywidgets as widgets

from Magics import macro as magics
from IPython.display import display
from Magics.macro import *

# Function to download ensemble forecast of cyclone tracks data
def download_tracks_forecast(start_date):
    '''
        start_date is the value of the starting date forecast widget, dtype: datetime
        the function returns the actual starting date of the forecast
    '''
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

# Function to import forecast storms file and load it in a dataframe
def create_storms_df():
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
    # Load cyclone dataframe with radius of winds excedding the 35 kts threshold
    df3 = pdbufr.read_bufr('data/tc_test_track_data.bufr',
        columns=("stormIdentifier", "longStormName", "ensembleMemberNumber", "latitude", "longitude",
                 "effectiveRadiusWithRespectToWindSpeedsAboveThreshold"))
    # Add the Wind speed at 10m column to the storms dataframe 
    df_storms["windSpeedAt10M"] = df1.windSpeedAt10M
    df_storms["timePeriod"] = df2.timePeriod
    df_storms["effectiveRadiusWithRespectToWindSpeedsAboveThreshold"] = df3.effectiveRadiusWithRespectToWindSpeedsAboveThreshold
    # Storms with number higher than 10 are not real storms (according to what Fernando said)
    drop_condition = df_storms.stormIdentifier < '11'
    df_storms = df_storms[drop_condition]
    return df_storms

# Function to compute the average coordinate of the ensemble forecasted track at a certain time step
def meanposit(knpf, rlatpf, rlonpf):
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

# Function that returns the list of coordinates and the average radius of wind speed wxcedding 35 kts for the mean forecast track
def mean_forecast_track(df_storm):
    
    # Create 2 dataframe with latitude and logitude coordinates for each ensemble as columns
    members = df_storm.ensembleMemberNumber.unique()
    df_lat_tracks = pd.DataFrame()
    df_lon_tracks = pd.DataFrame()
    df_radius_tracks = pd.DataFrame()
    for member in members:
        df_track = df_storm[df_storm.ensembleMemberNumber == member]
        df_track.reset_index(drop=True, inplace=True)
        df_lat_tracks[f'latitude{member}'] = df_track.latitude
        df_lon_tracks[f'longitude{member}'] = df_track.longitude
        df_radius_tracks[f'radius{member}'] = df_track.effectiveRadiusWithRespectToWindSpeedsAboveThreshold
    # Add date information for the average track
    df_track.reset_index(drop=True, inplace=True)
    start = datetime(df_track.iloc[0]['year'], df_track.iloc[0]['month'], df_track.iloc[0]['day'], df_track.iloc[0]['hour'])
    dates = start + timedelta(hours=6) * (df_track.index+1)
    
    # Cycle through the rows of df_lat_track and df_lon_tracks to compute the average track lat,lon
    mean_track_coord = []
    timesteps = []
    radii = []
    for t in range(len(df_lat_tracks)):
        lat = df_lat_tracks.iloc[t].dropna().to_numpy()
        lon = df_lon_tracks.iloc[t].dropna().to_numpy()
        date = dates[t].strftime("%d-%m-%Y %H:%M")
        rad = df_radius_tracks.iloc[t].dropna()
        radius = rad[rad > 0].mean()
        if len(lat) > 0:
            mean_lat_lon = meanposit(len(lat), lat, lon)
            mean_track_coord.append(mean_lat_lon)
            timesteps.append(date)
            radii.append(radius)
        
    return mean_track_coord, radii, timesteps

# Function to determine the correspondent cyclones between forecast and observed data
def storms_pairing(df_storms_forecast, df_storms_observed):
    for_storms = df_storms_forecast.stormIdentifier.unique()
    obs_storms = df_storms_observed.NAME.squeeze().unique().tolist()

    storms_pair = []

    for f_storm in for_storms:
        dff = df_storms_forecast[df_storms_forecast.stormIdentifier == f_storm]
        loc_f = (dff.iloc[0].latitude, dff.iloc[0].longitude)
        max_dist = 40075 # kms of the equator
        for o_storm in obs_storms:
            dfo = df_storms_observed[df_storms_observed.NAME.squeeze() == o_storm]
            loc_o = (dfo.iloc[-1].LAT.squeeze(), dfo.iloc[-1].LON.squeeze())
            # Compute the distance between two point on earth with the haversine distance (output in km)
            hs_dist = hs.haversine(loc_f, loc_o)
            if hs_dist < max_dist:
                pair = f"{f_storm}-{o_storm}"
                max_dist = hs_dist
        if 'pair' in locals():
            storms_pair.append(pair)
    
    return storms_pair
    
## FUNCTIONS TO PLOT CYCLONE DATA IN IPYLEAFLET ## 

# Function to create the list of (lat,lon) points, timesteps and radii for the forecast ensemble tracks
def forecast_tracks_locations(df_storm_forecast):
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
        radius = df_track.effectiveRadiusWithRespectToWindSpeedsAboveThreshold.tolist()
        locs = []
        tmtstps = []
        for i in range(len(latitude)):
            loc = (latitude[i], longitude[i])
            locs.append(loc)
            tmtstps.append(dates[i].strftime("%d-%m-%Y %H:%M"))
        locations.append(locs)
        timesteps.append(tmtstps)
        radii.append(radius)
    return locations, radii, timesteps

# Function to create the list of (lat, lon) points for the observed track
def observed_track_locations(df_storm_observed):
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
def plot_cyclone_tracks_ipyleaflet(cyclone):
    
    # storm data preparation for plotting
    df_storms_forecast = create_storms_df()
    df_storms_observed = pd.read_csv('data/ibtracs.ACTIVE.list.v04r00.csv', header=[0,1])
    code, name = cyclone.split('-')
    df_f = df_storms_forecast[df_storms_forecast.stormIdentifier == code]
    df_f.reset_index(drop=True, inplace=True)
    df_o = df_storms_observed[df_storms_observed.NAME.squeeze() == name]
    df_o.reset_index(drop=True, inplace=True)
    initial_lat_lon = (df_f.latitude.iloc[0], df_f.longitude.iloc[0])
    
    # create lists of locations
    locations_f, radii_f, timesteps_f = forecast_tracks_locations(df_f)
    locations_o, timesteps_o = observed_track_locations(df_o)
    locations_avg, radii_avg, timesteps_avg = mean_forecast_track(df_f)
    
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
            # name='Track %.02d' % i,
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
    circle_avg = []
    for avg in range(len(locations_avg)):
        marker = ipyleaflet.CircleMarker(
            location = locations_avg[avg],
            radius=1,
            color="black",
            popup=widgets.HTML(value=f'<b> {timesteps_avg[avg]} </b>')
        )
        if radii_avg[avg] > 0:
            circle = ipyleaflet.Circle(
                location=locations_avg[avg],
                radius=int(radii_avg[avg]),
                color="rgba(0, 0, 0, 0)",
                fill_color="red",
                popup=widgets.HTML(value=f'<b> {radii_avg[avg]*10**(-3):.2f} km </b>')
            )
            circle_avg.append(circle)
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
    circles_layer_avg = ipyleaflet.LayerGroup(layers=circle_avg, name='Average circle of wind speeds > 35 kts')
    layer_group_avg = ipyleaflet.LayerGroup(layers=[track_avg, markers_layer_avg], name='Average Forecast Track')
    
    tc_track_map.add_layer(layer_group_avg)
    tc_track_map.add_layer(circles_layer_avg)

    # Add forecasted ensemble tracks to the map
    layer_group_f = ipyleaflet.LayerGroup(layers=[tracks_layer_group, markers_layer_group], name='Forecasted Ensemble Tracks')

    tc_track_map.add_layer(layer_group_f)

    # Add layers widget to the map
    layers_control = ipyleaflet.LayersControl()
    tc_track_map.add_control(layers_control);

    return tc_track_map

# Function for interactive plotting of observed and forecasted tracks with magics
def plot_cyclone_tracks_magics(cyclone, lat_boundaries, lon_boundaries):
    toplot = []
    
    # storm data preparation for plotting
    df_storms_forecast = create_storms_df()
    df_storms_observed = pd.read_csv('data/ibtracs.ACTIVE.list.v04r00.csv', header=[0,1])
    code, name = cyclone.split('-') 
    df_f = df_storms_forecast[df_storms_forecast.stormIdentifier == code]
    members = df_f.ensembleMemberNumber.unique()
    df_o = df_storms_observed[df_storms_observed.NAME.squeeze() == name]
    lat_o = df_o.LAT.squeeze().to_numpy(dtype='float')
    lon_o = df_o.LON.squeeze().to_numpy(dtype='float')
    
    # settings of the geographical area
    bottom, up = lat_boundaries
    left, right = lon_boundaries
    area = mmap(
        subpage_map_projection="cylindrical",
        subpage_lower_left_longitude=int(left),
        subpage_lower_left_latitude=int(bottom),
        subpage_upper_right_longitude=int(right),
        subpage_upper_right_latitude=int(up),
    )
    toplot.append(area)
    
    # settings of the coastline
    coast = mcoast(
        map_coastline_land_shade = "on",
        map_coastline_land_shade_colour = "cream",
        map_coastline_sea_shade = "on",
        map_coastline_sea_shade_colour = "#70CEE2",
        # map_cities = "on",
        map_grid_line_style = "dash",
        map_grid_colour = "black",
        map_label = "on",
        map_label_colour = "brown",
        map_coastline_colour = "brown",
    )
    toplot.append(coast)
    
    # settings of the tracks colour
    colours = ["red", "blue", "green", "yellow", "purple", "orange", "cyan", "brown"]
    colour = 0
    
    # Plot the ensemble tracks of the forecast
    for member in members:

        df_track = df_f[df_f.ensembleMemberNumber == member]
        df_track.dropna(subset = ['latitude', 'longitude'], inplace=True)
        lon = df_track.longitude.to_numpy()
        lat = df_track.latitude.to_numpy()

        data = minput(
            input_type = 'geographical',
            input_x_values = lon,
            input_y_values = lat,
        )

        line=msymb(
            symbol_type='marker',
            symbol_marker_index = 28,
            symbol_colour = colours[colour],
            symbol_height = 0.20,
            symbol_text_font_colour = "black",       
            symbol_connect_line ='on'
        )

        colour += 1
        if colour == len(colours):
            colour = 0

        toplot.append(data)
        toplot.append(line)
    
    # Plot the observed track of the cyclone
    data_obs = minput(
        input_type = 'geographical',
        input_x_values = lon_o,
        input_y_values = lat_o,
    )
    toplot.append(data_obs)

    line_obs = msymb(
        symbol_type='marker',
        symbol_marker_index = 28,
        symbol_colour = "black",
        symbol_height = 0.20,
        symbol_text_font_colour = "black",       
        symbol_connect_line ='on'
    )
    toplot.append(line_obs)
    
    title = mtext(
        text_lines= [f"<font colour='navy'> <b> {cyclone} observed and forecasted trajectories </b> </font>",
                    # f"{start.strftime('%d %b %Y %H:%M')} - {end.strftime('%d %b %Y %H:%M')}"
                    ],
        text_justification= 'centre',
        text_font_size= 0.7,
        # text_font_style= 'bold',
        text_mode='title',
    )
    toplot.append(title)
    
    display(plot(toplot))