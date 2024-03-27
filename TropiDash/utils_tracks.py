# -*- coding: utf-8 -*-
# isort: off

import os
import pdbufr
import math
import ipyleaflet
import rasterio
import pandas as pd
import numpy as np
import xarray as xr
import rioxarray 
import ipywidgets as widgets
import branca.colormap as bc

from eccodes import *
from ecmwf.opendata import Client
from datetime import datetime, timedelta
from itertools import chain, tee
from scipy.optimize import root_scalar
from scipy.spatial import KDTree
from localtileserver import get_leaflet_tile_layer, TileClient

def download_tracks_forecast(start_date):
    """
    Downloads the forecast of the tropical cyclone tracks from ECMWF's open data dataset azure
    and saves them in a bufr file "data/tracks/{date of the forecast}.bufr".
    If today's forecast is not available yet, it downloads yesterday's forecast.
    
    start_date: datetime object
        Day of the forecast of which the data want to be downloaded

    Returns:
    start_date: datetime object
        Day of the forecast of which the data have been downloaded
    """
    # Check if the file already exists
    if os.path.exists(f"data/tracks/{start_date.strftime('%Y%m%d')}.bufr"):
        return start_date
    else:
        # Define the ECMWF server source of data
        client = Client(source="azure")
        try:
            # Download the data from the ECMWF server
            client.retrieve(
                date=int(start_date.strftime("%Y%m%d")),
                time=0,
                stream="enfo",
                type="tf",
                step=240,
                target=f"data/tracks/{start_date.strftime('%Y%m%d')}.bufr",
            );
            return start_date
        # Usually early in the morning the forecast of the current day is not available
        except:
            # Remove possible empty file and download the forecast of the previous day
            os.remove(f"data/tracks/{start_date.strftime('%Y%m%d')}.bufr")
            start_date = start_date - timedelta(days=1)
            # Download the data from the ECMWF server
            client.retrieve(
                date=int(start_date.strftime("%Y%m%d")),
                time=0,
                stream="enfo",
                type="tf",
                step=240,
                target=f"data/tracks/{start_date.strftime('%Y%m%d')}.bufr",
            );
            return start_date

def create_storms_df(start_date):
    """
    Creates a dataframe containing data of the cyclone tracks starting from the forecast file
    "data/tc_track_data_{date of the forecast}.bufr".

    start_date: datetime object
        Day of the forecast of the data selected by user
    
    Returns:
    df_storms: pandas DataFrame
        Dataframe containing the forecast data for the active cyclones.
    """
    # Load cyclone dataframe with Mean sea level pressure value
    df_storms = pdbufr.read_bufr(f"data/tracks/{start_date.strftime('%Y%m%d')}.bufr",
        columns=("stormIdentifier", "longStormName", "ensembleMemberNumber", "year", "month", "day", "hour", "latitude", "longitude",
                 "pressureReducedToMeanSeaLevel"))
    
    # Load cyclone dataframe with Wind speed at 10m value
    df1 = pdbufr.read_bufr(f"data/tracks/{start_date.strftime('%Y%m%d')}.bufr",
        columns=("stormIdentifier", "longStormName", "ensembleMemberNumber", "latitude", "longitude",
                 "windSpeedAt10M"))
    
    # Add the Wind speed at 10m column to the storms dataframe 
    df_storms["windSpeedAt10M"] = df1.windSpeedAt10M
    
    # Storms with stormIdentifier higher than 70 are not real storms
    drop_condition = df_storms.stormIdentifier < '70'
    df_storms = df_storms[drop_condition]

    # Drop storm containing only NaN values
    for cyclone in df_storms.stormIdentifier.unique():
        df_storm = df_storms[df_storms.stormIdentifier == cyclone]
        df_temp = df_storm.dropna(subset=['latitude', 'longitude'])
        if df_temp.empty:
            df_storms = df_storms[df_storms.stormIdentifier != cyclone]
    
    df_storms.reset_index(inplace=True, drop=True)

    # Check if the dataframe is empty return the empty dataframe
    if df_storms.empty:
        pass
    else:
        # Build the dataframe with the timeperiod column
        timeperiod = []
        start_date = datetime(df_storms.year[0], df_storms.month[0], df_storms.day[0], df_storms.hour[0])
        for cyclone in df_storms.stormIdentifier.unique():
            df_cyclone = df_storms[df_storms.stormIdentifier == cyclone]
            df_cyclone.reset_index(inplace=True, drop=True)
            members = df_cyclone.ensembleMemberNumber.unique()
            for member in members:
                df_track = df_cyclone[df_cyclone.ensembleMemberNumber == member]
                for i in range(len(df_track)):
                    timeperiod.append(6 * (i+1))
        
        # Add the timePeriod column to the storms dataframe 
        df_storms["timePeriod"] = timeperiod

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
    Computes the average forecast track of a tropical cyclone given its ensemble members together with the
    10-th, 25-th, 50-th, 75-th and 90-th percentiles for the mean sea level pressure at the center of the cyclone
    and the maximum sustained wind speed at 10m.
    
    df_storm: pandas DataFrame
        Dataframe containing the cyclone forecast data
        
    Returns:
    mean_track_coord: list
        List containing (lat,lon) coordinates of the average track
    timesteps: list
        List containing the cyclone's timesteps
    pressures: list
        List containing the cyclone's mean sea level pressure percentiles
    winds: list
        List containing the cyclone's maximum sustained wind speed at 10m percentiles
    """
    # Create empty dataframes to fill with each ensemble members information
    members = df_storm.ensembleMemberNumber.unique()
    df_lat_tracks = pd.DataFrame()
    df_lon_tracks = pd.DataFrame()
    df_prs_tracks = pd.DataFrame()
    df_wds_tracks = pd.DataFrame()

    # Cycle thorugh the ensemble members and save the information in the just created dataframes
    for member in members:
        df_track = df_storm[df_storm.ensembleMemberNumber == member]
        df_track.reset_index(drop=True, inplace=True)
        df_lat_tracks[f'latitude{member}'] = df_track.latitude
        df_lon_tracks[f'longitude{member}'] = df_track.longitude
        df_prs_tracks[f'pressure{member}'] = df_track.pressureReducedToMeanSeaLevel
        df_wds_tracks[f'wind{member}'] = df_track.windSpeedAt10M
    
    # Add date information for the average track
    df_track.reset_index(drop=True, inplace=True)
    start = datetime(df_track.iloc[0]['year'], df_track.iloc[0]['month'], df_track.iloc[0]['day'], df_track.iloc[0]['hour'])
    dates = start + timedelta(hours=6) * (df_track.index+1)
    
    # Cycle through the rows of df_lat_track and df_lon_tracks to compute the average track lat,lon and pressure and wind speed perecentiles
    mean_track_coord = []
    timesteps = []
    pressures = []
    winds = []
    
    for t in range(len(df_lat_tracks)):
        lat = df_lat_tracks.iloc[t].dropna().to_numpy()
        lon = df_lon_tracks.iloc[t].dropna().to_numpy()
        prs = df_prs_tracks.iloc[t].dropna().to_numpy()  * 10**-2 # Pa to hPa
        wds = df_wds_tracks.iloc[t].dropna().to_numpy()
        date = dates[t].strftime("%d-%m-%Y %H:%M")
        if len(lat) > 0:
            mean_lat_lon = meanposit(len(lat), lat, lon)
            pressures.append(np.percentile(prs, [10, 25, 50, 75, 90]))
            winds.append(np.percentile(wds, [10, 25, 50, 75, 90]))
            mean_track_coord.append(mean_lat_lon)
            timesteps.append(date)
        
    return mean_track_coord, timesteps, pressures, winds

def forecast_tracks_locations(df_storm_forecast):
    """
    Returns 4 lists containing all the cyclone's forecast ensemble tracks information for plotting in ipyleaflet Map
    
    df_storm_forecast: pandas.DataFrame
        Dataframe containing the cyclone forecast data
        
    Returns:
    locations: list
        List containing the lists of the ensemble tracks coordinates
    timesteps: list
        List containing the lists of the ensemble tracks timesteps
    pressures: list
        List containing the lists of the ensemble tracks mean sea level pressure at the storms centers
    wind_speeds: list
        List containing the lists of the ensemble tracks maximum sustained wind speed at 10m
    """
    # List all the ensemble members tracks for the cyclone forecast and define the empty lists
    members = df_storm_forecast.ensembleMemberNumber.unique()
    locations = []
    timesteps = []
    pressures = []
    wind_speeds = []

    # Cycle through the members of the forecast
    for member in members:
        df_track = df_storm_forecast[df_storm_forecast.ensembleMemberNumber == member]
        df_track.reset_index(drop=True, inplace=True)

        # Create a date column in dataframe
        start = datetime(df_track.iloc[0]['year'], df_track.iloc[0]['month'], df_track.iloc[0]['day'], df_track.iloc[0]['hour'])
        df_track['date'] = start + timedelta(hours=6) * (df_track.index+1)
        df_track.dropna(subset = ['latitude', 'longitude'], inplace=True)

        # Save the information from a single track in lists
        latitude = df_track.latitude.tolist()
        longitude = df_track.longitude.tolist()
        dates = df_track.date.tolist()
        pressure = df_track.pressureReducedToMeanSeaLevel.tolist()
        wind_speed = df_track.windSpeedAt10M.tolist()
        locs = []
        tmtstps = []
        press = []

        # Cycle through the list of location
        for i in range(len(latitude)):
            loc = (latitude[i], longitude[i])
            locs.append(loc)
            tmtstps.append(dates[i].strftime("%d-%m-%Y %H:%M"))
            press.append(pressure[i] * 10**-2) # Pa to hPa
        
        locations.append(locs)
        timesteps.append(tmtstps)
        pressures.append(press)
        wind_speeds.append(wind_speed)

    return locations, timesteps, pressures, wind_speeds

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
    # Convert the pandas dataframe columns to lists
    latitude = df_storm_observed.LAT.squeeze().tolist()
    longitude = df_storm_observed.LON.squeeze().tolist()
    dates = df_storm_observed.ISO_TIME.squeeze().tolist()
    timesteps = []
    locations = []

    # Cycle through the lists
    for i in range(len(latitude)):
        loc = (latitude[i], longitude[i])
        date = datetime.strptime(dates[i][:-3], ('%Y-%m-%d %H:%M'))
        locations.append(loc)
        timesteps.append(date.strftime('%d-%m-%Y %H:%M'))

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
    df.reset_index(drop=True, inplace=True)
    dates = datetime(df.year[0], df.month[0], df.day[0], df.hour[0]) + timedelta(hours=1) * df.timePeriod
    df["date"] = dates.dt.strftime("%Y%m%d")
    df["step"] = dates.dt.strftime("%H").astype("int") * 100
    df.drop(columns=['year','month','day','hour'])
    df = df.reindex(['lat','lon','number','date','step','wind','msl'], axis=1)
    return df

def strike_probability_map(df_storm_forecast):
    """
    Computes the strike probability map of a cyclone and saves it to raster file
    "data/tracks/pts_raster.tif"
    
    df_storm_forecast: pandas.DataFrame
        Dataframe containing the cyclone forecast tracks data
        
    Returns:
    strike_map: xarray.DataArray
        Dataarray containing the coordinates and the values of the strike probability map
    tif_path: str
        Filename of the raster file containing the strike probability map
    """
    # Reorganization of the dataframe to adjust it to pts algorithm
    storm_code = df_storm_forecast.stormIdentifier.unique()[0]
    df = storm_df_reorganization(df_storm_forecast.copy())
    df.dropna(subset = ['lat', 'lon'], inplace=True)
    df["id"] = df.number
    forecast_date = min(df.date)
    
    # Set of general parameters for the pts algorithm
    distance = 200.0e3
    overlap = 0.7
    dist_circle = distance_from_overlap(distance, overlap)
    basetime = datetime.strptime(min(df.date), '%Y%m%d')
    filter_wind = 0.0
    
    # K-d tree building
    lat_max = df.lat.max() + 10
    lat_min = df.lat.min() - 10
    lon_max = df.lon.max() + 10
    lon_min = df.lon.min() - 10
    
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
            track = tracks[(tracks.id == id) & (filter_wind <= tracks.wind)].sort_values("t")

            # special cases
            if track.shape[0] == 1:
                p = ll_to_ecef(track.lat.iat[0], track.lon.iat[0])
                pts.update(tree.query_ball_point(p, r=distance))
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
    
    tif_path = f"data/tracks/pts_raster_{storm_code}_{forecast_date}.tif"
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
    layer_group_avg: ipyleaflet.LayerGroup
        Layer group containing the average forecast track information
    stp_map: ipyleaflet.TileLayer
        Tile layer containing the strike probability map
    cmap_control: ipyleaflet.ColormapControl
        Colormap widget for the strike probability map
    """
    # storm data preparation for plotting
    df_f = df_storm_forecast[df_storm_forecast.ensembleMemberNumber.isin(ens_members)]
    df_f.reset_index(drop=True, inplace=True)
    initial_lat_lon = (df_storm_forecast.latitude.iloc[0], df_storm_forecast.longitude.iloc[0])
    
    # compute the strike probability map and open raster file
    df_storm = df_storm_forecast.copy()
    strike_map_xr, tif_path = strike_probability_map(df_storm)
    client = TileClient(tif_path)
    
    # create lists of locations
    locations_f, timesteps_f, pressures_f, winds_f = forecast_tracks_locations(df_f)
    locations_o, timesteps_o = observed_track_locations(df_storm_observed)
    locations_avg, timesteps_avg, pressures_avg, winds_avg = mean_forecast_track(df_storm_forecast)
    
    # Create the basemap for plotting
    tc_track_map = ipyleaflet.Map(
        center=initial_lat_lon,
        basemap=ipyleaflet.basemaps.Esri.WorldTopoMap,
        zoom = 3.0,
        # scroll_wheel_zoom=True,
    )

    # Define forecasted tracks polyline element for the map
    colours = ["red", "blue", "green", "yellow", "purple", "orange", "cyan", "brown"]
    tracks_layer = []
    colour = 0
    i = 0
    # Cycle on the ensembles of the forecast track
    for locs in locations_f:

        tmtstps = timesteps_f[i]
        press = pressures_f[i]
        wind = winds_f[i]
        
        # Define the ensemble track polyline for the map
        track = ipyleaflet.Polyline(
            locations=locs,
            color=colours[colour],
            fill=False,
            weight=2,
        )
        if len(locations_f) <= 5:
            # Define the markers element for each position of the cyclone ensemble forecast
            markers = []
            for j in range(len(locs)):
                marker = ipyleaflet.CircleMarker(
                    location=locs[j],
                    radius=1,
                    color=colours[colour],
                    popup=widgets.HTML(value=f"<center><b>VT: {tmtstps[j]}</b> </center>"
                                       f"<center> ({locs[j][0]:.2f}, {locs[j][1]:.2f}) </center>"
                                       f"<center> Pressure: {press[j]:.2f} hPa</center>" 
                                       f"<center>Wind speed: {wind[j]:.2f} m/s</center>")
                )
                markers.append(marker)
            markers_layer = ipyleaflet.LayerGroup(layers=markers)
            track_layer = ipyleaflet.LayerGroup(layers=[track, markers_layer])
            tracks_layer.append(track_layer)
        else:
            tracks_layer.append(track)
            
        colour += 1
        if colour == len(colours):
            colour = 0
        
        i += 1
        
    # Define average forecast polyline for the map
    track_avg = ipyleaflet.Polyline(
            locations=locations_avg,
            color="black",
            fill=False,
            weight=3,
            # name='Average Forecast Track',
        )
    
    # Define the markers element for each position of the average track
    marker_avg = []
    for avg in range(len(locations_avg)):
        marker = ipyleaflet.CircleMarker(
            location = locations_avg[avg],
            radius=1,
            color="black",
            popup=widgets.HTML(value=f"<center><b>VT: {timesteps_avg[avg]} </b> </center>"
                               f"<center> ({locations_avg[avg][0]:.2f}, {locations_avg[avg][1]:.2f}) </center>"
                               f"Percentiles: Pressure || Wind speed <br>"
                               f"10<sup>th</sup>: {pressures_avg[avg][0]:.1f} hPa || {winds_avg[avg][0]:.2f} m/s <br>"
                               f"25<sup>th</sup>: {pressures_avg[avg][1]:.1f} hPa || {winds_avg[avg][1]:.2f} m/s <br>"
                               f"50<sup>th</sup>: {pressures_avg[avg][2]:.1f} hPa || {winds_avg[avg][2]:.2f} m/s <br>"
                               f"75<sup>th</sup>: {pressures_avg[avg][3]:.1f} hPa || {winds_avg[avg][3]:.2f} m/s <br>"
                               f"90<sup>th</sup>: {pressures_avg[avg][4]:.1f} hPa || {winds_avg[avg][4]:.2f} m/s"
                               )
        )
        marker_avg.append(marker)
    
    # Define observed tracks polyline element for the map
    track_o = ipyleaflet.Polyline(
            locations=locations_o,
            color= "#ff00ff",
            fill=False,
            weight=2,
        )

    # Define the markers element for each position of the observed track
    marker_o = []
    for o in range(len(locations_o)):
        marker = ipyleaflet.CircleMarker(
            location = locations_o[o],
            radius=1,
            color="#ff00ff",
            popup=widgets.HTML(value=f"<b>VT: {timesteps_o[o]} </b>"
                               f"<center> ({locations_o[o][0]:.2f}, {locations_o[o][1]:.2f}) </center>"
                               )
            )
        marker_o.append(marker)
        
    # Define raster layer for strike probability map and its colormap widget
    # palette = ["#c4ff70", "#6ae24c", "#2a9134", "#137547", "#046335", "#2397d1", "#557ff3", "#143cdc", "#3910b4", "#1e0063"]
    palette = ["#8df52c", "#6ae24c", "#61bb30", "#508b15", "#057941", "#2397d1", "#557ff3", "#143cdc", "#3910b4", "#1e0063"]
    stp_map = get_leaflet_tile_layer(client, name = "Strike Probability Map", opacity = 0.8, palette = palette, nodata=0.0, max_zoom = 30)
    
    with rasterio.open(tif_path) as r:
        minv = "%.2f" % round(r.read(1).ravel().min(), 1)
        maxv = "%.2f" % round(r.read(1).ravel().max(), 1)
    
    cmap_control = ipyleaflet.ColormapControl(
                                caption = "Strike probability",
                                colormap = bc.StepColormap(palette),
                                value_min = float(minv),
                                value_max = float(maxv),
                                position = 'topright',
                                transparent_bg = True
                                )
    
    # Add observed track to the map
    markers_layer_o = ipyleaflet.LayerGroup(layers=marker_o)
    layer_group_o = ipyleaflet.LayerGroup(layers=[track_o, markers_layer_o], name='Observed Track')

    tc_track_map.add_layer(layer_group_o)

    # Add average forecast track to the map
    markers_layer_avg = ipyleaflet.LayerGroup(layers=marker_avg)
    layer_group_avg = ipyleaflet.LayerGroup(layers=[track_avg, markers_layer_avg], name='Average Forecast Track')
    
    tc_track_map.add_layer(layer_group_avg)

    # Add forecasted ensemble tracks to the map
    tracks_layer_group = ipyleaflet.LayerGroup(layers=tracks_layer, name='Forecasted Ensemble Tracks')
    tc_track_map.add_layer(tracks_layer_group)
    
    # Add strike probability map layer to the map and colormap
    tc_track_map.add_layer(stp_map)
    tc_track_map.add_control(cmap_control)

    # Add layers widget to the map
    layers_control = ipyleaflet.LayersControl()
    tc_track_map.add_control(layers_control);
    
    # Add posibility to open the map full screen
    tc_track_map.add_control(ipyleaflet.FullScreenControl())

    return(tc_track_map, [layer_group_avg, stp_map, cmap_control])