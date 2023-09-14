import numpy as np
import xarray as xr
import rioxarray

from eccodes import *
from datetime import datetime, timedelta
from itertools import chain, tee
from scipy.optimize import root_scalar
from scipy.spatial import KDTree

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