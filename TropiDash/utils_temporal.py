from bqplot import Lines, Bars, Figure, LinearScale, Axis, DateScale
import ipyleaflet
import ipywidgets as widgets
import pandas as pd
import numpy as np
import xarray as xr

def get_forecast_datasets(file_names):
    '''
    Loads the forecast of the variables in the variables list for the specified forecast date and returns a dataset 
    with the variables data.
    The variables in the list are:
    - Mean sea level pressure
    - Skin temperature
    - Total precipitation
    - Wind components at 10m height
    - Probability of 10m Wind gusts over 25 knots

    file_names: list
        list with the names of the variable files to load
    
    Returns:
    ds_msl: xarray.Dataset
        dataset for variable msl
    ds_skt: xarray.Dataset
        dataset for variable skt
    ds_tp: xarray.Dataset
        dataset for variable tp
    ds_wind : xarray.Dataset
        dataset for variable wind components and wind speed at 10m
    ds_fgg: xarray.Dataset
        dataset for variable 10fgg25
    '''
    # Load the datasets for the variables msl, skt and tp
    filename = file_names[0]
    ds_msl = xr.open_dataset(filename, engine="cfgrib")
    filename = file_names[1]
    ds_skt = xr.open_dataset(filename, engine="cfgrib")
    filename = file_names[2]
    ds_tp = xr.open_dataset(filename, engine="cfgrib")
    filename = file_names[3]
    ds_wind = xr.open_dataset(filename, engine="cfgrib")
    filename = file_names[4]
    ds_fgg = xr.open_dataset(filename, engine="cfgrib")

    # Compute the wind speed from the wind components
    ds_wind['ws10'] = np.sqrt(ds_wind.u10**2 + ds_wind.v10**2)

    # Convert skin temperature to Celsius and mean sea level pressure to hPa
    ds_skt['skt'] = ds_skt.skt - 273.15 # Convert to Celsius
    ds_msl['msl'] = ds_msl.msl / 100 # Convert to hPa

    return ds_msl, ds_skt, ds_tp, ds_wind, ds_fgg

def get_df_pos(ds_msl, ds_skt, ds_tp, ds_wind, ds_fgg, position):
    '''
    Extracts the data for the variables in the datasets for the specified position and returns the data in a 
    different dataframe for each variable.

    ds_msl: xarray.Dataset
        dataset for variable msl
    ds_skt: xarray.Dataset
        dataset for variable skt
    ds_tp: xarray.Dataset
        dataset for variable tp
    ds_wind : xarray.Dataset
        dataset for variable wind components and wind speed at 10m
    ds_fgg: xarray.Dataset
        dataset for variable 10fgg25
    position: tuple
        tuple with the latitude and longitude of the position
    
    Returns:
    df_msl: pandas.DataFrame
        DataFrame with the temporal evolution of msl
    df_skt: pandas.DataFrame
        DataFrame with the temporal evolution of skt
    df_tp: pandas.DataFrame
        DataFrame with the temporal evolution of tp
    df_wind: pandas.DataFrame
        DataFrame with the temporal evolution of ws10
    df_fgg: pandas.DataFrame
        DataFrame with the temporal evolution of 10fgg25
    '''
    # Extract latitude and longitude from the position
    point_lat = position[0]
    point_lon = position[1]

    # Extract the data for the specified position
    point_ds = ds_msl.sel(latitude=point_lat, longitude=point_lon, method="nearest")
    df_msl = pd.DataFrame({'date': point_ds.valid_time.values, 'msl': point_ds.msl.values})
    point_ds = ds_skt.sel(latitude=point_lat, longitude=point_lon, method="nearest")
    df_skt = pd.DataFrame({'date': point_ds.valid_time.values, 'skt': point_ds.skt.values})
    point_ds = ds_tp.sel(latitude=point_lat, longitude=point_lon, method="nearest")
    df_tp = pd.DataFrame({'date': point_ds.valid_time.values, 'tp': point_ds.tp.values})
    point_ds = ds_wind.sel(latitude=point_lat, longitude=point_lon, method="nearest")
    df_wind = pd.DataFrame({'date': point_ds.valid_time.values, 'ws10': point_ds.ws10.values})
    point_ds = ds_fgg.sel(latitude=point_lat, longitude=point_lon, method="nearest")
    df_fgg = pd.DataFrame({'date': point_ds.valid_time.values, '10fgg25': point_ds.fg10g25.values})

    return df_msl, df_skt, df_tp, df_wind, df_fgg

def get_plot(ds_msl, ds_skt, ds_tp, ds_wind, ds_fgg, position):
    '''
    Creates the plots for the temporal evolution of the variables msl, skt, tp and 10fgg25 for the specified position.

    df_msl: pandas.DataFrame
        DataFrame with the temporal evolution of msl
    df_skt: pandas.DataFrame
        DataFrame with the temporal evolution of skt
    df_tp: pandas.DataFrame
        DataFrame with the temporal evolution of tp
    df_wind: pandas.DataFrame
        DataFrame with the temporal evolution of ws10
    df_fgg: pandas.DataFrame
        DataFrame with the temporal evolution of 10fgg25
    position: tuple
        tuple with the latitude and longitude of the position
    
    Returns:
    plt_tp: bqplot.Figure
        plot with the temporal evolution of total precipitation
    plt_msl: bqplot.Figure
        plot with the temporal evolution of mean sea level pressure
    plt_skt: bqplot.Figure
        plot with the temporal evolution of skin temperature
    plt_ws: bqplot.Figure
        plot with the temporal evolution of wind speed at 10m
    plt_fgg: bqplot.Figure
        plot with the temporal evolution of wind gust over 25 knots
    '''
    # Extract variables temporal evolution for the specified position
    df_msl, df_skt, df_tp, df_wind, df_fgg = get_df_pos(ds_msl, ds_skt, ds_tp, ds_wind, ds_fgg, position)
    
    # Define the scales for the plots and the characteristics of the lines and bars
    tp = df_tp['tp'].values
    x_tp = df_tp['date'].values
    lines_tp = Bars(x=x_tp, y=tp, colors='blue', scales={'x':DateScale(min=np.datetime64(x_tp[0]), max=np.datetime64(x_tp[-1])), 
                                                         'y': LinearScale(min=0, max=float(tp.max()))})
    msl = df_msl['msl'].values
    x_msl = df_msl['date'].values
    lines_msl = Lines(x=x_msl, y=msl, colors='orange', scales={'x':DateScale(min=np.datetime64(x_msl[0]), max=np.datetime64(x_msl[-1])), 
                                                               'y': LinearScale(min=float(msl.min()), max=float(msl.max()))})
    skt = df_skt['skt'].values
    x_skt = df_skt['date'].values
    lines_skt = Lines(x=x_skt, y=skt, colors='red', scales={'x':DateScale(min=np.datetime64(x_skt[0]),max=np.datetime64(x_skt[-1])),
                                                            'y': LinearScale(min=float(skt.min()), max=float(skt.max()))})
    ws = df_wind['ws10'].values
    x_ws = df_wind['date'].values
    lines_ws = Lines(x=x_ws, y=ws, colors='purple', scales={'x':DateScale(min=np.datetime64(x_ws[0]),max=np.datetime64(x_ws[-1])),
                                                            'y': LinearScale(min=float(ws.min()), max=float(ws.max()))})
    x_fgg = df_fgg['date'].values
    fgg = df_fgg['10fgg25'].values
    lines_fgg = Bars(x=x_fgg, y=fgg, colors='green', scales={'x':DateScale(min=np.datetime64(x_fgg[0]),max=np.datetime64(x_fgg[-1])),
                                                             'y': LinearScale(min=0, max=float(fgg.max()))})
    
    # Define the axes for the plots
    ax_tp = Axis(label='[mm]', scale=LinearScale(min=0, max=float(tp.max())),  orientation="vertical", side="left")
    ax_msl = Axis(label='[hPa]', scale=LinearScale(min=float(msl.min()), max=float(msl.max())), orientation="vertical", side="left")
    ax_skt = Axis(label='[°C]', scale=LinearScale(min=float(skt.min()), max=float(skt.max())), orientation="vertical", side="left")
    ax_ws = Axis(label='[m/s]', scale=LinearScale(min=float(ws.min()), max=float(ws.max())), orientation="vertical", side="left")
    ax_fgg = Axis(label='[%]', scale=LinearScale(min=0, max=float(fgg.max())), orientation="vertical", side="left")
    
    if len(x_tp) < 7:
        stp = 1
    elif len(x_tp) < 13:
        stp = 2
    elif len(x_tp) < 19:
        stp = 3
    else:
        stp = 4
    ax_x_tp = Axis(scale=DateScale(min=np.datetime64(x_tp[0]), max=np.datetime64(x_tp[-1])), label='Time',
                   num_ticks=len(x_tp[::stp]), tick_format='%d/%mh%H', orientation="vertical", side="bottom")
    ax_x_msl = Axis(scale=DateScale(min=np.datetime64(x_msl[0]), max=np.datetime64(x_msl[-1])), label='Time',
                    num_ticks=len(x_msl[::stp]), tick_format='%d/%mh%H', orientation="vertical", side="bottom")
    ax_x_skt = Axis(scale=DateScale(min=np.datetime64(x_skt[0]), max=np.datetime64(x_skt[-1])), label='Time',
                    num_ticks=len(x_skt[::stp]), tick_format='%d/%mh%H', orientation="vertical", side="bottom")
    ax_x_ws = Axis(scale=DateScale(min=np.datetime64(x_ws[0]), max=np.datetime64(x_ws[-1])), label='Time',
                   num_ticks=len(x_ws[::stp]), tick_format='%d/%mh%H', orientation="vertical", side="bottom")
    ax_x_fgg = Axis(scale=DateScale(min=np.datetime64(x_fgg[0]), max=np.datetime64(x_fgg[-1])), label='Time',
                    num_ticks=len(x_fgg[::stp]), tick_format='%d/%mh%H', orientation="vertical", side="bottom")

    # Create the plots
    plt_tp = Figure(axes=[ax_x_tp, ax_tp], title='Total precipitation', marks=[lines_tp], animation_duration=10,
                    layout={"max_width": "600px", "max_height": "400px"})
    plt_msl = Figure(axes=[ax_x_msl, ax_msl], title='Mean Sea Level Pressure', marks=[lines_msl], animation_duration=10,
                     layout={"max_width": "600px", "max_height": "400px"})
    plt_skt = Figure(axes=[ax_x_skt, ax_skt], title='Skin Temperature', marks=[lines_skt], animation_duration=10,
                     layout={"max_width": "600px", "max_height": "400px"})
    plt_ws = Figure(axes=[ax_x_ws, ax_ws], title='10m Wind Speed', marks=[lines_ws], animation_duration=10,
                    layout={"max_width": "600px", "max_height": "400px"})
    plt_fgg = Figure(axes=[ax_x_fgg, ax_fgg], title='10m Wind Gust over 25 knots', marks=[lines_fgg], animation_duration=10,
                     layout={"max_width": "600px", "max_height": "400px"})
    
    return plt_tp, plt_msl, plt_skt, plt_ws, plt_fgg

def create_maps_s5(initial_lat_lon, locations_avg, file_names):
    '''
    Creates the map with the plots of the temporal evolution of the variables msl, skt, tp and 10fgg25 
    for the specified position.

    initial_lat_lon: tuple
        tuple with the latitude and longitude of the initial position of the average forecast track
    locations_avg: list
        list with the latitude and longitude of the average forecast track
    file_names: list
        list with the names of the variable files to load

    Returns:
    m: ipyleaflet.Map
        map with the plots of the temporal evolution of the variables msl, skt, tp and 10fgg25
    '''
    # Load the datasets for the variables msl, skt, tp and 10fgg25
    ds_msl, ds_skt, ds_tp, ds_wind, ds_fgg = get_forecast_datasets(file_names)

    # Get the plots for the initial position
    plt_tp, plt_msl, plt_skt, plt_ws, plt_fgg = get_plot(ds_msl, ds_skt, ds_tp, ds_wind, ds_fgg, initial_lat_lon)

    # Create the map
    m = ipyleaflet.Map(center=[initial_lat_lon[0], initial_lat_lon[1]+50],
                       basemap=ipyleaflet.basemaps.Esri.WorldTopoMap,
                       zoom = 2)
    
    # Create the widget with the plots
    item_layout = widgets.Layout(overflow_y='scroll', width='550px', height='350px', flex_flow='column', display='block')
    main_figure = widgets.Box(children=[plt_tp, plt_msl, plt_skt, plt_ws, plt_fgg], layout=item_layout,  width='550px', height='550px')
    widget_control1 = ipyleaflet.WidgetControl(widget=main_figure, position="bottomright")
    m.add(widget_control1)

    # Create the marker to be dragged by the user
    marker = ipyleaflet.Marker(location=initial_lat_lon, draggable=True, name = 'Position') 
    m.add_layer(marker)

    # Function to update the plots when the marker is dragged
    def update_plots(new_lat, new_lon):
        plt_tp, plt_msl, plt_skt, plt_ws, plt_fgg = get_plot(ds_msl, ds_skt, ds_tp, ds_wind, ds_fgg, (new_lat, new_lon))
        main_figure.children = [plt_tp, plt_msl, plt_skt, plt_ws, plt_fgg]

    # Event handler called whenever the marker is dragged
    def on_marker_drag(event):
        new_location = event['new']  # The new location of the marker
        update_plots(new_location[0], new_location[1])

    # Add the event listener to the marker
    marker.observe(on_marker_drag, names='location')

    # Add the average forecast track to the map
    avg_track_antpath = ipyleaflet.AntPath(locations = locations_avg, color = "red")
    m.add_layer(avg_track_antpath)

    # Add the full screen control to the map
    m.add_control(ipyleaflet.FullScreenControl())

    return m
    