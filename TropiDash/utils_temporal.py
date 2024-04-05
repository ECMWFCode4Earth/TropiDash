from bqplot import Lines, Bars, Figure, LinearScale, Axis, DateScale
import ipyleaflet
import ipywidgets as widgets
import pandas as pd
import numpy as np
import xarray as xr

def get_forecast_datasets(initial_date, stepsdict, variables):
    '''
    Loads the forecast of the variables in the variables list for the specified forecast date and returns a dataset 
    with the variables data.
    The variables in the list are:
    - Mean sea level pressure
    - Skin temperature
    - Total precipitation
    - 10m Wind gusts over 25 knots

    initial_date: str
        date of the forecast in format %Y%m%d
    stepsdict: dict
        dictionary with the forecasts steps for each variable
    variables: list
        list with the variables to load
    
    Returns:
    ds_3var: xarray.Dataset
        dataset with the variables msl, skt and tp
    ds_fgg: xarray.Dataset
        dataset with the variable 10fgg25
    '''
    # Load the datasets for the variables msl, skt and tp
    steps = stepsdict["base"]
    var = variables[0]
    filename = f"data/atm/{initial_date}/{var}_{initial_date}_time0_steps{steps[0]}-{steps[-1]}_oper_fc.grib"
    ds_msl = xr.open_dataset(filename, engine="cfgrib")
    var = variables[1]
    filename = f"data/atm/{initial_date}/{var}_{initial_date}_time0_steps{steps[0]}-{steps[-1]}_oper_fc.grib"
    ds_skt = xr.open_dataset(filename, engine="cfgrib")
    var = variables[2]
    filename = f"data/atm/{initial_date}/{var}_{initial_date}_time0_steps{steps[0]}-{steps[-1]}_oper_fc.grib"
    ds_tp = xr.open_dataset(filename, engine="cfgrib")
    
    # Create a new dataset with mean sea level pressure, skin temperature and total precipitation
    ds_3var = ds_msl.copy()
    ds_3var = ds_3var.assign(tp=ds_tp.tp)
    ds_3var = ds_3var.assign(skt=ds_skt.skt - 273.15) # Convert to Celsius
    ds_3var = ds_3var.assign(msl=ds_msl.msl / 100) # Convert to hPa
    
    # The dataset containting the wind gust over 25 knots has a different temporal period (forecasts starts 12 hours later)
    var = variables[3]
    steps = stepsdict["10fgg25"]
    filename = f"data/atm/{initial_date}/{var}_{initial_date}_time0_steps{steps[0]}-{steps[-1]}_enfo_ep.grib"
    ds_fgg = xr.open_dataset(filename, engine="cfgrib")

    return ds_3var, ds_fgg

def get_df_pos(ds_3var, ds_fgg, position):
    '''
    Extracts the data for the variables in the datasets ds_3var and ds_fgg for the specified position and returns
    the data in 2 different DataFrame.

    ds_3var: xarray.Dataset
        dataset with the variables msl, skt and tp
    ds_fgg: xarray.Dataset
        dataset with the variable 10fgg25
    position: tuple
        tuple with the latitude and longitude of the position
    
    Returns:
    df_3var: pandas.DataFrame
        DataFrame with the temporal evolution of variables msl, skt and tp
    df_fgg: pandas.DataFrame
        DataFrame with the temporal evolution of variable 10fgg25
    '''
    point_lat = position[0]
    point_lon = position[1]
    point_ds_3var = ds_3var.sel(latitude=point_lat, longitude=point_lon, method="nearest")
    point_ds_fgg = ds_fgg.sel(latitude=point_lat, longitude=point_lon, method="nearest")
    df_3var = pd.DataFrame({'date': ds_3var.valid_time.values, 'msl': point_ds_3var.msl.values, 'tp': point_ds_3var.tp.values, 
                            'skt': point_ds_3var.skt.values})
    df_fgg = pd.DataFrame({'date': ds_fgg.valid_time.values, '10fgg25': point_ds_fgg.fg10g25.values})

    return df_3var, df_fgg

def get_plot(ds_3var, ds_fgg, position):
    '''
    Creates the plots for the temporal evolution of the variables msl, skt, tp and 10fgg25 for the specified position.

    ds_3var: xarray.Dataset
        dataset with the variables msl, skt and tp
    ds_fgg: xarray.Dataset
        dataset with the variable 10fgg25
    position: tuple
        tuple with the latitude and longitude of the position
    
    Returns:
    plt_tp: bqplot.Figure
        plot with the temporal evolution of total precipitation
    plt_msl: bqplot.Figure
        plot with the temporal evolution of mean sea level pressure
    plt_skt: bqplot.Figure
        plot with the temporal evolution of skin temperature
    plt_fgg: bqplot.Figure
        plot with the temporal evolution of wind gust over 25 knots
    '''
    # Extract variables temporal evolution for the specified position
    df_3var, df_fgg = get_df_pos(ds_3var, ds_fgg, position)
    
    # Define the scales for the plots and the characteristics of the lines and bars
    xdata_3vars = df_3var['date'].values
    tp = df_3var['tp'].values
    msl = df_3var['msl'].values
    skt = df_3var['skt'].values
    lines_tp = Bars(x=xdata_3vars, y=tp, colors='blue',
                           scales={'x':DateScale(min=np.datetime64(xdata_3vars[0]), max=np.datetime64(xdata_3vars[-1])),
                                   'y': LinearScale(min=0, max=float(tp.max()))})
    lines_msl = Lines(x=xdata_3vars, y=msl, colors='orange',
                             scales={'x':DateScale(min=np.datetime64(xdata_3vars[0]), max=np.datetime64(xdata_3vars[-1])),
                                     'y': LinearScale(min=float(msl.min()), max=float(msl.max()))})
    lines_skt = Lines(x=xdata_3vars, y=skt, colors='red',
                             scales={'x':DateScale(min=np.datetime64(xdata_3vars[0]),max=np.datetime64(xdata_3vars[-1])),
                                     'y': LinearScale(min=float(skt.min()), max=float(skt.max()))})
    
    xdata_fgg = df_fgg['date'].values
    fgg = df_fgg['10fgg25'].values
    lines_fgg = Bars(x=xdata_fgg, y=fgg, colors='green',
                            scales={'x':DateScale(min=np.datetime64(xdata_fgg[0]),max=np.datetime64(xdata_fgg[-1])),
                                    'y': LinearScale(min=0, max=float(fgg.max()))})
    
    # Define the axes for the plots
    ax_tp = Axis(label='mm', scale=LinearScale(min=0, max=float(tp.max())),  orientation="vertical", side="left")
    ax_msl = Axis(label='hPa', scale=LinearScale(min=float(msl.min()), max=float(msl.max())), orientation="vertical", side="left")
    ax_skt = Axis(label='°C', scale=LinearScale(min=float(skt.min()), max=float(skt.max())), orientation="vertical", side="left")
    ax_fgg = Axis(label='%', scale=LinearScale(min=0, max=float(fgg.max())), orientation="vertical", side="left")
    
    ax_x_3var = Axis(scale=DateScale(min=np.datetime64(xdata_3vars[0]), max=np.datetime64(xdata_3vars[-1])), label='Time',
                     num_ticks=len(xdata_3vars[::3]), tick_format='%d/%mh%H', orientation="vertical", side="bottom")
    ax_x_fgg = Axis(scale=DateScale(min=np.datetime64(xdata_fgg[0]), max=np.datetime64(xdata_fgg[-1])), label='Time',
                     num_ticks=len(xdata_fgg[::3]), tick_format='%d/%mh%H', orientation="vertical", side="bottom")
    
    # Create the plots
    plt_tp = Figure(axes=[ax_x_3var, ax_tp], title='Total precipitation', marks=[lines_tp], animation_duration=10,
                    layout={"max_width": "600px", "max_height": "400px"})
    plt_msl = Figure(axes=[ax_x_3var, ax_msl], title='Mean Sea Level Pressure', marks=[lines_msl], animation_duration=10,
                     layout={"max_width": "600px", "max_height": "400px"})
    plt_skt = Figure(axes=[ax_x_3var, ax_skt], title='Skin Temperature', marks=[lines_skt], animation_duration=10,
                     layout={"max_width": "600px", "max_height": "400px"})
    plt_fgg = Figure(axes=[ax_x_fgg, ax_fgg], title='10m Wind Gust over 25 knots', marks=[lines_fgg], animation_duration=10,
                     layout={"max_width": "600px", "max_height": "400px"})
    
    return plt_tp, plt_msl, plt_skt, plt_fgg

def create_maps_s5(initial_lat_lon, initial_date, locations_avg, stepsdict, variables):
    '''
    Creates the map with the plots of the temporal evolution of the variables msl, skt, tp and 10fgg25 
    for the specified position.

    initial_lat_lon: tuple
        tuple with the latitude and longitude of the initial position of the average forecast track
    initial_date: str
        date of the forecast in format %Y%m%d
    locations_avg: list
        list with the latitude and longitude of the average forecast track
    stepsdict: dict
        dictionary with the forecasts steps for each variable
    variables: list
        list with the variables to load

    Returns:
    m: ipyleaflet.Map
        map with the plots of the temporal evolution of the variables msl, skt, tp and 10fgg25
    '''
    # Load the datasets for the variables msl, skt, tp and 10fgg25
    ds_3var, ds_fgg = get_forecast_datasets(initial_date, stepsdict, variables)

    # Get the plots for the initial position
    plt_tp, plt_msl, plt_skt, plt_fgg = get_plot(ds_3var, ds_fgg, initial_lat_lon)

    # Create the map
    m = ipyleaflet.Map(
            center=[initial_lat_lon[0], initial_lat_lon[1]+50],
            basemap=ipyleaflet.basemaps.Esri.WorldTopoMap,
            zoom = 2,
        )
    
    # Create the widget with the plots
    item_layout = widgets.Layout(overflow_y='scroll', width='450px', height='350px', flex_flow='column', display='block')
    main_figure = widgets.Box(children=[plt_tp, plt_msl, plt_skt, plt_fgg], layout=item_layout,  width='550px', height='550px')
    widget_control1 = ipyleaflet.WidgetControl(widget=main_figure, position="bottomright")
    m.add(widget_control1)

    # Create the marker to be dragged by the user
    marker = ipyleaflet.Marker(location=initial_lat_lon, draggable=True, name = 'Position') 
    m.add_layer(marker)

    # Function to update the plots when the marker is dragged
    def update_plots(new_lat, new_lon):
        plt_tp, plt_msl, plt_skt, plt_fgg = get_plot(ds_3var, ds_fgg, (new_lat, new_lon))
        main_figure.children = [plt_tp, plt_msl, plt_skt, plt_fgg]

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
    