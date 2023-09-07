import ipywidgets as widgets
import pandas as pd
from IPython.display import display
from Magics.macro import * 

from ecmwf.opendata import Client
import numpy as np
import ipyleaflet
from ipyleaflet import Map, Marker

import ipywidgets
import xarray as xr
from scipy.spatial.distance import cdist

import matplotlib.pyplot as plt

from bqplot import Lines, Figure, LinearScale, Axis, DateScale
import bqplot
from datetime import datetime

### FUNCTIONS DEF
def closest_point(point, points):
    """ 
        Find closest point from a list of points. 
        Input: 
            - point: reference point 
            - points: list of points 
        Output: 
            - closest point 
    """
    return points[cdist([point], points).argmin()]

def match_value(df, col1, x, col2):
    """ 
        Match value x from col1 row to value in col2.
        Input: 
            - df: dataframe
            - col1: string; column name 1 of the dataframe
            - x: pattern no match
            - col2: string; column name 2 of the dataframe
        Output: 
            - Value matched 
    """
    return df[df[col1] == x][col2].values

def handle_move(data_allpoints, steps_to_download):
    def on_move(event, **kwargs):
        pos = kwargs['location']
        data_plot = get_df_pos(data_allpoints, pos, steps_to_download)
        lines1.y = data_plot['tp'].y.to_numpy()
        lines2.y = data_plot['msl'].y.to_numpy()
        lines3.y = data_plot['skt'].y.to_numpy()
        lines4.y = data_plot['10fgg15'].y.to_numpy()
    return on_move
    
def get_allvars_allpoints(date_to_download):
    """
        Read grib file of each variable; transform it to pandas dataframe;
            store dataframes into a dictionary with the key being the variable names. 
        Input: 
            - date_to_download: string (YYYYMMDD) or integer (0/1/-1) indicating the 
                starting date of the forecast
        Output: 
            - out: dictionary containing the pandas dataframe of each variable
    """
    out = {}
    for var in ['tp', 'msl', 'skt', '10fgg15']:
        
        filename = 'data/' + var + '_' + str(date_to_download) + '.grib'
        ds = xr.open_dataset(filename, engine="cfgrib")

        df = ds.to_dataframe()
        lats = df.index.get_level_values("latitude")
        lons = df.index.get_level_values("longitude")
        
        if var == '10fgg15':
            vals = df['fg10g15']
        else:
            vals = df[var]
        # DATAFRAME OF ALL POINTS FROM MAP
        df_datapoints = pd.DataFrame({'Lat': lats, 'Lon':lons, 'Value': vals})
        df_datapoints['point'] = [(x, y) for x,y in zip(df_datapoints['Lat'], df_datapoints['Lon'])]
        out[var] = df_datapoints
    return(out)

def get_df_pos(data_allpoints, pos, steps_to_download):
    """
        For each variable, for a given point coordinate, returns the temporal evolution of the 
            variable in the given Lat, Lon coordinates.  
        Input: 
            - data_allpoints: dictionary containing a pandas dataframe on each key; key are the variable name
            - pos: lat lon tuple
        Output: 
            - dictionary containing a pandas dataframe for each variable
    """
    data_pos_allvars = {}
    for var in ['tp', 'msl', 'skt', '10fgg15']:

        df_datapoints = data_allpoints[var]
        
        df_selpoint = pd.DataFrame({'Lat': pd.Series(pos[0]), 'Lon':pd.Series(pos[1])})
        df_selpoint['point'] = [(x, y) for x,y in zip(df_selpoint['Lat'], df_selpoint['Lon'])]
        
        # IDENTIFY THE CLOSEST POINT
        df_selpoint['closest'] = [closest_point(x, list(df_datapoints['point'])) for x in df_selpoint['point']]

        # GET THE VALUE OF VAR FROM THE CLOSEST POINT FOR ALL TIMESTEPS:
        vec_values = [match_value(df_datapoints, 'point', x, 'Value') for x in df_selpoint['closest']]

        # RETURN DATAFRAME TO PLOT
        df_pos_var = pd.DataFrame({'x': pd.Series(steps_to_download), 'y':pd.Series(vec_values[0])})
        data_pos_allvars[var] = df_pos_var
    return(data_pos_allvars)

def get_initial_plot(data_allpoints, initial_date, initial_latlon, steps_to_download):
    '''Input: 
        - data_allpoints: df of global data
        - initial_latlon: list with the initial coordinates
        - steps_to_download: list of steps to download
        Output: 
        - p1, p2, p3, p4: one plot for each variable of: total precipitation, msl,
                          skin temperature, prob. of wind gust at 1m of > 15m/s'''
    global lines1, lines2, lines3, lines4
    
    data_initial_plot = get_df_pos(data_allpoints, initial_latlon, steps_to_download)

    
    xdata1 = data_initial_plot['tp'].x.to_numpy()
    ydata1 = data_initial_plot['tp'].y.to_numpy()
    lines1 = bqplot.Bars(x=[xdata1], y=[ydata1], 
              scales={"x": LinearScale(min=float(min(xdata1)), max=float(max(xdata1))),
                      "y": LinearScale(min=0, max=float(max(ydata1)))})
    ax_y1 = Axis(label='mm',
                 scale=LinearScale(min=0, max=float(max(ydata1))), 
                 orientation="vertical", side="left")

    xdata2 = data_initial_plot['msl'].x.to_numpy()
    ydata2 = data_initial_plot['msl'].y.to_numpy()
    lines2 = Lines(x=[xdata2], y=[ydata2], 
              scales={"x": LinearScale(min=float(min(xdata2)), max=float(max(xdata2))),
                      "y": LinearScale(min=950, max=1040)})
    ax_y2 = Axis(label='hPa', 
                 scale=LinearScale(min=0, max=float(max(ydata2))), 
                 orientation="vertical", side="left")
    
    
    xdata3 = data_initial_plot['skt'].x.to_numpy()
    ydata3 = data_initial_plot['skt'].y.to_numpy()
    lines3 = Lines(x=[xdata3], y=[ydata3], 
              scales={"x": LinearScale(min=float(min(xdata3)), max=float(max(xdata3))),
                      "y": LinearScale(min=0, max=float(max(ydata3)))})
    ax_y3 = Axis(label='ºC', 
                 scale=LinearScale(min=float(min(ydata3)), max=float(max(ydata3))), 
                 orientation="vertical", side="left")
    
    xdata4 = data_initial_plot['10fgg15'].x.to_numpy()
    ydata4 = data_initial_plot['10fgg15'].y.to_numpy()
    lines4 = bqplot.Bars(x=[xdata4], y=[ydata4], 
              scales={"x": LinearScale(min=float(min(xdata4)), max=float(max(xdata4))),
                      "y": LinearScale(min=0, max=float(max(ydata4)))})
    ax_y4 = Axis(label='%', 
                scale=LinearScale(min=0, max=float(max(ydata4))), 
                 orientation="vertical", side="left")

    titol = 'Hours from ' + str(datetime.strptime(initial_date, '%Y%m%d').date().strftime("%Y-%m-%d"))
    ax_x = Axis(label=titol,
            scale=LinearScale(min=float(steps_to_download[0]), max=float(steps_to_download[-2])),
            num_ticks=int(len(steps_to_download)*0.5))
    

    p1 = Figure(
            axes=[ax_x, ax_y1],
            title='Total precipitation',
            marks=[lines1],
            animation_duration=10,
            layout={"max_width": "350px", "max_height": "350px"},
        )

    p2 = Figure(
            axes=[ax_x, ax_y2], title='Mean sea level pressure',
            marks=[lines2], animation_duration=10,
            layout={"max_width": "350px", "max_height": "350px"},
        )

    p3 = Figure(
            axes=[ax_x, ax_y3], title='Skin temperature',
            marks=[lines3], animation_duration=10,
            layout={"max_width": "350px", "max_height": "350px"},
        )

    p4 = Figure(
            axes=[ax_x, ax_y4],title='Wind gust 10m > 15m/s',
            marks=[lines4], animation_duration=10,
            layout={"max_width": "350px", "max_height": "350px"},
        )

    return p1, p2, p3, p4

def download_data_s5(date_to_download, steps_to_download, steps_to_download2):
    ### DOWNLOAD DATA
    '''Input: 
        - date_to_download: string, date to download in the format %Y%m%d
        - steps_to_download: list of strings with the steps to download for tp, msp, skn
        - steps_to_download2: list of strings with the steps to download for wind prob.
        '''
    print('Downloading data...')
    client = Client("azure", beta=True) # ecwf: last five days
    for var in ['tp', 'msl', 'skt']:
        print(var)
        client.retrieve(
            date = date_to_download, #date start of the forecast
            time = 0,  # time start of the forecast or 12
            step = steps_to_download, #step of the forecast
            stream = "oper",
            type = "fc",
            levtype = "sfc",
            param = var,
            target = 'data/' + var + '_' + str(date_to_download) + '.grib',
        )

    var = '10fgg15'
    client.retrieve(
            date = date_to_download, #date start of the forecast
            time = 0,  # time start of the forecast or 12
            step = steps_to_download2, #step of the forecast
            stream = "enfo",
            type = "ep",
            param = var,
            target = 'data/' + var + '_' + str(date_to_download) + '.grib',
        )    

def load_data_s5(date_to_download, filename_avg_track):
    '''Input: 
    - date_to_download: string, date to download in the format %Y%m%d
    - filename_avg_track: string
    '''
    global data_allpoints
    print('Loading data...')
    grib_data = {}
    for var in ['tp', 'msl', 'skt', '10fgg15']:

        filename = 'data/' + var + '_' + str(date_to_download) + '.grib'
        grib_data[var] = mgrib(grib_input_file_name=filename)

    data_allpoints = get_allvars_allpoints(date_to_download)
    data_allpoints['skt']['Value'] = data_allpoints['skt']['Value'] - 273.15
    data_allpoints['msl']['Value'] = data_allpoints['msl']['Value']/100

    # LOAD AVG TRACK
    df_avg_track = pd.read_csv(filename_avg_track)
    return(data_allpoints, df_avg_track)

# CREATE MAP 

def map_s5(initial_latlon, initial_date, final_date, filename_avg_track):
    '''Input: 
    - initial_date, final_date: string of the date in format %Y%m%d
    - filename_avg_track: string
    Output: map
    '''
    global steps_to_download
    ### VARIABLE INITIALIZATION
    cyclone_days = pd.date_range(start=initial_date, end=final_date)
    date_to_download = initial_date
    steps_to_download = list(np.arange(0,240, 12)[0:2*len(cyclone_days)])
    print(steps_to_download)
    steps_to_download2 = [str(i)+'-'+str(i+24) for i in steps_to_download]
    # download data: 
    download_data_s5(date_to_download, steps_to_download, steps_to_download2)

    ### LOAD DATA
    data_allpoints, df_avg_track =load_data_s5(date_to_download, filename_avg_track)
    print('Printing map...')
    m = Map(
        center=[initial_latlon[0], initial_latlon[1]+50],
        basemap=ipyleaflet.basemaps.OpenStreetMap.France,
        zoom = 2,
    )

    p1, p2, p3, p4 = get_initial_plot(data_allpoints, initial_date, initial_latlon, steps_to_download)
    
    item_layout = ipywidgets.Layout(overflow_y='scroll', width='350px', height='350px',
                               flex_flow='column', display='block')
    main_figure = widgets.Box(children=[p1, p2, p3, p4], layout=item_layout,  width='350px', height='350px')

    #header = ipywidgets.HTML("<h1>Fictional World Weather</h1>", layout=ipywidgets.Layout(height='auto'))
    #header.style.text_align='center'
    #layout_plots = ipywidgets.AppLayout(header=header, center=main_figure)
    widget_control1 = ipyleaflet.WidgetControl(widget=main_figure, position="bottomright")
    #widget_control1 = ipyleaflet.WidgetControl(widget=layout_plots, position="bottomright")

    marker = Marker(location=initial_latlon, draggable=True, name = 'Position') 
    m.add_layer(marker)
    marker.on_move(handle_move(data_allpoints, steps_to_download))


    m.add(widget_control1)
    avg_track = ipyleaflet.AntPath(locations = df_avg_track.values.tolist(), color = "red")
    m.add_layer(avg_track)
    return m