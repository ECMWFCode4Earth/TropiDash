# TropiDash Documentation for User support

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

<!-- code_chunk_output -->

- [TropiDash Documentation for User support](#tropidash-documentation-for-user-support)
  - [Introduction](#introduction)
  - [Section 1 - Tropical Cyclone Tracks and Characteristics](#section-1---tropical-cyclone-tracks-and-characteristics)
    - [Data sources](#data-sources)
    - [Data structure](#data-structure)
    - [Plot](#plot)
    - [Functions](#functions)
  - [Section 2 - Atmospherical variables](#section-2---atmospherical-variables)
    - [Data sources](#data-sources-1)
    - [Data manipulation](#data-manipulation)
    - [Plot](#plot-1)
    - [Functions](#functions-1)
  - [Section 3 - Impacts variables](#section-3---impacts-variables)
    - [Data sources](#data-sources-2)
    - [Plot](#plot-2)
    - [Functions](#functions-2)
  - [Section 4 - Joint visualization](#section-4---joint-visualization)
  - [Data sources](#data-sources-3)
  - [Functions](#functions-3)
  - [Section 5 - Point-Wise Temporal Evolution of Atmospheric Variable](#section-5---point-wise-temporal-evolution-of-atmospheric-variable)
    - [Data sources](#data-sources-4)
    - [Plot](#plot-3)
    - [Functions](#functions-4)

<!-- /code_chunk_output -->


## Introduction

Welcome to the User Support Documentation for **TropiDash**, an innovative dashboard designed for the visualization and display of tropical cyclone data. Developed using Python within a Jupyter Notebook environment, TropiDash harnesses the power of coding to deliver dynamic insights into tropical cyclones hazard and impacts. As part of the Code for Earth initiative by the European Centre of Medium-Range Weather Forecast (ECMWF), this dashboard represents a collaborative effort to provide accessible and actionable information about tropical cyclones. TropiDash is organized into five distinct sections, each offering unique perspectives on tropical cyclones:

1. Tropical Cyclone Tracks and Characteristics 
2. Atmospheric Variables 
3. Impact Variables
4. Multi-Layer Map combining data from all previous sections
5. Point-Wise Temporal Evolution of Atmospheric Variables

Whether you're a weather enthusiast, a researcher, or a decision-maker, this documentation will guide you through using TropiDash effectively.
To support you for an effective use of the dashboard please also refer to the tutorials specific for each section in jupyter notebook. You can find them in the __[tutorials folder](https://github.com/ECMWFCode4Earth/TropiDash/tree/main/tutorials)__.

## Section 1 - Tropical Cyclone Tracks and Characteristics

In this section, the goal is to provide users with a comprehensive view of tropical cyclones by visualizing both observed and forecasted tracks and crucial additional information. Here, you'll find four main layers: the Ensemble Forecast Tracks, the Observed Track, the Average Forecast Track, and the Strike Probability Map. These layers collectively offer a comprehensive understanding of cyclone behavior, enabling you to assess the potential impact of these powerful weather phenomena. Whether you're tracking the historical path of a cyclone, examining the consensus of future forecasts, or gauging the likelihood of a strike in a particular area, this section equips you with the tools you need to stay informed about cyclones' positions and future developments.

### Data sources
The data for Section 1 is sourced from two datasets. The forecasts are derived from the open dataset provided by the European Centre of Medium-Range Weather Forecasts (ECMWF), accessed with the Azure client through the ecmwf.opendata library. You can find more information regarding the open-data __[here](https://www.ecmwf.int/en/forecasts/datasets/open-data)__. In parallel, the observed track data is obtained from the International Best Track Archive for Climate Stewardship (IBTrACS) dataset, maintained by the National Oceanic and Atmospheric Administration (NOAA). You can find more information regarding the IBTrACS dataset __[here](https://www.ncei.noaa.gov/products/international-best-track-archive)__. 

### Data structure
The ECMWF data contains the ensemble forecast tracks of the active cyclones with a temporal resolution of 6 hours. The information is saved in a BUFR file, which can be imported into Jupyter Notebook as pandas dataframes for detailed analysis and visualization. The dataframe comprises a wealth of information, including the storm identifier and name, ensemble member number, latitude and longitude coordinates of the cyclone center, and crucial temporal details such as the year, month, day, and hour of the forecast initiation. Additionally, each row of the dataframe specifies the time period from the forecast date, i.e., hours passed, allowing users to track cyclone development over time. Key meteorological parameters such as the central pressure in *Pascal* and maximum sustained wind speed in *Metre per second* within the cyclone system are also included. The presence of "false" storms in the BUFR file is worth noting, identifiable by storm identifiers with numbers greater than 70.

The IBTrACS data contains the observed track of the active cyclones with a temporal resolution of 3 hours. These observations are stored in CSV files, providing accessibility and ease of use for users. The CSV file is structured with multiple columns, each contributing essential information for comprehensive analysis. These columns include the storm identifier, year, cardinal number of the storm for that season, basin, subbasin (if applicable), storm name as designated by the relevant agency, timestamp in UTC, storm nature, latitude, and longitude coordinates of the cyclone center. Additionally, the data contains the maximum sustained wind speed reported by the World Meteorological Organization (WMO) for the current location, the responsible reporting agency, track type, distance to the nearest land from the current position, the closest predicted landfall location within the next 6 hours, and an interpolation flag. It's important to note that IBTrACS cyclone track data undergo post-processing after the cyclone dissipates, often resulting in incomplete information, particularly for the latest observations. 

It is important to acknowledge that within TropiDash, users may encounter a temporal discrepancy between the latest observed cyclone position in the IBTrACS dataset and the initial forecasted position within the ECMWF data. This divergence arises from the inherent nature of data acquisition and processing. The IBTrACS dataset relies on the meticulous collection and verification of observational data, necessitating time. Consequently, the most recent observations often pertain to roughly 2 to 3 days before the current date. In contrast, the ECMWF forecast data offers a forward-looking perspective, initiating from the present and projecting into the future. This time lag between observational and forecast data underscores the importance of considering the temporal context when interpreting and utilizing cyclone information within TropiDash.

### Plot
Considering the information available from the presented dataset, we organize the dashboard section to display the following layers over a geographical map:
1. *Ensemble forecast tracks*: this layer groups all the ensemble members of the forecast track for the selected cyclone. Suppose the user plots five or fewer ensemble members for each forecasted position. In that case, clicking on the location displays the predicted mean sea level pressure at the hurricane's center and the maximum sustained wind speed. With six or more ensemble members selected, this option is not available anymore because it makes the coding process heavier. 
2. *Average forecast track*: this layer displays the mean forecast track computed considering all the 52 ensemble members included in the forecast. Clicking on the mean locations, the user can also visualize the 10<sup>th</sup>, 25<sup>th</sup>, 50<sup>th</sup>, 75<sup>th</sup>, and 90<sup>th</sup> percentiles for mean sea level pressure and maximum sustained wind speed of the ensemble forecasted tracks. 
3. *Observed track*: this layer displays the observed track of the cyclone. 
4. *Strike probability map*: this layer displays the map of the strike probability of the cyclone. The strike probability is the probability that a tropical cyclone will pass within a 300 km radius of a given location and within a time window of 48 hours. It quickly assesses high-risk areas, although with some uncertainty in the exact timing or position. Please find more information regarding the strike probability __[here](https://charts.ecmwf.int/products/medium-tc-genesis?base_time=202309140000&layer_name=genesis_ts&projection=opencharts_global&valid_time=202309170000)__.

### Functions

To better understand how the data are processed and how the interactive plots are produced please refer to the python script containing the functions for this section: __[TropiDash/utils_tracks.py](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash/utils_tracks.py)__.

## Section 2 - Atmospherical variables

This Secton's goal is to provide the user the possibility to visualize forecasts of atmospheric variables related to cyclone formation and forecasting. The variables are plotted at user-selected time steps ahead from the user-selected date. The tutorial for this Section usage is provided at __[tutorials/Section2_Atmospheric_Variables_Tutorial.ipynb](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/tutorials/Section2_Atmospheric_Variables_Tutorial.ipynb)__.

### Data sources

The atmosperic data source is the open dataset provided by the European Centre of Medium-Range Weather Forecasts (ECMWF), accessed with the Azure client through [`ecmwf.opendata`](https://github.com/ecmwf/ecmwf-opendata) package. You can find more information regarding the open-data __[here](https://www.ecmwf.int/en/forecasts/datasets/open-data)__ 

Data products used here are based on the medium-range (high-resolution and ensemble) forecast model which is released 1 hour after the real-time dissemination schedule at 0.4 degrees resolution. The downloaded data in this section corresponds to the forecast released on the starting date of the cyclone. The temporal steps displayed are 1 every 12 hours for each day from the starting  date of the selected cyclone until the last day of the cyclone.

The variables displayed in this section are the accumulated precipitation, mean sea level pressure, skin temperature, wind speed and direction, and probability of wind gusts of more than 25 m/s at 10 meters.

### Data manipulation

The original data files downloaded from ECMWF's Open Data catalog are in GRIB format. This format has been changed to .tif through `xarray`, `rioxarray`, `cfgrib` and `rasterio` packages to be able to plot the data through `localtilesever` as an ipyleaflet's `Tile layer`. The only operation done on the data itself was to adjust the measuring units to better understand the data (e.g. from Pa to hPa).

### Plot

The plot is an interactive map deployed through [`ipyleaflet`](https://ipyleaflet.readthedocs.io/en/latest/index.html) and showing a list of layers. The layers can be selected and de-selected, the plot can be panned and zoomed and visualized in full-screen mode. Each variable shown is a **forecasted** variable. The available layers are:
 - Mean sea level pressure: ipyleaflet's `Tile layer` showing the global mean sea level pressure in hPa;
 - Skin temperature: ipyleaflet's `Tile layer` showing gloabl surface temperarure expressed in °C;
 - Total Precipitation: ipyleaflet's `Tile layer` showing worldwide total precipitation in m;
 - Probability of 10 metre wind gust of at least 25 m/s: ipyleaflet's `Tile layer` showing gloablly the probability expressed in percentage to have wind gust of at least 25 m/s at ten meters.
 - 10 metre wind component: ipyleaflet's `Velocity layer` showing the wind speed and direction.

A widget is shown to make the user select the forecasting step, meaning how many hours ahead from the selected date will the variable shown be forecasted.

Colorbars are provided in the top right corner to be able to interpret the data shown.

### Functions

The function used to manage the data download, manipulation and loading process as well as the plotting function are written and documented in: __[TropiDash/utils_atm.py](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash/utils_atm.py)__.

## Section 3 - Impacts variables

This Section's goal is to provide the user with a plot containing variables which can help their understanding of the possible impacts the selected cyclone may generate (e.g. in terms of population impacted), as well as providing risk maps and exposition indexes. The tutorial for this Section usage is provided at [tutorials/Section3_Impact_Variables_Tutorial.ipynb](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/tutorials/Section3_Impact_Variables_Tutorial.ipynb)

### Data sources

Data for this Section is downloaded from the following sources:
- Coastal hazard for different return periods: https://datacatalog.worldbank.org/search/dataset/0038579/Global-coastal-flood-hazard
- Cyclone hazard for different return periods: https://datacatalog.worldbank.org/search/dataset/0038577/Global-cyclone-hazard
- Population count (1 x 1 km cells): https://hub.worldpop.org/geodata/summary?id=24777
<!-- - Exposition Indexes: https://data.humdata.org/dataset/1efb6ee7-051a-440f-a2cf-e652fecccf73 -->

Population count raster was resampled to pixels of 10x10 km using QGIS, then passed through a function explained in [tutorials/X1_Download_population_data.ipynb](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/tutorials/X1_Download_population_data.ipynb) ("Fix nodata issue" section).

### Plot

The plot is an interactive map deployed through [`ipyleaflet`](https://ipyleaflet.readthedocs.io/en/latest/index.html) and showing a list of layers. The layers can be selected and de-selected, the plot can be panned and zoomed and visualized in full-screen mode. The available layers are:
 - Cyclone hazard: ipyleaflet's `Tile layer` showing cyclone hazard for a selected return period;
 - Coastal hazard: ipyleaflet's `Tile layer` showing coastal hazard for a selected return period;
 - Population: ipyleaflet's `Tile layer` showing population count in 10 km x 10 km cells;
 - Tsunamis, Coastal floods and Sea level rise Exposition Indexes: ipyleaflet's `Choropleth layer` showing each country's exposition index.

Two widgets are shown to make the user choose which return period to assign to the layers of cyclone and coastal hazards.

Colorbars are provided in the top right corner to be able to interpret the data shown.

### Functions

The function used to manage the data download, their loading processes as well as the plotting functions are written and documented in: __[TropiDash/utils_impacts.py](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash/utils_impacts.py)__.

## Section 4 - Joint visualization

To make the user able to visualize easier cyclone, atmospheric and impact variables this Section shows a selection of variables already shown in Section 1, 2 and 3, proiding all the available widgets. Below, the list of available data shown in this Section's plot:
- Cyclone's average forecasted track
- Strike probability map
- All atmospheric variables of Section 2
- Population count
- Cyclone and coastal hazard

The tutorial for this Section usage is provided at [tutorials/Section4_Joint_Visualization_tutorial.ipynb](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/tutorials/Section4_Joint_Visualization_tutorial.ipynb)

### Data sources

The sources for cylcone variables are listed in Section 1, the ones for atmospheric variables are listed in Section 2, and the ones for impacts variables are lister in Section 3.

### Plot

The plot is an interactive map deployed through [`ipyleaflet`](https://ipyleaflet.readthedocs.io/en/latest/index.html), which shows a combination of the layers shown in the previous Sections, so cyclone tracks and strike probability as well as atmospheric and hazard variables. 

### Functions

The plotting function is provided at __[TropiDashutils_section4.py](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash/utils_section4.py)__. It makes use of functions defined for the other Sections.

## Section 5 - Point-Wise Temporal Evolution of Atmospheric Variable

This section introduces the possibility to visualize the forecasted evolution of key atmospheric variables over time. Users are empowered to generate plots for a selection of critical variables, including accumulated precipitation, mean sea level pressure, skin temperature, and the probability of wind gusts exceeding 25 m/s at a 10-meter height. This functionality is designed to offer insights into the forecasted conditions at any user-specified location on the map, providing a tool for detailed weather analysis.

The tutorial for this Section is provided at [tutorials/Section5_Temporal_Evolution_Tutorial.ipynb](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/tutorials/Section4_Temporal_Evolution_Tutorial.ipynb)

### Data sources

The variables plotted in this section are sourced from the ECMWF's Open Data catalog, as detailed in [Section 2](#section-2---atmospherical-variables). This section includes data on accumulated precipitation, mean sea level pressure, skin temperature, and the probability of experiencing wind gusts exceeding 25 m/s at a height of 10 meters.

### Plot

Within this section, users are presented with an interactive background map that includes a movable marker. By positioning this marker, users can visually explore the temporal progression of the four atmospheric variables at the chosen location: accumulated precipitation, mean sea level pressure, skin temperature, and the probability of wind gusts surpassing 25 m/s at a height of 10 meters. This exploration spans from the inception to the conclusion of the selected cyclone event. A distinctive red-dashed line illustrates the average trajectory derived from ensemble forecasts.

As the marker is relocated across the map, the plots for these variables are dynamically updated, ensuring the most relevant data is always at the user's fingertips. Initially, the plot displayed will show the daily accumulated precipitation; further variables can be viewed by scrolling within the dedicated white box area.

This interactive experience is powered by [`ipyleaflet`](https://ipyleaflet.readthedocs.io/en/latest/index.html), a python library that delivers a user-friendly and engaging mapping interface.

### Functions
To better understand how the data are processed and how the interactive plots are produced please refer to the python script containing the functions for this section: __[TropiDash/utils_TemporalEvolution.py](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash/utils_temporal.py)__.
