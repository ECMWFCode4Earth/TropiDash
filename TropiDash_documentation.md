# TropiDash Documentation for User support

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

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

To better understand how the data are processed and how the interactive plots are produced please refer to the python script containing the functions for this section: __[utils_tracks.py](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash/utils_tracks.py)__.

## 2. Section 2 - Atmospherical variables

### Data sources

### Plot

### Functions

## 3. Section 3 - Impacts variables

### Data sources

### Plot

### Functions

## Section 5 - Point-Wise Temporal Evolution of Atmospheric Variable

introducció

### Data sources

Data source is the open dataset provided by the European Centre of Medium-Range Weather Forecasts (ECMWF), accessed with the Azure client through the ecmwf.opendata library. You can find more information regarding the open-data __[here](https://www.ecmwf.int/en/forecasts/datasets/open-data)__ 

The data products used here are based on the medium-range (high-resolution and ensemble) forecast model which is released 1 hour after the real-time dissemination schedule at 0.4 degrees resolution. The downloaded data in this section corresponds to the forecast released on the starting date of the cyclone. The temporal steps displayed are 1 every 12 hours for each day from the starting  date of the selected cyclone until the last day of the cyclone. The variables displayed in this section are the accumulated precipitation, mean sea level pressure, skin temperature, and probability of wind gusts of more than 25m/s at 10m. 

### Plot

In this section, users can find a background map, a marker and the temporal evolution on the marker position of the accumulated precipitation, mean sea level pressure, skin temperature, and probability of wind gusts of more than 25m/s at 10m, from the first day of the selected cyclone until the last one. The red-lashed line represents the average track of the tracks given by the ensembles. The plots are automatically uploaded every time the user moves the pointer location. Please notice that it may take a few seconds to upload the plots. The first variable displayed is the daily accumulated precipitation and the other ones appear on scrolling down on the white box. 

### Functions
To better understand how the data are processed and how the interactive plots are produced please refer to the python script containing the functions for this section: __[utils_tracks.py](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash/utils_tracks.py)__.
