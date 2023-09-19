# TropiDash: towards a comprehensive tropical cyclone hazard dashboard

ECMWF Code for Earth 2023 project

Run the Dashboard directly on Binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ECMWFCode4Earth/TropiDash/HEAD?urlpath=voila%2Frender%2FTropiDash%2FTropiDash_backbone.ipynb)

## Team
Filippo Dainelli - filippo.dainelli@polimi.it\
Paolo Colombo - paolo1.colombo@polimi.it\
Laura Paredes i Fortuny - paredes_laufor@externos.gva.es

## Project abstract

This project will produce a platform on Jupyter Notebooks able to visualize key meteorological parameters in plots and maps to understand tropical cyclone hazards evolution better. It will gather the currently used and most effective visualizations and reproduce them in a dashboard by applying interactive elements. The end user will be provided with sound documentation which will enable the platform usage and enhancement after the end of the project. The development of the project is envisioned following 4 major phases:
1. Existing context recognition
2. Platform design
3. Platform development
4. Ensure the future maintenance

## Repository content
- **TropiDash folder**: it contains the main script of the project ([TropiDash_backbone.ipynb](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash/TropiDash_backbone.ipynb)), the libraries created specifically for the project that are imported in TropiDash_backbone.ipynb (utils_NameSection.ipynb), and the folder 'data' in which data will be stored when running the main script.
- **tutorials** folder: it contains one tutorial for each Dashboard' Section (i.e. Sections 1 to 5) plus another about downloading population data used in Section 3. It also contains the folder 'data' where the data needed to run the tutorials is already downloaded.
- **TropiDash Documentation for User Support** ([TropiDash_documentation.md](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash_documentation.md)): detailed description of the data sources, data structure, displayed plots and specific functions of each Dashboard' Section (i.e. Sections 1 to 5).
- **Package requirements** ([requirements.txt](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/requirements.txt)): list of necessary packages to run TropiDash_backbone.ipynb in local. 
