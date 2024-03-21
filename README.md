# TropiDash: towards a comprehensive tropical cyclone hazard dashboard

TropiDash is a visualization tool designed to empower researchers, meteorologists, and enthusiasts with the ability to analyze and interpret tropical cyclone data effectively. Utilizing in web technologies and data visualization libraries, TropiDash provides a comprehensive overview of tropical cyclones, their paths, intensities, and potential impacts.
The project was developed as part of the ECMWF Code for Earth 2023 initiative. [Here](https://codeforearth.ecmwf.int/) you can fdin more details regarding Code for Earth. [Here](https://github.com/ECMWFCode4Earth/challenges_2023) you can find more details about the other projects of 2023.  

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

## Repository structure

- **TropiDash folder**: Main script ([TropiDash_backbone.ipynb](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash/TropiDash_backbone.ipynb)), utility libraries, and data storage folders.
- **Tutorials Folder**: Step-by-step guides for each dashboard section and data downloading instructions.
- **Documentation** ([TropiDash_documentation.md](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash_documentation.md)): Detailed data source, structure, and function explanations for user support.
- **Python Packages** ([environment.yml](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/environment.yml)): Environment file with essential packages for local deployment. 

## Getting started 

1. Clone the repository:

```bash
    git clone https://github.com/ECMWFCode4Earth/TropiDash
```

2. Create the conda environemnt with the required dependecies to run the dasboard:

```bash
    conda env create -f environment.yml
```

3. Launch the dashboard by running the [TropiDash_backbone.ipynb](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash/TropiDash_backbone.ipynb)

## Getting Involved

Contributions to TropiDash are welcome! Whether it's enhancing functionality, improving documentation, or reporting issues, your input helps us make TropiDash better for everyone.

## License

This project is licensed under the Apache-2.0 License - see the LICENSE file for details.

