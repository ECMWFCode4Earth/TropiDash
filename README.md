# TropiDash: towards a comprehensive tropical cyclone hazard dashboard

TropiDash is a visualization tool designed to empower researchers, meteorologists, and enthusiasts with the ability to analyze and interpret tropical cyclone data effectively. Utilizing web technologies and data visualization libraries, TropiDash provides a comprehensive overview of tropical cyclones, their paths, intensities, and potential impacts. The dashboard is developed as a Jupyter Notebooks allowing also less experienced user to access it without problems. 
The project was developed as part of the ECMWF Code for Earth 2023 initiative. [Here](https://codeforearth.ecmwf.int/) you can fdin more details regarding Code for Earth. [Here](https://github.com/ECMWFCode4Earth/challenges_2023) you can find more details about the other projects of 2023.

Run the Dashboard directly on Binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ECMWFCode4Earth/TropiDash/HEAD?urlpath=voila%2Frender%2FTropiDash%2FTropiDash_backbone.ipynb)

## Repository structure

- **TropiDash folder**: Main notebook ([TropiDash_backbone.ipynb](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash/TropiDash_backbone.ipynb)), utility libraries, and data storage folders.
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

## Dashboard Preview

Get a closer look at TropiDash in action! These GIFs showcase some of the dashboard's features, demonstrating how users can interact with tropical cyclone data for insightful analysis.

Launching TropiDash in Jupyterlab. \\
![voila_launch](https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/40418464-992b-41cc-987c-0dc2a4d86a29)

Update the forecast after selecting the cyclone to visualize and the forecast date to consider. \\
![dashboard_launch](https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/7bd117ba-c799-4159-a08f-9f8da977b263)

Check the average forecast track and the information at a specific location.
![section1](https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/7fa6d1dd-fc98-4251-9fda-8021846b14ae)

Look at the total precipitation and the probability of 10m wind gusts over 25 m/s associated with the cyclone.
![section2](https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/f8f09d3b-e745-4753-8f6b-17fb494a0730)

Check the coastal hazard risk over the east coast of North America.
![section3](https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/4370a5c2-f69a-4a5b-98d9-7b6672105107)

Consult simultaneously the average forecast track and the coastal hazard risk of the region on the cyclone's path.
![section4](https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/468524fb-bd82-434f-a4a2-604b853d9392)

Check the precipitation forecast in a specific point along the average forecast track.
![section5](https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/f26d26d8-6a13-4b43-a3b4-7680c6f10083)

## Getting Involved

Contributions to TropiDash are welcome! Whether it's enhancing functionality, improving documentation, or reporting issues, your input helps us make TropiDash better for everyone.

## License

This project is licensed under the Apache-2.0 License - see the LICENSE file for details.

## Contact

For collaboration or any queries regarding this project, please contact:

- Filippo Dainelli: [filippo.dainelli@polimi.it]
- Paolo Colombo: [paolo1.colombo@polimi.it]
- Laura Paredes i Fortuny: [paredes_laufor@externos.gva.es]

