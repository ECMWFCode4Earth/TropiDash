# TropiDash: towards a comprehensive tropical cyclone hazard dashboard

TropiDash is a visualization tool designed to empower researchers, meteorologists, and enthusiasts with the ability to analyze and interpret tropical cyclone data effectively. Utilizing web technologies and data visualization libraries, TropiDash provides a comprehensive overview of tropical cyclones, their paths, intensities, and potential impacts. The dashboard is developed as a Jupyter Notebook allowing also less experienced user to access it without problems. 
The project was developed as part of the ECMWF Code for Earth 2023 initiative. You can find more details regarding the Code for Earth initiative at their [official website](https://codeforearth.ecmwf.int/). You can find more details about the other projects of 2023 at this [GitHub repository](https://github.com/ECMWFCode4Earth/challenges_2023) .

For the best experience with TropiDash, we recommend downloading and running the dashboard locally. This ensures optimal performance and access to the full range of features. Detailed instructions for setting up the dashboard on your local machine are provided in the [Getting Started section](#getting-started). Alternatively, for quick access and convenience, TropiDash is also available online via FIX THIS, allowing for immediate interaction without the need for local installation.

## Repository Structure

- **TropiDash folder:** Main notebook ([TropiDash_backbone.ipynb](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash/TropiDash_backbone.ipynb)), utility libraries, and data storage folders.
- **Tutorials Folder:** Step-by-step guides for each dashboard section and data downloading instructions.
- **Documentation ([TropiDash_documentation.md](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash_documentation.md)):** Detailed data source, structure, and function explanations for user support.
- **Python Packages ([environment.yml](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/environment.yml)):** Environment file with essential packages for local deployment. 

## Getting Started 

1. **Install Conda:** If you don't have Conda, download and install it from the official [Conda website](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
2. **Open Your Terminal:** Access your terminal or command prompt.
3. **Navigate to Your Desired Directory:** Use the `cd` command to navigate to the folder where you want to store TropiDash.
4. **Clone the repository:** Execute the clone command. You should have `git` installed, if you don't check the [Installation guide](https://git-scm.com/book/it/v2/Per-Iniziare-Installing-Git).
```bash
    git clone https://github.com/ECMWFCode4Earth/TropiDash
```
5. **Create and Activate the Conda Environment:** This ensures all dependecies are properly installed. Execute the 2 following commands:
```bash
    conda env create -f environment.yml
```
```bash
    conda activate tropidash
```
6. **Launch Jupyter Lab:** Navigate to the TropiDash directory in your terminal and execute:
```bash
    jupyter lab
```
7. **Open the main Notebook:** In Jupyter Lab, navigate to the folder containing [TropiDash_backbone.ipynb](https://github.com/ECMWFCode4Earth/TropiDash/blob/main/TropiDash/TropiDash_backbone.ipynb) and open it.
8. **Run the Dashboard:** Execute the dashboard by clicking on the Voilà button within Jupyter Lab. 
<img width="413" alt="Screenshot 2024-03-21 at 19 21 43" src="https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/93672154-9a4c-4c36-b067-56b2edebef22">

## Dashboard Preview

Get a closer look at TropiDash in action! These GIFs showcase some of the dashboard's features, demonstrating how users can interact with tropical cyclone data for insightful analysis.

*Launching TropiDash in Jupyterlab.*

![voila_launch](https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/40418464-992b-41cc-987c-0dc2a4d86a29)

*Update the forecast after selecting the cyclone to visualize and the forecast date to consider.*

![dashboard_launch](https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/7bd117ba-c799-4159-a08f-9f8da977b263)

*Check the average forecast track and the information at a specific location.*
![section1](https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/7fa6d1dd-fc98-4251-9fda-8021846b14ae)

*Look at the total precipitation and the probability of 10m wind gusts over 25 m/s associated with the cyclone.*
![section2](https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/f8f09d3b-e745-4753-8f6b-17fb494a0730)

*Check the coastal hazard risk over the east coast of North America.*
![section3](https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/4370a5c2-f69a-4a5b-98d9-7b6672105107)

*Consult simultaneously the average forecast track and the coastal hazard risk of the region on the cyclone's path.*
![section4](https://github.com/ECMWFCode4Earth/TropiDash/assets/54897571/468524fb-bd82-434f-a4a2-604b853d9392)

*Check the precipitation forecast in a specific point along the average forecast track.*
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

