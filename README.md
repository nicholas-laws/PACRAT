# PACRAT — <i><u>P</u>lasma-<u>A</u>ssisted <u>C</u>ombustion <u>R</u>eaction Pathways <u>A</u>nalysis <u>T</u>ool</i>

**PACRAT**, <i><u>P</u>lasma-<u>A</u>ssisted <u>C</u>ombustion <u>R</u>eaction Pathways <u>A</u>nalysis <u>T</u>ool</i>, is an automated reaction pathways analysis tool for the numerical modeling of plasma-assisted combustion kinetics. By reading any arbitrary plasma and/or combustion mechanism, PACRAT can identify the relevant species pathways and compute the total element flux for each pathway given an element to track, time-dependent reaction rate, temperature, electron temperature, and number densities data. These pathways are visualized in a network diagram, where each pathway is associated a gradient color corresponding to the element flux.

---

## Table of Contents
1. [Installation](#installation)
2. [Usage](#usage)
3. [License](#license)
4. [Contact Information](#contact-information) 

---

## Installation

### Prerequisites

Before you begin, ensure you have the following installed on your system:

- **Python 3.8 or Higher**: PACRAT is built using Python. You can download the latest version from [python.org](https://www.python.org/downloads/).
- **Git**: Required to clone the repository. Download and installation instructions are available at [git-scm.com](https://git-scm.com/downloads).

### Setup

In order to install PACRAT on your machine, please follow these instructions:

#### 1. Clone the PACRAT Repository

Open your terminal or command prompt and run the following command to clone the repository:

```bash
git clone https://github.mit.edu/MIT-APG/PAC_Reaction_Pathways_Analysis.git](https://github.com/nicholas-laws/PACRAT.git
```

#### 2. Navigate into the PACRACT directory
```bash
cd PAC_Reaction_Pathways_Analysis
```

#### 3. Install dependencies
Install the required dependencies using the provided requirements.txt file.
```bash
pip install -r requirements.txt
```

#### 4. Verify installation
To ensure all dependencies are installed correctly, run the following command:
```bash
python -c "import numpy; import pandas; import cantera; import matplotlib; import networkx; import scipy; import fuzzywuzzy; import netgraph; print('All dependencies are installed correctly.')"
```

#### 5. Download the `data` folder
By accessing the APG Dropbox, you will find a folder labeled `Nick Thesis - March 13 2024`. Please download the `data.zip` file and move this file into the `PAC_Reaction_Pathways_Analysis` directory. In the `PAC_Reaction_Pathways_Analysis`, run the following command:
```bash
unzip data.zip; rm -i *.zip
```

### Troubleshooting
* Missing Packages: If you encounter an ImportError, ensure that the package is listed in requirements.txt and installed correctly.

* Version Conflicts: If a package version causes compatibility issues, consider adjusting the version in requirements.txt and reinstalling.

## Usage

### Overview

**PACRAT** is a Python-based tool designed for comprehensive reaction pathways analysis in plasma-assisted combustion kinetics. It allows researchers to model, analyze, and visualize complex chemical reaction networks to better understand the underlying combustion processes.

### Running the Main Analysis Script

To perform a reaction pathway analysis, you'll primarily interact with the `main.py` script located in the project's `src` directory. This script orchestrates the reaction pathway analysis by loading relevant plasma and combustion mechanisms, time-dependent species data, processing reaction pathways, computing element flux, and visualizing the results.

#### Basic Command

```bash
cd src
python main.py
```

In the `main.py` script, you will have access these parametric inputs:
* `element_of_interest`: The element of interest you wish to track for your reaction pathways analysis.
* `combustion_mechanism`: The combustion kinetic mechanism you wish to consider in Cantera `.yaml` format.
* `combustion_data`: The time-dependent temperature, species mass fractions, and reaction rates for the combustion kinetic mechanism.
* `reverse_rate_coeffs`: The time-dependent reverse reaction rate coefficients for all reactions that are not pressure-dependent.
* `combustion_cutoff`: The element flux cutoff for all combustion kinetics reaction pathways.
* `plasma_mechanism`: The plasma kinetic mechanism you wish to consider in `.INP` format.
* `plasma_data`: The time-dependent temperature, electron temperature, electron number density, species number densities, and reaction rate coefficients for the plasma kinetic mechanism. 
* `plasma_cutoff`: The element flux cutoff for all plasma kinetics reactionn pathways

This repository is written in the form of a tutorial, having a variety of associated combustion and plasma kinetics data. In `main.py`, associated combustion and plasma kinetics data from a stoichiometric PAC numerical experiment run with an initial temperature of 1000 K. This associated data is denoted by the file naming convention of `result_a00_b05_c00040_X.csv`, where `X` can be `PD`, representing the time-dependent plasma data, `adj`, representing the time adjusted combustion kinetics on the ZDPlasKin timescale, and `reverse_rr`, representing the reverse reaction rate coefficients for all non-pressure dependent reactions in the combustion kinetics mechanism.

### Generating data for analysis

After running a numerical experiment with the `ZD_Plasma_Reactor` and extracting combustion and plasma data, time-dependent reverse reaction rate coefficients for combustion kinetics must be computed. Additionally, you must adjust the combustion data to be on the ZDPlasKin timescale.

#### Computing the reverse reaction rate coefficients for combustion kinetics

In the `data` folder, you will find `reverse_reaction_rates.py`. This script initializes a well-stirred reactor in Cantera in order to compute and collect the reverse reaction rate coefficients into a specified output `.CSV` from a pre-existing numerical combustion experiment. In the `reverse_reaction_rates.py` script, you will have access to the `reactor` function. This function allows you to manipulate the well-stirred reactor through the user-selected parameters:

* `mechanism`: The combustion kinetics mechanism of interest.
* `reactor_temperature`: The initial temperature of combustion kinetics.
* `reactor_pressure`: The constant pressure of the well-stirred reactor.
* `inlet_concentrations`: The inlent concentrations of the main combustion species.
* `input_temp`: The time-dependent temperature data associated with the numerical combustion experiment.
* `input_time`: The time data associated with the numerical combustion experiment.
* `file_name`: The output file name of choice.

#### Adjusting the timescale of the combustion kinetics
In the `data` folder, you will find `adjust_timescale.py`. This script iterates through a numerical combustion experiment, checking if the time value tracked in the combustion data is also tracked in the plasma data. The user will have access to these parameters:
* `plasma_data`: The overall time-dependent number densities and reaction rate coefficients for all plasma kinetics data.
* `comb_data`: The overall time-dependent mass fraction data for all combustion reactions in the numerical model.
* `file_name`: The user-selected output filename for the adjusted combustion data.

## License

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

This project is licensed under the [MIT License](LICENSE). See the [LICENSE](LICENSE) file for details.

## Contact Information

If you have any questions, suggestions, or need support with **PACRAT**, feel free to reach out to me directly:

### **Author**

**Name:** Nicholas Laws

**Role:** Undergraduate Research Assistant
 
**Affiliation:** MIT Aerospace Plasma Group

### **Email**

For direct inquiries, please email me at:  
📧 [nrl49@cornell.edu](mailto:nrl49@cornell.edu)

### **LinkedIn**

Professional networking and collaborations:  
🔗 [https://www.linkedin.com/in/nicholasrlaws/](https://www.linkedin.com/in/nicholasrlaws/)

### **PACRAT Repository**

For issues, feature requests, or contributions, please use the GitHub Issues page:  
🐞 [https://github.mit.edu/MIT-APG/PAC_Reaction_Pathways_Analysis/issues](https://github.mit.edu/MIT-APG/PAC_Reaction_Pathways_Analysis/issues)

---
