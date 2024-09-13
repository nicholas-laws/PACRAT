# Created by: Nicholas Laws
# Date: 2024

# Import relevant functions.
from comb_pathways import *
from plasma_pathways import *
from plasma_chemistry_objs import *
import os

# Get the directory of the currently running script.
dir = os.path.dirname(os.path.abspath(__file__))
# Define numerical constants.
R = 8.31432e3 
P0 = 101.325e3
FWHM = 20e-9
sigma = FWHM/(2*np.sqrt(2*np.log(2)))
# User-defined element of interest to be tracked.
element_of_interest = "O"
# Define combustion related data to parse.
combustion_mechanism = os.path.join(dir, "data/combustion/gri30.yaml") 
combustion_data = pd.read_csv(os.path.join(dir, "data/combustion/result_a00_b05_c00104_adj.csv"))
reverse_rate_coeffs = pd.read_csv(os.path.join(dir, "data/combustion/result_a00_b05_c00104_reverse_rr.csv")) 
combustion_temperature = combustion_data["T"].to_numpy()

# User-defined cutoff for combustion-related element flux.
combustion_cutoff = 1e-16

# Collect species and reaction relevant to the numerical combustion model.
combustion_species = get_species(combustion_mechanism)
combustion_reactions = get_reactions(combustion_mechanism)

# Generate the combustion-related reaction pathways and element flux.
combustion_flux_dict = get_flux_dict(combustion_reactions, element_of_interest, reverse_rate_coeffs, combustion_temperature, R, P0, combustion_data, combustion_cutoff)

# Define plasma related data to parse.
plasma_mechanism = os.path.join(dir, "data/plasma/kinet_CH4_PAC_V5_nrg.inp") 
plasma_data = pd.read_csv(os.path.join(dir, "data/plasma/result_a00_b05_c00104_PD.csv"))

# Read the .INP plasma mechanism of interest.
V5_plasma = read_inp(plasma_mechanism)

# Collect the species, neutrals, and reactions to be tracked from the .INP plasma mechanism of interest.
V5_species = find_species(V5_plasma)
V5_neutrals = find_neutrals(V5_species)
V5_reactions = find_reactions(V5_plasma)

# Create reaction objects for the .INP plasma mechanism of interest.
V5_reaction_objs = reaction_objects(V5_reactions, V5_species)

# Generate the plasma-related reaction pathways.
plasma_pathways_dict = collect_pathways(V5_reaction_objs, element_of_interest)

# User-defined cutoff for plasma-related element flux.
plasma_cutoff = 0

# User-defined boolean to determine if pulse period or interpulse is plotted (True = pulse period, False = interpulse period).
plot_pulse_period = True

# Generate the plasma-related element flux.
plasma_flux_dict = get_plasma_flux_dict(plasma_pathways_dict, V5_reaction_objs, element_of_interest, plasma_data, plasma_cutoff, V5_neutrals, plot_pulse_period, sigma)

# Combine the combustion and plasma-related pathways into a single dictionary.
PAC_pathways = coupled_pathways(combustion_flux_dict, plasma_flux_dict)

# Generate the network plot for reaction pathways.
coupled_network_plotter(PAC_pathways)