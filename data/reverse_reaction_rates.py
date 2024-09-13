# Created by: Nicholas Laws
# Date: 2024

## Imports ##
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from comb_pathways import *
import cantera as ct
import math
import networkx as nx
import matplotlib.colors as mcolors
from fuzzywuzzy import fuzz
import re
import os

"""
A function to generate the time-dependent reverse reaction rate coefficients for the user-selected combustion mechanism.
"""
def reactor(mechanism, reactor_temperature, reactor_pressure, inlet_concentrations, input_temp, input_time, file_name):
    """
    Parameters:
    -----------
    mechanism : str
        A string representing the filename of the user-selected combustion mechanism.
    reactor_temperature : float
        A float representing the user-selected temperature of the WSR.
    reactor_pressure : float
        A float representing the user-selected pressure of the WSR.
    inlet_concentrations : Dict
        A dictionary where the key is the inlet species and the value is the concentration of these species.
    input_temp : np.array
        A numpy array representing the time-dependent temperature of the system.
    input_time : np.array
        A numpy array representing the time associated with the combustion data.
    file_name : str
        A string representing the output filename.
    """
    # Get the directory of the currently running script.
    dir = os.path.dirname(os.path.abspath(__file__))
    # Initialize a Cantera gas object using the user-selected combustion mechanism.
    gas = ct.Solution(mechanism)
    # Define initial inlet gas conditions.
    gas.TPX = reactor_temperature, reactor_pressure, inlet_concentrations
    # Initialize a variable that represents the reactions present in the user-selected combustion mechanism.
    rxns = ct.Reaction.list_from_file(mechanism, gas)
    # Collect the reaction equations for all reactions in the combustion mechanism.
    rxn_strings = [rxn.equation for rxn in rxns]
    # Initialize an empty dictionary where the key is the reaction equation and the value is an empty list.
    reverse_rate_coeff = {rxn.equation: [] for rxn in rxns}
    # Iterate through the time array.
    for j, t_j in enumerate(input_time):
        # Set the temperature and pressure at time t_j.
        gas.TP = input_temp[j], reactor_pressure
        # Collect the reverse rate coefficients at time t_j.
        temp_reverse_rate_coefficient = gas.reverse_rate_constants
        # Place the reverse reaction rate coefficient at time-step i in the overall dictionary each iteration, for each reaction tracked.
        for i, rxn in enumerate(rxns):
            reverse_rate_coeff[rxn_strings[i]].append(temp_reverse_rate_coefficient[i])
    # Filter pressure-dependent reactions out.
    arr_lens = np.array([len(reverse_rate_coeff[key]) for key in list(reverse_rate_coeff.keys())])
    reverse_rate_coeff_data = {}
    for key in list(reverse_rate_coeff.keys()):
        temp_len = len(reverse_rate_coeff[key])
        if temp_len > np.min(arr_lens):
            pass
        else:
            reverse_rate_coeff_data[key] = reverse_rate_coeff[key]
    # Save the reverse reaction coefficient data for all non-pressure-dependent combustion mechanism reactions. 
    reverse_rate_coeff_data = pd.DataFrame(reverse_rate_coeff_data)
    csv_file = os.path.join(dir, f'data/output/{file_name}.csv') 
    reverse_rate_coeff_data.to_csv(csv_file, index=False)

# Get the directory of the currently running script.
dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Define combustion related data to parse. 
mechanism = os.path.join(dir, "data/combustion/gri30.yaml") 
comb_data = os.path.join(dir, "data/combustion/result_a00_b05_c00104_adj.csv") 
# Define initial inlet gas conditions.
reactor_temperature = comb_data["T"].to_numpy()[0]  # K
reactor_pressure = ct.one_atm  # atm
inlet_concentrations = {"CH4": comb_data["X_ch4"].to_numpy()[0], "O2": comb_data["X_o2"].to_numpy()[0], "N2": comb_data["X_n2"].to_numpy()[0]}
# Define the output filename.
file_name = "result_a00_b05_c00104_reverse_rr"
reactor(mechanism, reactor_temperature, reactor_pressure, inlet_concentrations, comb_data["T"].to_numpy(), comb_data["t"].to_numpy(), file_name)