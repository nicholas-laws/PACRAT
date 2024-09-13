# Created by: Nicholas Laws
# Date: 2024

## Imports
import pandas as pd
import time
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
from cantera_validation import *
import random
from IDT import *
import os

# Get the directory of the currently running script.
dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Define the filepath to the combustion mechanism. 
mechanism = os.path.join(dir, "data/combustion/gri30.yaml") 
# Initialize a Cantera gas object using the user-selected combustion mechanism.
gas = ct.Solution(mechanism)
# Define initial inlet gas conditions.
reactor_temperature = 1000.0  # K
reactor_pressure = ct.one_atm  # atm
inlet_concentrations = {"CH4": 1.0/10.52, "O2": 2.0/10.52, "N2": 7.52/10.52}
gas.TPX = reactor_temperature, reactor_pressure, inlet_concentrations
# Define reactor parameters.
residence_time = 0.1 # s
reactor_volume = 30.5 * (1e-2) ** 3  # m^3
# Define the maximum simulation runtime.
max_simulation_time = 1 # s
# Initialize the WSR.
fuel_air_mixture_tank = ct.Reservoir(gas)
exhaust = ct.Reservoir(gas)
stirred_reactor = ct.IdealGasConstPressureMoleReactor(gas, energy="off", volume=reactor_volume)
# Initialize a mass flow controller Cantera object where the mass flow rate is a function of the WSR mass and the user-selected residence time.
mass_flow_controller = ct.MassFlowController(
    upstream=fuel_air_mixture_tank,
    downstream=stirred_reactor,
    mdot=stirred_reactor.mass / residence_time,
)
# Initialize a pressure regular Cantera object.
pressure_regulator = ct.PressureController(
    upstream=stirred_reactor, downstream=exhaust, primary=mass_flow_controller
)
# Define the maximum time-step size.
max_step_size = 1.0e-2
# Define the WSR network.
reactor_network = ct.ReactorNet([stirred_reactor])
# Create a SolutionArray to store relevant data.
states = ct.SolutionArray(gas, extra=["t"])
# Initialize a variable that represents the reactions present in the user-selected combustion mechanism.
rxns = ct.Reaction.list_from_file(mechanism, gas)
# Collect the reaction equations for all reactions in the combustion mechanism.
rxn_strings = [rxn.equation for rxn in rxns]
# Initialize an empty dictionary where the key is the reaction equation and the value is an empty list.
reverse_rate_coeff = {rxn.equation: [] for rxn in rxns}
# Start the stopwatch to determine runtime.
tic = time.time()
# Set simulation start time to 0s.
t = 0
counter = 0
pressure = []
# Set the maximum time-step size to the user-defined value.
reactor_network.max_time_step = max_step_size
# Loop until the maximum time is reached, progressing by one time-step each iteration.
while t <= max_simulation_time:
    # Produce a time-step.
    t = reactor_network.step()
    # Extract the state of the reactor.
    states.append(stirred_reactor.thermo.state, t=t)
    temp_reverse_rate_coefficient = gas.reverse_rate_constants
    temp_pressure = gas.P
    pressure.append(temp_pressure)
    # Place the reverse reaction rate coefficient at time-step i in the overall dictionary each iteration, for each reaction tracked.
    for i, rxn in enumerate(rxns):
        reverse_rate_coeff[rxn_strings[i]].append(temp_reverse_rate_coefficient[i])
# Stop the stopwatch to determine runtime.
toc = time.time()
# Print out the simulation runtime.
print(f"Simulation runtime: {toc-tic:3.2f}s")
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
csv_file = os.path.join(dir, 'data/output/WSR.csv') 
reverse_rate_coeff_data.to_csv(csv_file, index=False)
states.save(os.path.join(dir, 'data/output/WSR.csv'), basis="mole", overwrite=True)
